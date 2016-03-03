#! /usr/bin/python
# -*- coding: utf-8 -*-
"""A program to calculate aggregation number distribution over time, of
both surfactants and solvent/counterions/substrate"""
# Import necessary utilities
from __future__ import print_function
import argparse
import multiprocessing as mp
import datetime
import sys
import time
import cProfile
from re import search
from math import sqrt
from decimal import getcontext, Decimal
from copy import copy
from Queue import Empty
# Global variables to speed up individual functions
coordinatePattern = r'(?<!bo)x \([0-9]{1,}x[0-9]{1,}\):'
boxPattern = r' {0,}box \([0-9]{1,}x[0-9]{1,}\):'
atomPattern = r'x\[ {,}[0-9]{1,}\]\=\{'
framePattern = r'[a-z, A-Z, 0-9, _]{1,}\.[a-z, A-Z, 0-9]{1,} frame [0-9]*:'
allEntryPattern = r'Time [0-9]{1,}\\n'
coordinateLength = 0
frameLength = 0

# Find the average number per list in a list of lists
# E.g if there was a list [[1, 5, 6], [3, 4, 5]]
# the function would return [4, 4] to the file
def average_value_list(filename, function):
    avs = []
    av_sum = 0
    num = 0.0
    for state in function:
        for number in state:
            av_sum += number
            num += 1
        if num != 0:
            av = av_sum/num
            avs.append(av)
        else:
            avs.append(0)
    return avs

# It's easier and less complicated (for now) to just do all this in the display
# function using find_aggregation_numbers()
# Prints the average aggregation number to the specified file
def average_micelle_print(args, output):
    listData = []
    for frame in output:
        # Sorts dict by aggregation number, with frame number being last
        # Removes the frame number to allow max sorting
        ags = sorted(frame.items())
        ags = ags[0:len(ags)-1]
        # Find sum
        av_sum = 0
        num = 0.0
        for pair in ags:
            av_sum += pair[0]*pair[1]
            num += pair[1]
        av = av_sum/num
        # The index is the frame number
        listData.append(av)
    function_print(args.output, listData)
    return 0

# Expects a list where the index is the independent variable
# and the value is the dependent variable.
def function_print(filename, function):
    with open(filename, 'w') as f:
        for value in function:
            print(str(function.index(value))+', '+str(value), file=f)
    return 0

# Account for periodic boundary conditions
def pbc(positionA, positionB, boxSize, cutoff):
    newPositionA = []
    newPositionB = []
    bigger = [0, 0, 0]
    if positionA[0] >= positionB[0]:
        xdist = positionA[0]-positionB[0]
        bigger[0] = 0
    else:
        xdist = positionB[0]-positionA[0]
        bigger[0] = 1
    if positionA[1] >= positionB[1]:
        ydist = positionA[1]-positionB[1]
        bigger[1] = 0
    else:
        ydist = positionB[1]-positionA[1]
        bigger[1] = 1
    if positionA[2] >= positionB[2]:
        zdist = positionA[2]-positionB[2]
        bigger[2] = 0
    else:
        zdist = positionB[2]-positionA[2]
        bigger[2] = 1
    dists = [xdist, ydist, zdist]
    bc = boxSize-cutoff
    for dist in dists:
        index = dists.index(dist)
        if abs(dist) >= bc:
            if bigger[index] == 0:
                newPositionA.append(positionA[index]-boxSize)
                newPositionB.append(positionB[index])
            else:
                newPositionB.append(positionB[index]-boxSize)
                newPositionA.append(positionA[index])
        else:
            newPositionA.append(positionA[index])
            newPositionB.append(positionB[index])
    return newPositionA, newPositionB

# Find the straight-line distance between two points
def distance(positionA, positionB, boxSize, cutoff):
    positionA, positionB = pbc(positionA, positionB, boxSize, cutoff)
    # For positions [x, y, z] in nanometers
    xdist = abs(positionA[0])-abs(positionB[0])
    ydist = abs(positionA[1])-abs(positionB[1])
    zdist = abs(positionA[2])-abs(positionB[2])
    dist = sqrt(xdist**2 + ydist**2 + zdist**2)
    return dist

# Get a position in Angstroms from a line
def extract_pos(line, xtc=False):
    # if a tpr file: "   x[  388]={ 7.55358e+00,  4.29939e+00,  4.35216e+00}"
    # if an xtc, move everything to the right by 3 chars
    #      x[11550]={ 1.13700e+00,  4.32100e+00,  3.07000e-01}
    # In the above case, where the atom number is in the ten thousands, subtract
    # the line length from 58 to account for it.
    expand = abs(len(line)-58)
    if xtc:
        xstring = line[17+expand:28+expand]
        ystring = line[31+expand:42+expand]
        zstring = line[45+expand:56+expand]
    else:
        xstring = line[14+expand:25+expand]
        ystring = line[28+expand:39+expand]
        zstring = line[42+expand:53+expand]
    x = float(xstring[0:7])*(10**float(xstring[8:11]))*10 # in Angstroms
    y = float(ystring[0:7])*(10**float(ystring[8:11]))*10 # in Angstroms
    z = float(zstring[0:7])*(10**float(zstring[8:11]))*10 # in Angstroms
    return (x, y, z)

# Find the next frame (line number of first coordinate, byte of first coord)
# from a file
# We have to start at the beginning of the file every time.
def find_next_instance(fileName, pattern, start=0):
    with open(fileName, 'r') as f:
        enum = enumerate(f)
        for lineNumber, line in enum:
            if lineNumber < start:
                continue
            check = search(pattern, line)
            if check:
                return lineNumber
    return None

# Given a certain amount of lines and a starting position, read that many
# lines from a file and return them as a list
def get_lines(fileName, length, start=0):
    lines = []
    with open(fileName, 'r') as f:
        f.seek(0)
        enum = enumerate(f)
        for lineNumber, line in enum:
            if lineNumber < start:
                continue
            if lineNumber == start+length:
                return lines
            if lineNumber >= start:
                lines.append(line)
    return lines

def find_next_frame_list(lines, pattern, start=0):
    enum = lines[start:]
    for line in enum:
        check = search(pattern, line)
        if check:
            return enum.index(line)
    return None

# Returns the fist distance found below cutoff
def molecule_cutoff(moleculeA, moleculeB, cutoff, boxSize):
    for atomA in moleculeA['atoms']:
        for atomB in moleculeB['atoms']:
            dist = distance(atomA, atomB, boxSize, cutoff)
            if dist <= cutoff:
                return dist
    return None

# Factory function for constructing the parsing object
def get_parser():
    # Argument parsing object
    parser = argparse.ArgumentParser(description='Calculate number of \
        monomers in micelle.', \
        formatter_class=argparse.RawDescriptionHelpFormatter)

    # Get name of .trr file
    parser.add_argument('--file', '-f', required=False, action='store', \
                        default=None)
    # Get name of .tpr file
    parser.add_argument('--tpr', '-s', required=False, action='store', \
                        default=None, help='A text .tpr file. Because \
                        trajectory files and .tpr files have different \
                        formats, this program will crash if you specify a \
                        trajectory file here.')
    # Get names of molecules in file
    parser.add_argument('--names', '-nms', required=True, action='store', \
                        help='The names of the types of molecules, in order of\
                        appearance in the .xtc or .tpr file. Call the micelle-\
                        former MCL. If you want to skip something call it \
                        SKIP. You can call multiple things MCL and SKIP.', \
                        type=str, nargs='+', default=['MCL'])
    # Get number of amphiphilic molecules in solution
    parser.add_argument('--nmols', '-nm', required=True, action='store', \
                        help='The number of copies of molecules in \
                        solution. Use the same order as the --names flags.', \
                        type=int, nargs='+')
    # Get number of atoms in amphiphilic molecule
    parser.add_argument('--nums', '-n', required=True, action='store', \
                        help='The number of atoms in the molecules. Use the \
                        same order as the --names flags.', type=int, nargs='+')
    # Get distance cutoff
    parser.add_argument('--cutoffs', '-c', required=True, action='store', \
                        default=[4.0], help='The radius of the circle around \
                        an atom that determines whether another atom is \
                        close to it. Use same order as the --names flags.', \
                        type=float, nargs='+')
    # Skip this many frames in a trajectory file
    #parser.add_argument('--frameskip', '-df', required=False, \
    #                    action='store', help='The number of frames to skip \
    #                    between every analysis step on a trajectory. \
    #                    Default 0 to analyze every frame. You can use trjconv \
    #                    and gmxdump to output a trajectory with frames left \
    #                    out for speed, or you can do the same with \
    #                    frameskip.', type=int, default=0)
    # Get name of output file
    parser.add_argument('--output', '-o', required=False, action='store', \
                        help='The name of the output file for the micelle \
                        average aggregation numbers.', type=str, \
                        default='aggregation.csv')
    # Get distance cutoff
    #parser.add_argument('--boxsize', '-b', required=True, action='store', \
    #                    default=50.0, help='The length of one side of the \
    #                    simulation box in Angstroms.', type=float)
    # Get number of frames to scan in .xtc file
    parser.add_argument('--frames', '-fr', required=False, action='store', \
                        help='The number of frames in the trr or xtc file.', \
                        type=int, default=-1)
    # Get number of ps per frame
    parser.add_argument('--ps_per_frame', '-ppf', required=False, \
                        action='store', help='The number of picoseconds \
                        contained in one frame of data. Defaults to 1 to \
                        display frame number only.', type=int, default=10)
    # Get counterion file prefix
    parser.add_argument('--counterion', '-ci', required=False, action='store',\
                        type=str, default='counterion', \
                        help='The prefix used before all counterion \
                        aggregation files.')
    # Get whether MCL forms micelles or not
    parser.add_argument('--nomicelle', '-no', required=False, action='store', \
                        help='Whether the MCL must be in a micelle for counte-\
                        rion distribution to be calculated. To analyze the \
                        distribution around lone atoms too, use either 1 or \
                        True for this flag. Defaults to False.', \
                        type=bool, default=False)
    # Get the chunk size (the number of frames to be processed at one time
    parser.add_argument('--chunk', '-ch', required=False, action='store', \
                        help='The number of frames to process at one time in \
                        memory. Be aware that processing too many frames at \
                        one time will fill up your computer\'s RAM! ', \
                        type=int, default=10000)
    # Get the complete output file
    parser.add_argument('--toutput', '-to', required=False, action='store', \
                        help='The name for the log file. Includes the complete\
                        output of all data.', type=str, default='allout.txt')
    # Get the number of threads
    parser.add_argument('--threads', '-t', required=False, action='store', \
                        help='The number of parallel processes to create to \
                        scan the .xtc file.', type=int, default=4)
    # Get the number of threads
    parser.add_argument('--size', '-sz', required=False, action='store', \
                        help='The number of monomers required for an aggregate \
                        to be counted as a micelle. Default = 2.', type=int, default=2)
    parser.epilog = '''
Here are some examples of proper usage, ranging from simple to complex.

$ aggregate -f md.trr --names MCL NA --nmols 5 5 --nums 15 1 -o
micelle_size.txt -ci Sodium -c 4.0 4.0

    This command will analyze a trajectory file, assemble a list of micelles,
analyze the sodium ion distribution around the micelles, and print everything
out in two files: micelle_size.txt and SodiumNA.txt. micelle_size.txt will
contain the average aggregation number for every time t in picoseconds.
SodiumNA.txt will contain the average number of sodium ions within 4 Angstroms
of a micelle for every time t in picoseconds.

Here is a slightly more complicated example. We want to see the average amount
of water around a certain type of molecule every two picoseconds. To do this we
call the molecule in question MCL as usual but also include the flag
--nomicelle True.

$ aggregate -f md.trr --names MCL SOL NA --nmols 15 2000 15 --nums 32 3 1
-o micelle_size.txt -ci Counterion -c 4.0 4.0 5.0 -ppf 2 -no True

    This command will output the usual micelle_size.txt, but CounterionNA.txt
will be the average for every molecule MCL, not just the ones in micelles.
Additionally, instead of every picosecond this command will correctly analyze
data output every 2 picoseconds in md.trr.

Here is a fairly complex example involving many types of counterions,
substrates, and solvent:

$ aggregate -f md.trr --names MCL SOL NA OBN KET --nmols 50 34584 85 2 3
--nums 77 3 1 30 5 -ci Counterion -o micelle_size.txt -c 4.0 4.5 5.0 3.9 4.0
-ppf 100 --nomicelle True --chunk 10000 -fr 100000 -t 12 -fs 2

    This run of the program will create 12 threads to deal with
two substrates, OBN and KET, a counterion, NA, solvent, SOL, and a
micelle-forming molecule, MCL. The execution of the program will be broken up
into chunks of 10000 frames each, or 1000000 ps each, for a total of
10,000,000 ps or 10 microseconds of simulation time. Additionally, the program
will only analyze every other frame, saving time.

In order to begin using this program you will need to create an ASCII version
of your Gromacs trajectory or topology file using trjconv and dump.
Typically the following will work:
$ gmx trjconv -f md.xtc -o aggregate.xtc -pbc atom -ur rect
$ gmx dump -f aggregate.xtc > aggregate.txt
$ aggregate ...
'''
    return parser

# Loop through every molecule, use the cutoff radius to calculate
# micelle distribution
def calculate_micelles(solution, cutoffs, boxSize, names, nomicelle):
    for moleculeA in solution['molecules']:
        if moleculeA['name'] != 'MCL':
            continue
        for moleculeB in solution['molecules']:
            if moleculeB['name'] != 'MCL':
                continue
            if moleculeA == moleculeB:
                continue
            index = names.index('MCL')
            dist = molecule_cutoff(moleculeA, moleculeB, cutoffs[index], \
                                   boxSize)
            if dist != None:
                add_new_micelle(solution, moleculeA, moleculeB)
            continue
    # If we need to analyze every MCL molecule, just add every one to its
    # own micelle if it is not in one.
    if nomicelle == True or nomicelle == 1:
        for molecule in solution['molecules']:
            if molecule['name'] != 'MCL':
                continue
            for micelle in solution['micelles']:
                if molecule in micelle['molecules']:
                    break
            else:
                solution['micelles'].append({'micelleIndex':\
                    get_micelle_index(), 'molecules':[molecule], \
                    'near':{'all':[]}})
                continue
    return


# Loop through every molecule and use the cutoff radius to determine
# counterion distribution
def calculate_near(solution, cutoffs, boxSize, names):
    for moleculeA in solution['molecules']:
        if moleculeA['name'] == 'MCL':
            # Haven't found an ion/substrate yet, so move on.
            continue
        for moleculeB in solution['molecules']:
            if moleculeB['name'] != 'MCL':
                # Only check against surfactants in micelles
                continue
            theMicelle = None
            for micelle in solution['micelles']:
                if moleculeB in micelle['molecules']:
                    theMicelle = micelle
                    break
            if theMicelle != None:
                index = names.index(moleculeA['name'])
                dist = molecule_cutoff(moleculeA, moleculeB, \
                                       cutoffs[index], boxSize)
                if dist != None:
                    if moleculeA not in theMicelle['near']['all']:
                        theMicelle['near']['all'].append(moleculeA)
                        if moleculeA['name'] in theMicelle['near'].keys():
                            theMicelle['near'][moleculeA['name']] = \
                                theMicelle['near'][moleculeA['name']] + 1
                        else:
                            theMicelle['near'][moleculeA['name']] = 1
                    continue
    return

def add_dicts(dictA, dictB):
    newDict = {}
    for key in dictA.keys():
        newDict[key] = dictA[key]
    for key in dictB.keys():
        if key in newDict.keys():
            newDict[key] = newDict[key] + dictB[key]
        else:
            newDict[key] = dictB[key]
    return newDict

def add_new_micelle(solution, moleculeA, moleculeB):
    micelleA = None
    micelleB = None
    for micelle in solution['micelles']:
        if moleculeA in micelle['molecules']:
            micelleA = micelle
            break
    for micelle in solution['micelles']:
        if moleculeB in micelle['molecules']:
            micelleB = micelle
            break
    if micelleA != None and micelleA == micelleB:
        return
    if (micelleA, micelleB) == (None, None):
        solution['micelles'].append({'micelleIndex':get_micelle_index(), \
                                    'molecules':[moleculeA, moleculeB], \
                                    'near':{'all':[]}})
    elif micelleA != None and micelleB != None:
        solution['micelles'].remove(micelleA)
        solution['micelles'].remove(micelleB)
        new_near = add_dicts(micelleA['near'], micelleB['near'])
        solution['micelles'].append({'micelleIndex':get_micelle_index(), \
                                    'molecules':micelleA['molecules']+\
                                    micelleB['molecules'], \
                                    'near':new_near})
    elif micelleA != None:
        micelleA['molecules'].append(moleculeB)
    elif micelleB != None:
        micelleB['molecules'].append(moleculeA)
    return

def get_solution_index():
    return 's' + str(datetime.datetime.now().microsecond)

def get_micelle_index():
    return datetime.datetime.now().microsecond

def printall(solution, names, header=''):
    text = str(header) + '\n'
    for micelle in solution['micelles']:
        out = 'Micelle #' + str(micelle['micelleIndex']) + ': '
        amount = out + str(len(micelle['molecules'])) + '\n'
        text = text + amount
    out = ''
    for counterion in names:
        if counterion in ['MCL', 'SKIP']:
            continue
        out = out + 'Counterion-'+counterion+':\n'
        for micelle in solution['micelles']:
            if counterion in micelle['near'].keys():
                out = out + '    ' + str(micelle['near'][counterion]) + \
                    ' near Micelle #' + str(micelle['micelleIndex']) + '\n'
    final = text + out
    atoms = 0
    for molecule in solution['molecules']:
        atoms = atoms + len(molecule['atoms'])
    atoms = 'Number of atoms: ' + str(atoms) + '\n'
    molecules = 'Number of molecules: '+str(len(solution['molecules'])) + '\n'
    return final + molecules + atoms

def populate_frame(args, lines, xtc=False, positions=None):
    solution = {'micelles':[], 'molecules':[], 'near':{'all':[]}, \
        'solutionIndex':get_solution_index()}
    linesProcessed = 0
    # The line numbers of switching to a different molecule
    if positions == None:
        positions = get_positions(args)
    atoms = []
    frameIndex = 0
    for line in lines:
        if linesProcessed == positions[frameIndex]:
            frameIndex = frameIndex + 1
        if args.names[frameIndex] == 'SKIP':
            linesProcessed = linesProcessed + 1
            continue
        #print(line)
        pos = extract_pos(line, xtc)
        atoms.append(pos)
        linesProcessed = linesProcessed + 1
        if len(atoms) == args.nums[frameIndex]:
            molecule = {'name':args.names[frameIndex], 'atoms':atoms}
            solution['molecules'].append(molecule)
            atoms = []
    #number = {}
    #for molecule in solution['molecules']:
    #    if molecule['name'] in number.keys():
    #        number[molecule['name']] = number[molecule['name']] + 1
    #    else:
    #        number[molecule['name']] = 1
    return solution

# Scan a tpr file
def scan_tpr(args):
    global framePattern
    startLine = find_next_instance(args.tpr, framePattern)
    linesToProcess = find_coordinate_length(args)
    lines = []
    get_lines(args.file, linesToProcess, start=startLine)
    solution = populate_frame(args, lines, xtc=False)
    calculate_micelles(solution, args.cutoff, args.boxsize, args.names, \
                       args.nomicelle)
    calculate_near(solution, args.cutoff, args.boxsize, args.names)
    print(printall(solution, '.tpr file statistics'))
    return

# Returns a list of dicts: [{'frame':frameno, aggregation_number:number of
# micelles in frame with that number}, ]
def find_aggregation_numbers_list(frames, nomicelle, cutoff):
    aggregation = {}
    c = 0
    for frame in frames:
        data = {'MCL':True, 'time':frame['time']}
        for micelle in frame['micelles']:
            ag_number = len(micelle['molecules'])
            if ag_number == 1 and nomicelle == False:
                continue
            if ag_number < cutoff and nomicelle == False:
                continue
            if ag_number in data:
                data[ag_number] = data[ag_number] + 1
            else:
                data[ag_number] = 1
        aggregation[frame['time']] = data
        c = c + 1
    return aggregation

# Returns a list of dicts: [{'frame':frameno, aggregation_number:number of
# micelles in frame with that number}, ]
def find_aggregation_numbers(frame, nomicelle):
    data = {'time':frame['time'], 'MCL':True}
    for micelle in frame['micelles']:
        ag_number = len(micelle['molecules'])
        if ag_number == 1 and nomicelle == False:
            continue
        if ag_number in data:
            data[ag_number] = data[ag_number] + 1
        else:
            data[ag_number] = 1
    return data

# This function only applies to trajectory files; it will fail for
# tpr files
def find_frame_count(args, output=True):
    global framePattern, frameLength
    suggested = args.frames
    counter = 0
    with open(args.file, 'r') as f:
        f.seek(0)
        while True:
            try:
                f.next()
            except StopIteration:
                break
            counter = counter + 1
            if suggested != -1 and counter/frameLength == suggested:
                if output: print('\n', end='')
                return suggested
            if counter%1000 == 0 and output:
                print('\rframeCounter: '+str(counter), file=sys.stdout, end='')
    if output: print('\n', end='')
    return counter/frameLength

# Returns how many atoms per frame we should process, neglecting the
# fluff of box sizes, etc. included in the frame data by Gromacs
def find_coordinate_length(args):
    length = 0
    for index, name in enumerate(args.names):
        length = length + (args.nums[index]*args.nmols[index])
    return length


def find_frame_length(args, pattern):
    # Because every frame is a set distance apart, just subtract
    frameOne = find_next_instance(args.file, pattern, start=0)
    frameTwo = find_next_instance(args.file, pattern, start=1)
    return frameTwo-frameOne

# Get the number of atoms per type of molecule by name
def get_positions(args):
    positions = []
    linesForName = []
    for index, name in enumerate(args.names):
        nums, nmols = args.nums[index], args.nmols[index]
        linesForName.append(nums*nmols)
        positions.append(sum(linesForName))
    return positions

def populate(args, frameCount, extCounter):
    global coordinatePattern, framePattern, coordinateLength, frameLength, \
           atomPattern
    myCounter = 0
    frames = []
    positions = get_positions(args)
    # We may or may not be starting at 0, so if not, multiply by the number
    # of frameswe have already completed, adding from the first frame
    # to neglect the fluff Gromacs adds to the beginning of a trr file
    # Also we do not need to f.seek(0) because the find_next_instance
    # and get_lines will do that for us.
    startLine = 0
    if extCounter != 0:
        startLine = find_next_instance(args.file, framePattern, 0)
        # Again, we may not be the first frame. Multiply frameLength
        # by previous number of frames to find our chunk's position
        startLine = extCounter*frameLength
    # If counter started at 1 like frameCount, this would be <=
    frameLine = 0
    while myCounter < frameCount:
        if myCounter == 0:
            # The first coordinate is 1 line after the coordinate header
            #
            frameLine = find_next_instance(args.file, coordinatePattern, \
                                           startLine)
        else:
            # Don't waste time looking for frames, just use the standard
            # length of a frame (but make sure to subtract the one you added
            # earlier!)
            frameLine = startLine + frameLength
            #frameLine = find_next_instance(args.file, coordinatePattern, \
            #                               startLine)
        if frameLine == None:
            print('Forcing break')
            # There are no more frames left.
            break
        if myCounter == 0:
            frameLine = frameLine + 1
        # Start looking from current position in file to skip previously-
        # calculated frames
        startLine = frameLine
        # Frame-skipping will eventually be operational
        #if args.frameskip != 0:
        #    if myCounter%args.frameskip == 0:
        #        myCounter = myCounter + 1
        #        continue
        frameLines = get_lines(args.file, coordinateLength+2, start=frameLine-2)
        solution = populate_frame(args, frameLines[2:], xtc=True, \
                                  positions=positions)
        solution['time'] = (extCounter+myCounter)*args.ps_per_frame
        solution['boxSize'] = extract_pos(frameLines[0], xtc=True)[2]
        frames.append(solution)
        myCounter = myCounter + 1
    return frames

def analyze(args, frames):
    calculateCounter = 0
    for solution in frames:
        calculate_micelles(solution, args.cutoffs, solution['boxSize'], \
                           args.names, args.nomicelle)
        calculate_near(solution, args.cutoffs, solution['boxSize'], args.names)
        calculateCounter = calculateCounter + 1
    # Return data for processing by other functions
    return frames

def scan_chunk(args, frameCount, extCounter):
    chunk = populate(args, frameCount, extCounter)
    data = analyze(args, chunk)
    return data

def scan_chunks(args, frameCount, out, extCounter):
    if frameCount%args.chunk == 0:
        pass
    else:
        # If not this must be the last chunk, smaller than the others.
        #data = cProfile.runctx('scan_chunk(args, frameCount, extCounter)', \
        #                    globals(), locals())
        data = scan_chunk(args, frameCount, extCounter)
        out.put(data)
        return
    for i in range(0, frameCount/args.chunk):
        #data = cProfile.runctx('scan_chunk(args, args.chunk, (extCounter+(i*args.chunk)))', \
        #                    globals(), locals())
        data = scan_chunk(args, args.chunk, (extCounter+(i*args.chunk)))
        out.put(data)
    return

def average(data):
    mine = copy(data)
    for key in mine.keys():
        if key in ['MCL', 'time']:
            mine.pop(key)
    av_sum = 0
    num = 0.0
    items = []
    for key in mine.keys():
        items.append((key, mine[key]))
    for ag_number, amount in items:
        # ag_number is aggregation number of a particular micelle
        # amount is the number of micelles with that aggregation
        # number
        av_sum = av_sum + ag_number*amount
        num = num + amount
    if num == 0:
        return 0
    else:
        return av_sum/num
    return 0

def find_mmap_size(theMap, start):
    after = 0
    theMap.seek(start)
    while theMap.tell() != theMap.size()-1:
        theMap.read_byte()
        after = after + 1
    return after

def start_processes(arguments):
    frameCount = find_frame_count(arguments, output=True)
    processes = []
    outQueue = mp.Queue()
    counter = 0
    numberChunks = 0
    # This should always work because the parser will return an int number
    # of frames per chunk
    if arguments.chunk < 1:
        return
    if frameCount%arguments.chunk == 0:
        # Number of frames is exact multiple of arguments.chunk
        # or equals arguments.chunk
        # Use integer division
        numberChunks = frameCount/arguments.chunk
    else:
        # number of frames is not exact multiple of arguments.chunk
        # or is less than arguments.chunk
        if frameCount < arguments.chunk:
            numberChunks = 0
        else:
            # Use floor division. and add 1 for the extra chunk
            numberChunks = frameCount//arguments.chunk
    if arguments.threads == 0:
        return
    if numberChunks == arguments.threads:
        for i in range(0, arguments.threads):
            p = mp.Process(target=scan_chunks, args=(arguments, \
                           arguments.chunk, outQueue, counter))
            p.order = i
            processes.append(p)
            p.start()
            counter = counter + arguments.chunk
    elif numberChunks > arguments.threads and \
         numberChunks%arguments.threads == 0:
        chunksPerProcess = numberChunks/arguments.threads
        for j in range(0, arguments.threads):
            p = mp.Process(target=scan_chunks, args=(arguments, \
                           arguments.chunk*chunksPerProcess, outQueue, counter))
            p.order = j
            processes.append(p)
            p.start()
            counter = counter + arguments.chunk*chunksPerProcess
    elif numberChunks > arguments.threads and \
         numberChunks%arguments.threads !=0:
        remainingChunks = numberChunks%arguments.threads
        chunksPerProcess = numberChunks//arguments.threads
        for k in range(0, arguments.threads):
            if remainingChunks == 0:
                p = mp.Process(target=scan_chunks, args=(arguments, \
                               arguments.chunk*chunksPerProcess, outQueue, \
                               counter))
                p.order = k
                processes.append(p)
                p.start()
                counter = counter + arguments.chunk*chunksPerProcess
            else:
                p = mp.Process(target=scan_chunks, args=(arguments, \
                               arguments.chunk*(chunksPerProcess+1), \
                               outQueue, counter))
                p.order = k
                processes.append(p)
                p.start()
                remainingChunks = remainingChunks - 1
                counter = counter + arguments.chunk*(chunksPerProcess+1)
    elif numberChunks < arguments.threads:
        for l in range(0, numberChunks):
            p = mp.Process(target=scan_chunks, args=(arguments, \
                           arguments.chunk, outQueue, counter))
            p.order = l
            processes.append(p)
            p.start()
            counter = counter + arguments.chunk
    # Figure out the remainder
    size = frameCount%arguments.chunk
    if size != 0:
        p = mp.Process(target=scan_chunks, \
                       args=(arguments, size, outQueue, counter))
        p.order = numberChunks
        processes.append(p)
        p.start()
        counter = counter + size
    # Now that we have our processes list set up, we need to call
    # start() on each process in the list. The processes will pipe data to
    # the output Queue for writing output files
    print('Number of processes: ' + str(len(processes)))
    print('Extra process will scan the remaining ' + str(size) + ' frame(s).')
    print('Total Frames to Process: '+str(counter))
    #for p in processes:
    #    p.start()
    return outQueue


def scan_xtc(arguments):
    global framePattern, coordinatePattern, frameLength, coordinateLength, \
           allEntryPattern
    dataDict = {}
    out = start_processes(arguments)
    display(arguments, out)
    return

# {'frame':frame number, aggregation_number:average, counterion_name:average
# near micelle}
def display(args, outQueue):
    frameCount = find_frame_count(args, output=False)
    frameCounter = 0
    averageAggs = {}
    allData = {}
    getcontext().prec = 5
    time0 = 0
    time1 = 0
    while frameCounter < frameCount:
        time0 = time.time()
        data = outQueue.get()
        frameCounter = frameCounter + len(data)
        time1 = time.time()
        delta = time1-time0
        aggregations = find_aggregation_numbers_list(data, args.nomicelle, args.size)
        for i in data:
            averageAggs[i['time']] = {}
            averageAggs[i['time']]['MCL'] = average(aggregations[i['time']])
            nums = {}
            for micelle in i['micelles']:
                for molecule in micelle['near']['all']:
                    if molecule['name'] in nums.keys():
                        nums[molecule['name']] = nums[molecule['name']] + 1
                    else:
                        nums[molecule['name']] = 1
            for name in nums.keys():
                numMicelles = 0
                for micelle in i['micelles']:
                    if args.nomicelle == True:
                        numMicelles = numMicelles + 1
                    elif not len(micelle['molecules']) < args.size:
                        numMicelles = numMicelles + 1
                if numMicelles != 0:
                    averageAggs[i['time']][name] = float(nums[name])/numMicelles
                else:
                    averageAggs[i['time']][name] = 0
            allData[i['time']] = i
        percent = Decimal(len(allData))/frameCount*100
        try:
            remaining = (delta/len(data))*(frameCount-frameCounter)
        except ZeroDivisionError:
            remaining = frameCount
        print('\rPercent complete: ' + str(percent) + \
            '. Estimated completion time: ' + \
            str(time.ctime(time.time()+remaining)), \
            end='', file=sys.stdout)
        sys.stdout.flush()
    print('', file=sys.stdout)
    print('Writing data to files...', file=sys.stdout)
    output = open(args.output, 'w+b')
    for timel in sorted(averageAggs.keys()):
        output.write(str(timel) + ', ' + str(averageAggs[timel]['MCL'])+'\n')
    toutput = open(args.toutput, 'w+b')
    for timel in sorted(allData.keys()):
        toutput.write(printall(allData[timel], args.names, \
                      header='Time '+str(timel)))
    output.close()
    toutput.close()
    for name in args.names:
        if name in ['MCL', 'SKIP']:
            continue
        counterionFile = open(args.counterion+name+'.csv', 'w+b')
        for timel in sorted(averageAggs.keys()):
            try:
                counterionFile.write(str(timel) + ', ' + \
                    str(averageAggs[timel][name]) + '\n')
            except KeyError:
                counterionFile.write(str(timel) + ', 0\n')
        counterionFile.close()
    print('Done.', file=sys.stdout)
    return


# Main method
def main():
    global coordinateLength, frameLength
    arguments = get_parser().parse_args()
    coordinateLength = find_coordinate_length(arguments)
    frameLength = find_frame_length(arguments, framePattern)
    print('frameLength is: '+str(frameLength))
    # Ensure we get at least one file
    assert (arguments.file, arguments.tpr) != (None, None), 'You must provide \
                                              at least one .tpr, .trr, .xtc \
                                              file.'
    # Check xtclines
    if arguments.file != None:
        scan_xtc(arguments)
    # Check tprlines
    if arguments.tpr != None:
        scan_tpr(arguments)
    return

if __name__ == '__main__':
    main()

