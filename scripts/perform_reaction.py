#!/usr/bin/python3


# python script to perform a reaction (as is done in a rs@md simulation)
#
# uses:
# - MDAnalysis for IO
# - numpy for computations etc
# - pandas for data handling
#



import numpy as np
import pandas as pd
import argparse
import os

import MDAnalysis as mda


def read_reactions_file(*, file):
    print('... reading reaction template from file', file)
    # open file for reading
    lines = []
    with open(file, 'r') as inFile:
        lines = inFile.readlines()

    # go through line by line
    currentDirective = ""
    transitionTable = []
    translations = []

    for line in lines:
        line = line.strip()
        ## interpret line and save accordingly
        if not line:
            continue;
        elif line[0] == '#':
            ## comment line --> skip
            continue;
        elif '[' in line and ']' in line:
            ## begin of new directive
            currentDirective = line[line.find('[')+1:line.find(']')].strip()
        else:
            ## first: check if part of line is a comment
            if '#' in line:
                index = line.find("#");
                line = line[:index]
            ## second: parse content of line
            elif "products" in currentDirective:
                ## # molNr      molName     atomName    atomNr  origin->molNr   origin->atomNr
                inputs = line.split()
                transitionTable.append([int(inputs[0])-1, int(inputs[3])-1, inputs[1], inputs[2], int(inputs[4])-1, int(inputs[5])-1])
            elif "translations" in currentDirective:
                inputs = line.split()
                translations.append([int(inputs[0])-1, int(inputs[1])-1, int(inputs[2])-1, int(inputs[3])-1, float(inputs[4])])
    return transitionTable, translations



def perform_reaction(*, coordinatesFile, reactandIDs, transitionTable, translations, inputTranslations):
    print('... performing transitions')
    
    ## read coordinates of atoms
    data = pd.read_csv(coordinatesFile, skiprows=2, skipfooter=1, index_col=None,\
            sep=' ', skipinitialspace=True, names=['res', 'atom', 'id', 'x', 'y', 'z', 'vx', 'vy', 'vz'], engine='python')
    data['resid'] = [int(''.join(filter(str.isdigit, tmp))) for tmp in data['res'] ]
    data['resName'] = [''.join(filter(str.isalpha, tmp)) for tmp in data['res'] ]
    dataNew = data.copy()


    ## get reactand molecules & drop them from dataNew
    reactandIXs = []
    for id in reactandIDs:
        reactandIXs.append( data['res'] == id )

    allIXs = []
    for ix in reactandIXs:
        allIXs += np.argwhere(ix).flatten().tolist()
    allIXs = np.array(allIXs, dtype=int).flatten()

    dataNew.drop( index=allIXs, inplace=True )


    ## perform transitions, i.e. add new atom to newData
    highestMolID = int(''.join(filter(str.isdigit, data['res'].iloc[-1])))
    productIDs = [] 
    for transition in transitionTable:  ## contains: newResIx   newAtomIx   newResName   newAtomName   oldResIx   oldAtomIx
        newResIx, newAtomIx, newResName, newAtomName, oldResIx, oldAtomIx = tuple(transition)
        newAtom = dict()
        ## assign resid/resname
        newAtom['res'] = f'{(newResIx+highestMolID+1):>5d}{newResName}'
        newAtom['resid'] = newResIx + highestMolID + 1
        newAtom['resName'] = newResName
        ## assign newAtomName
        newAtom['atom'] = newAtomName
        ## assign newAtomIx
        newAtom['id'] = newAtomIx + 1
        ## assign positions
        newAtom['x'] = data[reactandIXs[oldResIx]]['x'].iloc[oldAtomIx]
        newAtom['y'] = data[reactandIXs[oldResIx]]['y'].iloc[oldAtomIx]
        newAtom['z'] = data[reactandIXs[oldResIx]]['z'].iloc[oldAtomIx]
        ## assign velocities
        newAtom['vx'] = data[reactandIXs[oldResIx]]['vx'].iloc[oldAtomIx]
        newAtom['vy'] = data[reactandIXs[oldResIx]]['vy'].iloc[oldAtomIx]
        newAtom['vz'] = data[reactandIXs[oldResIx]]['vz'].iloc[oldAtomIx]
        ## append new atom to new topology
        dataNew = dataNew.append(newAtom, ignore_index=True)
        productIDs.append(newResIx + highestMolID + 1)
    productIDs = np.unique(productIDs)


    ## perform translations
    for i, transl in enumerate(translations):
        ## contains: molix  atomix  molix  atomix  value
        mol1ix, atom1ix, mol2ix, atom2ix, value = tuple(transl)
        mol1id = productIDs[mol1ix]
        mol2id = productIDs[mol2ix]
        atom1id = atom1ix + 1
        atom2id = atom2ix + 1
        print('... performing translations for atom:')
        print( dataNew.loc[ (dataNew['resid'] == mol1id) & (dataNew['id'] == atom1id) ] )
        print('... towards / away from atom')
        print( dataNew.loc[ (dataNew['resid'] == mol2id) & (dataNew['id'] == atom2id) ] )
        if inputTranslations != None:
            value = inputTranslations[i]
            print('    using command line input translation value:', "%.3f" % value)
        ## get position vectors:
        atom1Pos = np.zeros((3))
        atom2Pos = np.zeros((3))
        for i, p in enumerate(['x', 'y', 'z']):
            atom1Pos[i] = dataNew.loc[ (dataNew['resid'] == mol1id) & (dataNew['id'] == atom1id), p ]
            atom2Pos[i] = dataNew.loc[ (dataNew['resid'] == mol2id) & (dataNew['id'] == atom2id), p ]
        ## compute distance & connection axis
        print('    distance before:', "%.3f" % np.linalg.norm(atom2Pos - atom1Pos))
        connection = atom2Pos - atom1Pos
        connection /= np.linalg.norm(connection)
        connection = connection * value
        ## translate atom
        for i, p in enumerate(['x', 'y', 'z']):
            dataNew.loc[ (dataNew['resid'] == mol1id) & (dataNew['id'] == atom1id), p] += connection[i]
        ## check distance
        for i, p in enumerate(['x', 'y', 'z']):
            atom1Pos[i] = dataNew.loc[ (dataNew['resid'] == mol1id) & (dataNew['id'] == atom1id), p ]
            atom2Pos[i] = dataNew.loc[ (dataNew['resid'] == mol2id) & (dataNew['id'] == atom2id), p ]
        print('    distance after: ', "%.3f" % np.linalg.norm(atom2Pos - atom1Pos)) 

    ## sort topology alphabetically
    dataNew.sort_values(by=['resName', 'id'], axis='index', inplace=True)

    return dataNew
   


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='a python script to perform a \
                    reaction, i.e. change the topology.')
    parser.add_argument('-i', dest='reactionFile', type=str, required=True, 
                    help='input file for the reaction template')
    parser.add_argument('-f', dest='coordinatesFile', type=str, required=True,
                    help='input coordinates file')
    parser.add_argument('--reactands', dest='reactands', type=str, nargs='+', required=True,
                    help='IDs of reactand molecules')
    parser.add_argument('--fileType', dest='fileType', type=str, required=False, 
                    help='which file type of topology to read (gmx)')
    parser.add_argument('--translations', dest='translations', type=float, nargs='+', required=False,
                    help='overwrite translation values from reaction template file')
    parser.add_argument('-o', dest='output', type=str, required=False,
                    help='output file name (reacted.gro)')

    args = parser.parse_args()

    fileType = 'gmx'
    if args.fileType != None:
        fileType = args.fileType

    outputfile = 'reacted.gro'
    if args.output != None:
        outputfile = args.output

    transitions, translations = read_reactions_file(file=args.reactionFile)
    data = perform_reaction(coordinatesFile=args.coordinatesFile, reactandIDs=args.reactands, \
                    transitionTable=transitions, translations=translations, inputTranslations=args.translations)

    ## get first two lines and last line of coordinatesFile
    firstlines = []
    lastline = []
    with open(args.coordinatesFile, 'r') as old:
        tmp = old.readlines()
        firstlines = tmp[:2]
        lastline = tmp[-1]

    ## write to file
    with open(outputfile, 'w') as output:
        output.write(firstlines[0])
        output.write(firstlines[1])
        for i, row in data.iterrows():
            output.write(f"{row['resid']:>5}{row['resName']:<5s}{row['atom']:>5s}{row['id']:5d}{row['x']:8.3f}{row['y']:8.3f}{row['z']:8.3f}{row['vx']:8.4f}{row['vy']:8.4f}{row['vz']:8.4f}\n")
        output.write(lastline)

    print('... reacted coordinates have been written to ', outputfile)