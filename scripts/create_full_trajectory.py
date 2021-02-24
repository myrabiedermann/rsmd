#!/usr/bin/python3


# python script to create a .xyz trajectory from a rs@md simulation
#
# uses:
# - MDAnalysis for IO
# - numpy for computations etc
#


import numpy as np
import argparse
import os

import MDAnalysis as mda
import MDAnalysis.transformations as mdatrans
import MDAnalysis.topology.guessers as mdaguess
import helper_functions as helper



def create_full_trajectory(*, nCycles, path, filetype, output):
    ## some setup
    getFilenameTop, getFilenameTrj = helper.make_get_filename(FILETYPE=filetype, PATH=path)
    getReactiveAtomIndices = helper.make_get_reactive_atoms(FILETYPE=filetype, PATH=path)
    
    frameCounter = 0

    write_frame = helper.make_get_write(filename=output)
    FILE = open(output, 'w')

    ## get initial data from cycle 0
    universe = mda.Universe( getFilenameTop(0), getFilenameTrj(0) )
    atomNames = universe.atoms.names
    resNames = universe.atoms.resnames
    resIDs = universe.atoms.resnums
    atomOrder = universe.atoms.ix
    dt = universe.trajectory.dt
    ## important: sort topology such that resnames are in alphabetical order (because rs@md does this, too...)
    # first: sort resnames alphabetically and adapt resIDs, atomNames
    sortedIndices = resNames.argsort(kind='stable')
    atomNames = atomNames[sortedIndices]
    resNames = resNames[sortedIndices]
    resIDs = resIDs[sortedIndices]
    atomOrderRearrangedAlphabetically = atomOrder[sortedIndices]
    # second: renumber residues
    counterResID = 1
    newResIDs = resIDs.copy()
    for ix in range(len(resIDs)-1):
        newResIDs[ix] = counterResID
        if resIDs[ix] != resIDs[ix+1]:
            counterResID += 1
    resIDs = newResIDs

	## create initial coordination file where all residues have been unwrapped and save it to file
    #for res in universe.atoms.residues:
    #    ag = res.atoms
    #    transform = mdatrans.unwrap(ag)
    #    universe.trajectory.add_transformations(transform)
    if not hasattr(universe, 'bonds'):
        guessed_bonds = mdaguess.guess_bonds(universe.atoms, universe.atoms.positions, universe.trajectory[0].dimensions)
        universe.add_TopologyAttr('bonds', guessed_bonds)
    transform = mdatrans.unwrap(universe.atoms)
    universe.trajectory.add_transformations(transform)
    outputinitial = 'initial' + '.' + output.split('.')[-1]
    initialFile = open(outputinitial, 'w')
    box = universe.trajectory[0].dimensions
    positions = universe.atoms.positions
    sortedPositions = positions[sortedIndices]
    write_frame(filestream=initialFile, title=f'rs@md t={0:9.2f} step= {0}', resnames=resNames, resnums=resIDs, names=atomNames, positions=sortedPositions, box=box)
    initialFile.close()
    print(f'-> initial frame has been written to {outputinitial}')

    ## write trajectory to file frame by frame
    for ts in universe.trajectory:  
        positions = universe.atoms.positions
        box = ts.dimensions
        ## sort positions before writing them
        sortedPositions = positions[sortedIndices]
        write_frame(filestream=FILE, title=f'rs@md t={frameCounter*dt:9.2f} step=  {frameCounter}', resnames=resNames, resnums=resIDs, names=atomNames, positions=sortedPositions, box=box)
        frameCounter += 1

    ## loop through all files and record remaining data
    firstReactiveCycle = True
    for cycle in np.arange(1, nCycles+1):
        topfile = getFilenameTop(cycle)
        trjfile = getFilenameTrj(cycle)

        if os.path.isfile( trjfile ):
            ## get reaction infos
            reactantsIx, productsIx = getReactiveAtomIndices(cycle)
            
            ## attention if firstReactiveCycle: their might have been a change in atomOrder between cycles 0 -> cycle
            ## due to reordering molecules in alphabetical order
            ## need to account for that by translating reactantsIx accordingly
            if firstReactiveCycle:
                reactantsIx = np.array( [np.argwhere(atomOrderRearrangedAlphabetically == rix)[0,0] for rix in reactantsIx] )
                firstReactiveCycle = False                

            ## important: you need to go through all transitions in an ordered fashion with respect to the product indices 
            ## (from small to larger iy)
            sortedIndices = productsIx.argsort()
            reactantsIx = reactantsIx[sortedIndices]
            productsIx = productsIx[sortedIndices]
            
            ## apply reactions to atomOrder
            ## i.e. change positions of reactants to products

            ## ... first: get entries at reactant positions and remove them from list
            reactantEntries = [ atomOrder[x] for x in reactantsIx ]
            for entry in reactantEntries:
                atomOrder = np.delete(atomOrder, np.argwhere(atomOrder==entry), axis=0)
           
            ## ... second: put entries back at new positions given by product positions
            for iy, entry in zip(productsIx, reactantEntries):
                atomOrder = np.insert(atomOrder, iy, entry, axis=0)

            ## get sorted indices for new atomOrder
            sortedIndices = atomOrder.argsort()

            ## import trajectory: .tpr, .gro/.xtc/...
            universe = mda.Universe(topfile, trjfile)

            ## write trajectory to file frame by frame
            for ts in universe.trajectory[1:]:  
                positions = universe.atoms.positions
                box = ts.dimensions
                
                ## sort positions before writing them
                sortedPositions = positions[sortedIndices]
                write_frame(filestream=FILE, title=f'rs@md t={frameCounter*dt:9.2f} step=  {frameCounter}', resnames=resNames, resnums=resIDs, names=atomNames, positions=sortedPositions, box=box)
                frameCounter += 1

        else:
            continue

    FILE.close()

    print(f'-> a total of {frameCounter} frames have been written to {output}')



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='a python script to create a full, per element \
                    trajectory in .xyz format from a rs@md simulation.')
    parser.add_argument('--nCycles', dest='nCycles', type=int, required=True, 
                    help='# of cycles to read')
    parser.add_argument('--path', dest='path', type=str, required=False,
                    help='path of the working directory (.)')
    parser.add_argument('--fileType', dest='fileType', type=str, required=False, 
                    help='which file type of topology to read (gmx)')
    parser.add_argument('-o', dest='output', type=str, required=False,
                    help='output file name (trajectory.gro)')

    args = parser.parse_args()

    path = '.'
    if args.path != None:
        path = args.path
    
    fileType = 'gmx'
    if args.fileType != None:
        fileType = args.fileType

    outputfile = 'trajectory.gro'
    if args.output != None:
        outputfile = args.output

    create_full_trajectory(nCycles=args.nCycles, path=path, filetype=fileType, output=outputfile)
