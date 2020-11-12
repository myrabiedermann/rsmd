#!/usr/bin/python3


# python script to compute radial distribution function 
#                   and cumulative number
#
# uses:
# - MDAnalysis for IO, selection of atoms / molecules, distance computations
# - numpy for computations etc
#
# notes:
# - if reference (i) or selection (j) consists of molecules, the center of mass is used as position
# - dimensions: distances -> Angstroem


import numpy as np
import argparse
import os

import MDAnalysis as mda
import helper_functions as helper


def compute_rdf(*, TOPFILE, TRJFILE, REF, SEL, BINWIDTH, BINMAX=None):

    print('computing rdf & cn from {}, {}'.format(TOPFILE, TRJFILE))
    print('... with reference = {}, selection = {}'.format(REF, SEL), flush=True )
    print(f'... bin width = {BINWIDTH}')
    print(f'... max bin = {BINMAX}')

    # import system & trajectory: .tpr, .gro/.xtc/...
    universe = mda.Universe(TOPFILE, TRJFILE)
    
    # set reference and selection
    reference = universe.select_atoms(REF)
    selection = universe.select_atoms(SEL)

    # rdf setup
    edges = np.arange(BINWIDTH/2, BINMAX + BINWIDTH/2, BINWIDTH)
    rdf, edges = np.histogram( [0], bins=edges )

    # make function to get distMatrix:
    distMatrix = helper.make_getDistMatrix( reference, selection )

    selectionDensity = 0
    nReferenceAtoms = reference.n_atoms

    # loop through trajectory:
    print('... collecting data', flush=True)
    nFrames = universe.trajectory.n_frames
    for ts in universe.trajectory: 
        dist = distMatrix( ts )
        rdf += np.histogram( np.ravel(dist), bins=edges )[0]    ## use only first return (--> histogram)

        selectionDensity += (dist.shape[1] / ts.volume)
        
    if nReferenceAtoms != dist.shape[0]:
        nReferenceAtoms = dist.shape[0]

    print('... analysing data', flush=True)

    selectionDensity /= nFrames
    sliceVolumina = (4./3.) * np.pi * (np.power(edges[1:],3) - np.power(edges[:-1], 3))

    # compute cumulative number 
    # pseudo: cn = sum([bin_height[i] for i in bin_indexes_to_integrate])
    # rdf needs to be only normalised by nReferenceAtoms and nFrames !
    cn = np.cumsum( rdf / (nReferenceAtoms * nFrames) )
    
    # rdf normalisation
    rdf = rdf / (nReferenceAtoms * nFrames * sliceVolumina * selectionDensity ) 

    return edges, rdf, cn


def compute_sum_rdf(*, nCycles, filepath, fileType, reference, selection, BINWIDTH, BINMAX):

    counter = 0
    counterPartLength = 0

    edges = np.arange(BINWIDTH/2, BINMAX + BINWIDTH/2, BINWIDTH)
    rdf = np.zeros( len(edges)-1 )
    cn = np.zeros( len(edges)-1 )
    radii = 0.5 * (edges[1:] + edges[:-1])

    getFilenameTop, getFilenameTrj = helper.make_get_filename(FILETYPE=fileType, PATH=filepath)

    tmp_rdf = np.zeros( len(edges)-1 )
    tmp_cn = np.zeros( len(edges)-1 )

    # loop through all files and record data
    for cycle in np.arange(0, nCycles+1):
        topfile = getFilenameTop(cycle)
        trjfile = getFilenameTrj(cycle)

        if os.path.isfile( trjfile ):
            # first: add previous data with correct weight
            counter += counterPartLength
            rdf += 1.0 * counterPartLength * tmp_rdf
            cn += 1.0 * counterPartLength * tmp_cn 
            print()
            # second: record new data
            counterPartLength = 1
            tmp_edges, tmp_rdf, tmp_cn = compute_rdf(TOPFILE=topfile, TRJFILE=trjfile, REF=reference, SEL=selection, BINWIDTH=BINWIDTH, BINMAX=BINMAX)
        else:
            counterPartLength += 1
            continue
    # add final data:
    counter += counterPartLength
    rdf += 1.0 * counterPartLength * tmp_rdf
    cn += 1.0 * counterPartLength * tmp_cn 

    rdf /= counter
    cn /= counter

    return edges, rdf, cn



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='a python script to compute the radial distribution function \
                    & cumulative number between a reference group and a selection group of atoms/molecules from a rs@md simulation.')
    parser.add_argument('--nCycles', dest='nCycles', type=int, required=True, 
                    help='# of cycles to read')
    parser.add_argument('--path', dest='PATH', type=str, required=False,
                    help='path of the working directory (.)')
    parser.add_argument('--fileType', dest='fileType', type=str, required=False, 
                    help='which file type of topology to read (gmx)')
    parser.add_argument('-o', dest='FILENAME', type=str, required=False,
                    help='output file name (analysis-rdf.data)')
    parser.add_argument('--ref', dest='REF', type=str, nargs='+', required=True,
                    help='reference specifier')
    parser.add_argument('--sel', dest='SEL', type=str, nargs='+', required=False,
                    help='selection specifier (is set equal to refernce if omitted)')
    parser.add_argument('--bin', dest='binwidth', type=float, required=False, 
                    help='bin width (0.02 A)')
    parser.add_argument('--binmax', dest='binmax', type=float, required=False, 
                    help='bin max (10.0 A)')

    args = parser.parse_args()

    PATH = '.'
    if args.PATH != None:
        PATH = args.PATH
    
    fileType = 'gmx'
    if args.fileType != None:
        fileType = args.fileType

    FILENAME = 'analysis-rdf.data'
    if args.FILENAME != None:
        FILENAME = args.FILENAME

    REF = ' '.join(args.REF)
    if args.SEL != None:
        SEL = ' '.join(args.SEL)
    else:
        SEL = REF

    binwidth = 0.02
    if args.binwidth != None:
        binwidth = args.binwidth
   
    binmax = 10.0
    if args.binmax != None:
        binmax = args.binmax


    edges, rdf, cn = compute_sum_rdf(nCycles=args.nCycles, filepath=PATH, fileType=fileType, reference=REF, selection=SEL, BINWIDTH=binwidth, BINMAX=binmax)


    # save data to file:
    radii = 0.5 * (edges[1:] + edges[:-1]) 
    with open(FILENAME, 'w') as output:
        output.write('#\n')
        output.write("# computation of radial distribution function and cumulative number as a function of radius\n")
        output.write('#\n')
        output.write('{:8} \t {:8} \t {:8}\n'.format('# r (Angs)', 'g(r)', 'cn'))
        for radius, gofr, cumnr in zip(radii, rdf, cn):
            output.write( f'{radius:8.2f} \t {gofr:8.3f} \t {cumnr:8.3f}\n' )

    print()
    print( f'g(r) & cumulative number data written to {FILENAME!r}' )

