#!/usr/bin/python3


# python script to compute the molecule evolution during a rs@md simulation
#
# uses:
# - numpy for computations etc
# - pandas for data storage
#
# notes:
# - dimensions: distances -> Angstroem



import numpy as np
import pandas as pd
import argparse
import os



def make_read_topology(*, FILETYPE):
    if FILETYPE == 'gmx' or FILETYPE == 'GMX' or FILETYPE == 'gromacs':
        
        def read_topology_gmx(cycle, file):
            df = pd.DataFrame(columns=['cycle'], index=[0])
            with open(file, 'r') as FILE:
                lines = [line.rstrip() for line in FILE]
                directiveMolecules = False
                for line in lines:
                    if 'molecules' in line:
                        directiveMolecules = True
                        continue
                    if ';' in line[:3] or '#' in line[:3]:
                        continue
                    if directiveMolecules:
                        content = line.split()
                        df[content[0]] = int( content[1] )
            return df

        return read_topology_gmx

    else:
        print("not-implemented-error: given FILETYPE = ", FILETYPE, 'not recognised or not implemented')



def make_get_filename(*, FILETYPE, PATH):
    if FILETYPE == 'gmx' or FILETYPE == 'GMX' or FILETYPE == 'gromacs':

        def get_filename_gmx(cycle):
            return PATH + f'/{cycle}.top'

        return get_filename_gmx

    else:
        print("not-implemented-error: given FILETYPE = ", FILETYPE, 'not recognised or not implemented')




def read_moleculeEvolution(*, path, nCycles, fileType='gmx'):

    ### read molecule types and numbers for all cycles ###
    moleculeData = pd.DataFrame(columns=['cycle'])

    get_filename = make_get_filename(FILETYPE=fileType, PATH=path)
    read_topology = make_read_topology(FILETYPE=fileType)

    for cycle in np.arange(nCycles):
        filename = get_filename(cycle)
        if not os.path.isfile( filename ):
            continue
        tmp = read_topology(cycle, filename)
        tmp['cycle'] = cycle
        moleculeData = moleculeData.append(tmp, ignore_index=True, sort=True)

    moleculeData.fillna(0, inplace=True)

    return moleculeData




if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='a python script to compute the molecule evolution from a rs@md simulation.')
    parser.add_argument('--path', dest='PATH', type=str, required=False,
                    help='path of the working directory')
    parser.add_argument('--nCycles', dest='nCycles', type=int, required=True, 
                    help='# of cycles to read')
    parser.add_argument('--fileType', dest='fileType', type=str, required=False, 
                    help='which file type of topology to read (gmx)')
    
    
    args = parser.parse_args()


    fileType = 'gmx'
    if args.fileType != None:
        fileType = args.fileType


    PATH = '.'
    if args.PATH != None:
        PATH = args.PATH


    data = read_moleculeEvolution(path=PATH, nCycles=args.nCycles, fileType=fileType)

    data.to_csv(f'{PATH}/moleculeEvolution.data', sep='\t')
    print(f'saved data to {PATH}/moleculeEvolution.data')
    
