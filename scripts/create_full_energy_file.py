#!/usr/bin/python3


# python script to create a data set of energies from a rs@md simulation
#
# uses:
# - the md engine to convert energy data files
# - numpy for computations etc
#


import numpy as np
import pandas as pd
import argparse
import os
import subprocess

import helper_functions as helper



def create_full_energies(*, nCycles, path, engine, filetype, types, output):
    print('reading energy data ...')
    ## some setup
    getFilenameEnergy = helper.make_get_filename_energy(FILETYPE=filetype, PATH=path)
    if types == None:
        if filetype == 'gmx':
            energyKeys = ['Potential', 'Kinetic-En.', 'Temperature', 'Box-X', 'Volume', 'Density']
        else:
            print('not-implemented-error: no default energy types were implemented')
            exit(1)
    else:
        energyKeys = types
        print(f'... using energy keys {energyKeys}')
    inputPipe = ''
    for key in energyKeys:
        inputPipe += f'{key}\n'
    inputPipe += '0\n'
    energyData = pd.DataFrame()
    lastTime = 0

    ## loop through all files and record remaining data
    for cycle in np.arange(0, nCycles+1):
        energyfile = getFilenameEnergy(cycle)

        if os.path.isfile( energyfile ):
            ## convert energy file to readable format
            completed= subprocess.run(['gmx', 'energy', '-f', energyfile, '-o', 'tmp.xvg', '-nobackup'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, input=inputPipe.encode() )
            
            ## read data
            tmp = np.loadtxt('tmp.xvg', comments=['#', '@'])

            ## get column specifier
            with open('tmp.xvg', 'r') as file:
                lines = file.readlines()
                cols = []
                for line in lines:
                    if '@ s' in line and 'legend' in line:
                        cols.append( line[line.find('"')+1:line.rfind('"')] )
            
            ## remove first line to omit t=0 value
            tmpData = pd.DataFrame(data=tmp[1:, :], columns=['time']+cols)
            
            ## set time correctly, then append to energyData
            tmpData['time'] = tmpData['time'] + lastTime
            if energyData.empty:
                energyData = tmpData
            else:
                energyData = energyData.append(tmpData, ignore_index=True, sort=False)
            lastTime = energyData['time'].iloc[-1]
            
        else:
            continue

    ## save data to file
    print('... done')
    energyData.to_csv(output, sep='\t')
    print(f'-> energy data for {lastTime} ps have been written to {output}')




if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='a python script to create a full data set \
                    of the energies during a rs@md simulation.')
    parser.add_argument('--nCycles', dest='nCycles', type=int, required=True, 
                    help='# of cycles to read')
    parser.add_argument('--path', dest='path', type=str, required=False,
                    help='path of the working directory (.)')
    parser.add_argument('--fileType', dest='fileType', type=str, required=False, 
                    help='which file type of topology to read (gmx)')
    parser.add_argument('--engine', dest='engine', type=str, required=False,
                    help='which md engine to use for file conversion (gmx)')
    parser.add_argument('--types', dest='types', type=str, required=False, nargs='+',
                    help='which energy terms to include in data set')
    parser.add_argument('-o', dest='output', type=str, required=False,
                    help='output file name (energies.data)')

    args = parser.parse_args()

    path = '.'
    if args.path != None:
        path = args.path
    
    fileType = 'gmx'
    if args.fileType != None:
        fileType = args.fileType
    
    engine = 'gmx'
    if args.engine != None:
        engine = args.engine

    outputfile = 'energies.data'
    if args.output != None:  
        outputfile = args.output

    create_full_energies(nCycles=args.nCycles, path=path, engine=engine, filetype=fileType, types=args.types, output=outputfile)