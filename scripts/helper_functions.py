#!/usr/bin/python3

import MDAnalysis as mda
import MDAnalysis.lib.distances as mdaDist # #distance_array
import numpy as np


def make_get_reactive_atoms(*, FILETYPE, PATH):
    if FILETYPE == 'gmx' or FILETYPE == 'GMX' or FILETYPE == 'gromacs':
       
        def get_indices(cycle):
            reactandIndices, productIndices = ([], [])

            with open(PATH + f'/{cycle}.reactands.ndx', 'r') as FILE:
                content = FILE.read()
                content = content[content.find(']')+1:]
                reactandIndices = [int(x)-1 for x in content.split()]   ## subtract -1 to get index instead of ID
        
            with open(PATH + f'/{cycle}.products.ndx', 'r') as FILE:
                content = FILE.read()
                content = content[content.find(']')+1:]
                productIndices = [int(x)-1 for x in content.split()]
            return np.array(reactandIndices), np.array(productIndices)
        
        return get_indices

    else: 
        print("not-implemented-error: given filetype = ", FILETYPE, 'not recognised or not implemented')
        exit(1)



def make_get_filename(*, FILETYPE, PATH):
    if FILETYPE == 'gmx' or FILETYPE == 'GMX' or FILETYPE == 'gromacs':

        def get_filename_top(cycle):
            return PATH + f'/{cycle}-md.gro'
        def get_filename_trj(cycle):
            return PATH + f'/{cycle}-md.xtc'
        
        return get_filename_top, get_filename_trj

    else:
        print("not-implemented-error: given filetype = ", FILETYPE, 'not recognised or not implemented')
        exit(1)



def make_get_filename_energy(*, FILETYPE, PATH):
    if FILETYPE == 'gmx' or FILETYPE == 'GMX' or FILETYPE == 'gromacs':

        def get_filename_edr(cycle):
            return PATH + f'/{cycle}-md.edr'
            
        return get_filename_edr
    
    else: 
        print("not-implemented-error: given filetype =", FILETYPE, 'not recognised or not implemented')
        exit(1)



def containsMolecules( selection ):
    # if n_atoms != residues.n_atoms --> more atoms in residues than in selection, selection doesn't consist of molecules
    if selection.n_atoms != selection.residues.n_atoms:
        return False
    # if a residue contains only one atom --> selection doesn't consist of molecules
    elif [ res.atoms.n_atoms for res in selection.residues ] == [1] * selection.n_residues:
        return False
    # else --> selection consists of molecules
    else:
        return True



def make_getDistMatrix( reference, selection):
    # if reference/selection contains molecules -> use CoM, else -> use positions
    if containsMolecules(reference) and containsMolecules(selection):
        print('... note: using center of mass positions for reference and selection')
        print('... warning: remember that center of mass computation only works correctly with whole (unwrapped) molecules !')

        def getDistMatrix_com_com( frame ):
            box = frame.dimensions
            ref_coor = np.array([res.atoms.center_of_mass() for res in reference.residues], dtype=np.float32)
            sel_coor = np.array([res.atoms.center_of_mass() for res in selection.residues], dtype=np.float32)
            dist = mdaDist.distance_array(ref_coor, sel_coor, box) 
            return dist
        
        return getDistMatrix_com_com

    elif containsMolecules(reference) and not containsMolecules(selection):
        print('... note: using center of mass positions for reference')
        print('... warning: remember that center of mass computation only works correctly with whole (unwrapped) molecules !')
        
        def getDistMatrix_com_direct( frame ):
            box = frame.dimensions
            ref_coor = np.array([res.atoms.center_of_mass() for res in reference.residues], dtype=np.float32)
            sel_coor = selection.positions
            dist = mdaDist.distance_array(ref_coor, sel_coor, box) 
            return dist

        return getDistMatrix_com_direct

    elif not containsMolecules(reference) and containsMolecules(selection):
        print('... note: using center of mass positions for selection')
        print('... warning: remember that center of mass computation only works correctly with trajectory of whole molecules !')

        def getDistMatrix_direct_com( frame ):
            box = frame.dimensions
            ref_coor = reference.positions
            sel_coor = np.array([res.atoms.center_of_mass() for res in selection.residues], dtype=np.float32)
            dist = mdaDist.distance_array(ref_coor, sel_coor, box) 
            return dist

        return getDistMatrix_direct_com

    elif not containsMolecules(reference) and not containsMolecules(selection):
        print('... note: not using center of mass positions at all')

        def getDistMatrix_direct_direct( frame ):
            box = frame.dimensions
            ref_coor = reference.positions
            sel_coor = selection.positions
            dist = mdaDist.distance_array(ref_coor, sel_coor, box) 
            return dist

        return getDistMatrix_direct_direct



def make_get_write(*, filename):
    if filename[-3:] == 'gro':
        print('setting frame writer to .gro format')
        return write_gro
    else:
        print('setting frame writer to .xyz format')
        return write_xyz

# .xyz file format contains:
# nAtoms
# title line
# one line per atom with:
# - name
# - position (in A, x y z in 3 columns, free format)
def write_xyz(*, filestream, title, names, positions, box=None):
    if( len(names) != len(positions) ):
        raise RuntimeError("length of names and positions doesn't match")

    filestream.write( f'{len(names)}\n' )
    filestream.write( f'{title}\n' )

    for name, pos in zip(names, positions):
        filestream.write( f'{name:>8} {pos[0]:11.3f}{pos[1]:11.3f}{pos[2]:11.3f}\n' )


# .gro file format contains: 
# title line
# nAtoms
# one line per atom with: 
# - residue number (5 positions, integer)
# - residue name (5 characters)
# - atom name (5 characters)
# - atom number (5 positions, integer)
# - position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
# - velocity (in nm/ps (or km/s), x y z in 3 columns, each 8 positions with 4 decimal places) 
# box dimensions
def write_gro(*, filestream, title, names, positions, box):
    if( len(names) != len(positions) ):
        raise RuntimeError("length of names and positions doesn't match")

    filestream.write( f'{title}\n' )
    filestream.write( f' {len(names)}\n' )

    for i, (name, pos) in enumerate(zip(names, positions)):
        filestream.write( f'{i+1:>5}{name:<5}{name:>5}{i+1:>5}{pos[0]*0.1:8.3f}{pos[1]*0.1:8.3f}{pos[2]*0.1:8.3f}\n' )
    filestream.write( f'   {box[0]*0.1:8.5f}   {box[1]*0.1:8.5f}   {box[2]*0.1:8.5f}\n' )
