__author__='Martin Hutereau'
__version__='3'
###############################################################################
import os
import math
import re
import sys
from datetime import date
###############################################################################
'''
Python3 script to convert a CP2K restart file to .cif format by pulling
from the &COORD and &CELL blocks.
'''
#==============================================================================
#==============================================================================
#Read file until given string; returns False if string not found.

def read_until(file,targetstring):

    while True:
        line=file.readline()
        if targetstring in line:
            return line
        elif not line:
            return False

#==============================================================================
#==============================================================================
#Reads cell parameters from a CP2K CELL block. Assumes current position is 
#&CELL

def cell_params(file):

    avec=[float(i) for i in file.readline().split()[1:4]]
    bvec=[float(i) for i in file.readline().split()[1:4]]
    cvec=[float(i) for i in file.readline().split()[1:4]]

    a=math.sqrt(avec[0]*avec[0]+avec[1]*avec[1]+avec[2]*avec[2])
    b=math.sqrt(bvec[0]*bvec[0]+bvec[1]*bvec[1]+bvec[2]*bvec[2])
    c=math.sqrt(cvec[0]*cvec[0]+cvec[1]*cvec[1]+cvec[2]*cvec[2])

    #Store angles in radians for now
    alpha=math.acos(
        (bvec[0]*cvec[0]+bvec[1]*cvec[1]+bvec[2]*cvec[2])/(b*c)
        )
    beta=math.acos(
        (avec[0]*cvec[0]+avec[1]*cvec[1]+avec[2]*cvec[2])/(a*c)
        )
    gamma=math.acos(
        (avec[0]*bvec[0]+avec[1]*bvec[1]+avec[2]*bvec[2])/(b*a)
        )

    return [a,b,c,alpha,beta,gamma]

#==============================================================================
#==============================================================================
#Grab coordinates from restart file, check whether they are Cartesian or not.
#Assumes current file position is &COORD.

def coords(file):

    atoms=[]
    scaled=False
    while True:
        line=file.readline()
        try:
            l,x,y,z=line.split()
            atoms.append([l,float(x),float(y),float(z)])
        except ValueError:
            if 'SCALED' in line:
                scaled=True #coordinates are fractional
            elif '&END COORD' in line:
                break
    return atoms,scaled

#==============================================================================
#==============================================================================
#Read the given CP2K restart file and extract the cell parameters and atomic
#information. In a default restart file, the &SUBSYS block is not repeated and
#contains, in order, &CELL (optionally followed by &CELL_REF) and &COORD.

def cp2k_reader(filename):

    with open(filename,'r') as fin:

        #Read until &SUBSYS block and extract cell parameters
        read_until(fin,'&SUBSYS\n')
        fin.readline()
        cell_parameters=cell_params(fin)

        #Read until &COORD block and extract atom types and coordinates
        read_until(fin,'&COORD\n')
        atoms,scaled=coords(fin)

    return atoms,cell_parameters,scaled

#==============================================================================
#==============================================================================
#Builds the matrix for converting from cartesian to fractional coordinates
#for the given cell parameters

def transform_matrix(cell_params):

    a,b,c,xx,yy,zz=cell_params
    alpha=xx
    beta=yy
    gamma=zz
    omega=1-math.cos(alpha)*math.cos(alpha)-math.cos(beta)*math.cos(beta)\
    -math.cos(gamma)*math.cos(gamma)\
    +2*math.cos(alpha)*math.cos(beta)*math.cos(gamma)
    
    V=a*b*c*math.sqrt(
        1-math.cos(alpha)*math.cos(alpha)-math.cos(beta)*math.cos(beta)\
        -math.cos(gamma)*math.cos(gamma)\
        +2*math.cos(alpha)*math.cos(beta)*math.cos(gamma)
        )

    #Standard crystallographic transform from cartesian to fractional can be
    #carried out using the matrix defined below.
    t_matrix=[
    [
    1/a,
    -math.cos(gamma)/(a*math.sin(gamma)),
    (b*c*math.cos(gamma)*(math.cos(alpha)-math.cos(beta)\
    *math.cos(gamma))/math.sin(gamma)-b*c*math.cos(beta)*math.sin(gamma))/V
    ],
    [
    0,
    1/(b*math.sin(gamma)),
    -(a*c*(math.cos(alpha)\
        -math.cos(beta)*math.cos(gamma)))/(V*math.sin(gamma))
    ],
    [0,0,a*b*math.sin(gamma)/V]
    ]

    return t_matrix

#==============================================================================
#==============================================================================
#A simple numpyless function for multiplying a 3-vector (Python list) by a 3x3
#matrix (3x3 Python list)

def matrix_multiply(matrix,vector):

    transformed_vector=[0,0,0]
    for q in range(3):
        transformed_vector[q]=matrix[q][0]*vector[0]+matrix[q][1]*vector[1]\
        +matrix[q][2]*vector[2]

    return transformed_vector


#==============================================================================
#==============================================================================
#Makes use of translational symmetry to move all atoms inside the unit cell

def transform_coords(atoms,cell_params,scaled):
    
    cart_to_frac=transform_matrix(cell_params)
    for a in atoms:
        #If atoms are in cartesian coordinates, transform them first
        if not scaled:
            a[1],a[2],a[3]=matrix_multiply(cart_to_frac,a[1:4])
            #Make sure all atoms are inside the unit cell
        for i in range(3):
            while a[i+1]<0:
                a[i+1]+=1
            while a[i+1]>1:
                a[i+1]-=1

    return atoms

#==============================================================================
#==============================================================================
#This function writes all the types and positions to a cif file. This will
#contain the cell parameters, atom types, and positions, all in P1 symmetry.

def write_cif(outfile,cell_params,atoms):

    with open(outfile,'w') as file:

        file.write('#This file was created on {}\n'.format(date.today()))

        file.write('\ndata_image0\n'\
        +'_symmetry_space_group_name_H-M    "P 1"\n'\
        +'_symmetry_int_tables_number       1\n'\
        +'\n'\
        +'loop_\n'\
        +' _symmetry_equiv_pos_as_xyz\n'\
        +" 'x, y, z'\n"\
        +'\n')

        file.write('_cell_length_a\t{:.4f}\n'.format(cell_params[0])\
        +'_cell_length_b\t{:.4f}\n'.format(cell_params[1])\
        +'_cell_length_c\t{:.4f}\n'.format(cell_params[2])\
        +'_cell_angle_alpha\t{:.3f}\n'.format(cell_params[3]*180/math.pi)\
        +'_cell_angle_beta\t{:.3f}\n'.format(cell_params[4]*180/math.pi)\
        +'_cell_angle_gamma\t{:.3f}\n'.format(cell_params[5]*180/math.pi)\
        +'\n')

        file.write('loop_\n'\
        +' _atom_site_label\n'\
        +' _atom_site_type_symbol\n'\
        +' _atom_site_fract_x\n'\
        +' _atom_site_fract_y\n'\
        +' _atom_site_fract_z\n')

        #Use enumerate to generate unique atom labels
        for i,a in enumerate(atoms):
            file.write(
                ' {} {} {:.4f} {:.4f} {:.4f}\n'.format(
                    a[0]+str(i),a[0],a[1],a[2],a[3]
                    )
                )

#==============================================================================
#==============================================================================
input_file=sys.argv[1]
try:
    atoms,cell_parameters,scaled=cp2k_reader(input_file)
    atoms=transform_coords(atoms,cell_parameters,scaled)

    #Try to infer a reasonable name for the cif file based on input file name
    if input_file[:-10]=='-1.restart':
        new_filename=input_file[:-10]+'.cif'
    elif '.' in input_file:
        new_filename='.'.join(input_file.split('.')[:-1])+'.cif'
    else:
        new_filename=input_file+'.cif'
    write_cif(new_filename,cell_parameters,atoms)
    print('\n.cif file {} has been written.'.format(new_filename))
except Exception as e:
    print("Error encountered in file {}. Continuing\n".format(input_file))
    print(e)
