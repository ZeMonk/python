###############################################################################
import os
import math
import sys
from datetime import date
###############################################################################
'''
Python3 script to convert a CP2K restart file to .cif format by pulling
from the &COORD and &CELL blocks. Much of the code can be made more elegant
using Numpy, but this introduces a dependency.
'''
###############################################################################

#Read the passed file until the given target string.

def read_until(file,targetstring):

    while True:
        line=file.readline()
        if targetstring in line:
            return

#==============================================================================
#==============================================================================
#Grab cell parameters. Assumes position in restart file is currently at the 
#cell parameters block &CELL.

def cell_params(file):

    #Store the cell vectors as lists of floats.
    avec=[float(i) for i in file.readline().split()[1:4]]
    bvec=[float(i) for i in file.readline().split()[1:4]]
    cvec=[float(i) for i in file.readline().split()[1:4]]

    #Calculate the magnitudes of the cell vectors.
    a=math.sqrt(avec[0]*avec[0]+avec[1]*avec[1]+avec[2]*avec[2])
    b=math.sqrt(bvec[0]*bvec[0]+bvec[1]*bvec[1]+bvec[2]*bvec[2])
    c=math.sqrt(cvec[0]*cvec[0]+cvec[1]*cvec[1]+cvec[2]*cvec[2])

    #Calculate cell angles, in degrees, from cell vectors and their magnitudes.
    # alpha=math.acos(
    #     (bvec[0]*cvec[0]+bvec[1]*cvec[1]+bvec[2]*cvec[2])/(b*c)
    #     )*180/math.pi
    # beta=math.acos(
    #     (avec[0]*cvec[0]+avec[1]*cvec[1]+avec[2]*cvec[2])/(a*c)
    #     )*180/math.pi
    # gamma=math.acos(
    #     (avec[0]*bvec[0]+avec[1]*bvec[1]+avec[2]*bvec[2])/(b*a)
    #     )*180/math.pi
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
#Extract the atom coordinates and check whether they are Cartesian or not. This
#assumes that the current file position is the coordinate block &COORD.

def coords(file):

    atoms=[]
    scaled=False #The default in CP2K restart is to print Cartesian coordinates.

    while True:
        line=file.readline()

        try:
            l,x,y,z=line.split() #Atom label, x coord, y coord, z coord.
            atoms.append([l,float(x),float(y),float(z)])

        except ValueError:
            if 'SCALED' in line:
                scaled=True
                #When the coordinates are fractional, CP2K indicates as much 
                #with the flag SCALED T in the coordinates block.

            elif '&END COORD' in line:
                break
                #Continue reading until the end of the coordinates block.

    return atoms,scaled

#==============================================================================
#==============================================================================
#Function for reading a CP2K restart file and pulling the atom coordinates
#and cell parameters.

def cp2kreader(filename):

    with open(filename,'r') as file:

        read_until(file,'&SUBSYS\n') #Read file until cell params block.
        
        #CP2K users may sometimes include a reference cell during calculations.
        #Check whether this is present and if so, skip it.
        check=file.readline()
        if 'CELL_REF' in check:
            read_until(file,'&CELL\n')

        cell_parameters=cell_params(file) #Grab cell params.

        read_until(file,'&COORD\n') #Read file until coordinates block
        atoms,scaled=coords(file) #Grab the coordinates. These may be scaled
        #or not, as described by the variable 'scaled'.

    return atoms,cell_parameters,scaled

#==============================================================================
#==============================================================================
#If the coordinates printed in the CP2K restart file are Cartesian, then they
#need to be converted to fractional using the corresponding transformation
#matrix. This function generates that matrix from the cell parameters.

def cart_to_frac_matrix(cell_params):

    a,b,c,alpha,beta,gamma=cell_params

    #Standard crystallograph transformation matrix involves 2 quantities,
    #omega and the cell volume, which are both dependent on the cell params.

    omega=1-math.cos(alpha)*math.cos(alpha)-math.cos(beta)*math.cos(beta)\
    -math.cos(gamma)*math.cos(gamma)\
    +2*math.cos(alpha)*math.cos(beta)*math.cos(gamma)
    
    V=a*b*c*math.sqrt(
        1-math.cos(alpha)*math.cos(alpha)-math.cos(beta)*math.cos(beta)\
        -math.cos(gamma)*math.cos(gamma)\
        +2*math.cos(alpha)*math.cos(beta)*math.cos(gamma)
        )

    #The transformation matrix is as defined below.
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
#A simple function for multiplying a 3x1 vector by a 3x3 matrix without using
#Numpy.

def matrix_multiply(matrix,vector):

    transformed_vector=[0,0,0]

    for q in range(3):
        transformed_vector[q]=matrix[q][0]*vector[0]+matrix[q][1]*vector[1]\
        +matrix[q][2]*vector[2]

    return transformed_vector


#==============================================================================
#==============================================================================
#Function which shifts all atoms to lie within [0,1] in fractional coordinates.
#If these are Cartesian, also transforms them to fractional.

def convert_coordinates(atoms,cell_params,scaled):
    
    t_matrix=cart_to_frac_matrix(cell_params)

    for a in atoms:
        if not scaled:
            a[1],a[2],a[3]=matrix_multiply(t_matrix,a[1:4])

        for i in range(3):
            
            while a[1+i]<0.0:
                a[1+i]+=1

            while a[1+i]>1.0:
                a[1+i]-=1

    #Let user know if coordinates have been transformed.
    if not scaled:
        print('Coordinates have been converted from Cartesian to fractional.')

    return atoms

#==============================================================================
#==============================================================================
#Function for writing the cell parameters and atomic coordinates in a cif file.

def makecif(outfile,cell_params,atoms):

    with open(outfile,'w') as file:

        file.write('#This file was created on {}\n'.format(date.today()))

        #Pre-position block. CP2K typically doesn't make use of symmetry,
        #so the output symmetry is P1.
        file.write('\ndata_image0\n'\
        +'_symmetry_space_group_name_H-M "P 1"\n'\
        +'_symmetry_int_tables_number 1\n'\
        +'\n'\
        +'loop_\n'\
        +' _symmetry_equiv_pos_as_xyz\n'\
        +" 'x, y, z'\n"\
        +'\n')

        #Cell parameter block.
        file.write('_cell_length_a {:.4f}\n'.format(cell_params[0])\
        +'_cell_length_b {:.4f}\n'.format(cell_params[1])\
        +'_cell_length_c {:.4f}\n'.format(cell_params[2])\
        +'_cell_angle_alpha {:.4f}\n'.format(cell_params[3]*180/math.pi)\
        +'_cell_angle_beta {:.4f}\n'.format(cell_params[4]*180/math.pi)\
        +'_cell_angle_gamma {:.4f}\n'.format(cell_params[5]*180/math.pi)\
        +'\n')

        #Position block.
        file.write('loop_\n'\
        +' _atom_site_label\n'\
        +' _atom_site_type_symbol\n'\
        +' _atom_site_fract_x\n'\
        +' _atom_site_fract_y\n'\
        +' _atom_site_fract_z\n')

        for i,a in enumerate(atoms):
            label=a[0]+str(i)
            #Each atom is given a label number+atom type

            file.write(
                ' {} {} {:.4f} {:.4f} {:.4f}\n'.format(
                    label,a[0],a[1],a[2],a[3]
                    )
                )

#==============================================================================
#==============================================================================
#Attempt to recognise path and file names from command line inputs.

if len(sys.argv)>2:
    print('''Unexpected command line arguments given. Only the first argument
(a path or file name) will be processed, the rest will be ignored.\n''')

path=sys.argv[1]
filename=False
#Check whether first command line argument is indeed a path.
try:
    os.chdir(path)
except NotADirectoryError:
    filename=sys.argv[1]

#If a file name is passed, this becomes the only file the script runs on.
if filename:
    targets=[filename]

#If a directory path is passed, the script will attempt to convert every
#restart file it finds.
else:
    targets=[f for f in os.listdir() if f[-10:]=='-1.restart']

for f in targets:
    try:
        atoms,cell_parameters, scaled=cp2kreader(f)
        atoms=convert_coordinates(atoms,cell_parameters,scaled)
        
        #Name of new file is derived from original.
        new_filename=f[:-10]+'.cif'
        makecif(new_filename,cell_parameters,atoms)

        print('\nCif file {} has been written using CP2K restart file {}.'\
            .format(new_filename,f))

    except Exception as e:
        print(
            "Error '{}' encountered in file {}. Skipping file...\n".format(e,f)
            )
