__author__='Martin Hutereau'
__version__='2'
###############################################################################
import sys
import math
from datetime import date
###############################################################################
'''
Python3 script to convert an xyz file into a cif file using some input
cell parameters. If the file contains multiple coordinate blocks, by default
the final one is used, but this can be chosen as an optional parameter.
'''
#==============================================================================
#==============================================================================
#Function for reading a single xyz coordinate block. This assumes proper xyz
#format of natom, blank line, then a number of coordinate lines equal to natom.

def coords(file):

    #Read number of atoms
    natoms=int(file.readline())
    file.readline()

    #Read atom types and coordinates
    atoms=[]
    for i in range(natoms):
        l=file.readline().split()
        atoms.append([str(l[0]),float(l[1]),float(l[2]),float(l[3])])

    return atoms

#==============================================================================
#==============================================================================
#Function for reading an xyz file and returning the requested set of 
#coordinates. By default, returns the final one.

def readxyz(filename,step=0):

    with open(filename,'r') as fin:

        #Grab the number of atoms and check how many lines are present in
        #the file.
        natoms=int(fin.readline())
        i=1

        #Continue until EOF, which returns '', which evaluates as False.
        while fin.readline():
            i+=1
        #Check to make sure that each block of coordinates is complete, round
        #down and warn the user otherwise.
        if not i%natoms:
            print('Issue detected with .xyz file. Expecting an integer'+\
                'multiple of {} ({} atoms) lines, but found {} instead.'+\
                'This may cause issues.'.format(natoms+2,natoms,i))
            numsteps=(i-i%natoms)/natoms
        else:
            numsteps=i/(natoms+2)

        #Return to start of file.
        fin.seek(0)

        #Read ahead to requested coordinate block and extract it
        if step>numsteps:
            print("The requested step is higher than number of frames found.")
        if step:
            for j in range(int(step-1)*(natoms+2)):
                fin.readline()
            atoms=coords(fin)
        else:
            print(int(numsteps-1)*(natoms+2),natoms,numsteps)
            for j in range(int(numsteps-1)*(natoms+2)):
                fin.readline()
            atoms=coords(fin)

        return atoms

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
#Function which shifts all atoms to lie within [0,1] in fractional coordinates.
#If these are Cartesian, also transforms them to fractional.

def convert_coordinates(atoms,cell_params):
    
    cart_to_frac=transform_matrix(cell_params)

    for a in atoms:
        a[1],a[2],a[3]=matrix_multiply(cart_to_frac,a[1:4])
        for i in range(3):
            while a[1+i]<0.0:
                a[1+i]+=1
            while a[1+i]>1.0:
                a[1+i]-=1

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
#Check number of command line inputs and convert angles to radians.

if len(sys.argv)==8:
    cell_params=sys.argv[2:]
    cell_params=[float(p) for p in cell_params]
    cell_params[3]*=math.pi/180
    cell_params[4]*=math.pi/180
    cell_params[5]*=math.pi/180

    #By default, grab the final set of coordinates
    step=0

elif len(sys.argv)==9:
    cell_params=sys.argv[2:8]
    cell_params=[float(p) for p in cell_params]
    cell_params[3]*=math.pi/180
    cell_params[4]*=math.pi/180
    cell_params[5]*=math.pi/180

    #Read second argument as step number
    step=int(sys.argv[8])
else:
    print('Unexpected number of arguments received. Expecting file_name + 6'+\
        'cell parameters + optional step number,'+\
        ' but received {}.File conversion will not take place.'.format(
            len(sys.argv)-1
            )
        )
    sys.exit()

#Grab the requested coordinates (by default, the final ones found in file).
atoms=readxyz(sys.argv[1],step)
#Convert to fractional coordinates
atoms=convert_coordinates(atoms,cell_params)

outf=str(sys.argv[1][:-4])+'.cif'
write_cif(outf,cell_params,atoms)

if step:
    print('\nCif file has been generated from xyz file {}'.format(sys.argv[1])+\
        ' using the coordinates from step {}.'.format(step))
else:
    print('\nCif file has been generated from xyz file {}'.format(sys.argv[1])+\
        ' using the coordinates from the last step.')
