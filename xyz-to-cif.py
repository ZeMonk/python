import sys
import math
from datetime import date

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
            numsteps=i/natoms

        #Return to start of file.
        fin.seek(0)

        #Read ahead to requested coordinate block and extract it
        if step:
            for j in range(int(step-1)*(natoms+2)):
                fin.readline()
            atoms=coords(fin)
        else:
            for j in range(int(numsteps-1)*(natoms+2)):
                fin.readline()
            atoms=coords(fin)

        return atoms

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

def convert_coordinates(atoms,cell_params):
    
    t_matrix=cart_to_frac_matrix(cell_params)

    for a in atoms:
        a[1],a[2],a[3]=matrix_multiply(t_matrix,a[1:4])

        for i in range(3):
            while a[1+i]<0.0:
                a[1+i]+=1
            while a[1+i]>1.0:
                a[1+i]-=1

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
makecif(outf,cell_params,atoms)

if step:
    print('\nCif file has been generated from xyz file {}'.format(sys.argv[1])+\
        ' using the coordinates from step {}.'.format(step))
else:
    print('\nCif file has been generated from xyz file {}'.format(sys.argv[1])+\
        ' using the coordinates from the last step.')