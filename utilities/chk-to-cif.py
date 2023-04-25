import sys
import math
from datetime import date

'''
Short utility script to convert a .chk output file from a Yaff job to cif
format. Keeps the information about charges and force field atom type labels
in the output cif.

Use with python3 some/example/path/targetfile
'''

#==============================================================================
#==============================================================================
#Dictionary of atomic numbers and their string labels, along with the reverse.

a_number_dict={
'H':1,'He':2,'Li':3,'Be':4,'B':5,'C':6,'N':7,'O':8,'F':9,'Ne':10,
'Na':11,'Mg':12,'Al':13,'Si':14,'P':15,'S':16,'Cl':17,'Ar':18,'K':19,'Ca':20,
'Sc':21,'Ti':22,'V':23,'Cr':24,'Mn':25,'Fe':26,'Co':27,'Ni':28,'Cu':29,'Zn':30,
'Ga':31,'Ge':32,'As':33,'Se':34,'Br':35,'Kr':36,'Rb':37,'Sr':38,'Y':39,'Zr':40,
'Nb':41,'Mo':42,'Tc':43,'Ru':44,'Rh':45,'Pd':46,'Ag':47,'Cd':48,'In':49,
'Sn':50,'Sb':51,'Te':52,'I':53,'Xe':54,'Cs':55,'Ba':56,'La':57,'Ce':58,'Pr':59,
'Nd':60,'Pm':61,'Sm':62,'Eu':63,'Gd':64,'Tb':65,'Dy':66,'Ho':67,'Er':68,
'Tm':69,'Yb':70,'Lu':71,'Hf':72,'Ta':73,'W':74,'Re':75,'Os':76,'Ir':77,'Pt':78,
'Au':79,'Hg':80,'Tl':81,'Pb':82,'Bi':83,'Po':84,'At':85,'Rn':86,'Fr':87,
'Ra':88,'Ac':89,'Th':90,'Pa':91,'U':92,'Np':93,'Pu':94,'Am':95,'Cm':96,'Bk':97,
'Cf':98,'Es':99,'Fm':100,'Md':101,'No':102,'Lr':103,'Rf':104,'Db':105,'Sg':106,
'Bh':107,'Hs':108,'Mt':109,'Ds':110,'Rg':111,'Cn':112,'Nh':113,'Fl':114,
'Mc':115,'Lv':116,'Ts':117,'Og':118
}
rdict=dict([(v,k) for k,v in a_number_dict.items()])

#==============================================================================
#==============================================================================
#Function to read lines until the target string is found. Line is returned
#if found, else False returned if EOF is reached.

def readuntil(fin,target):

    while True:
        line=fin.readline()
        if target in line:
            return line
        elif not line:
            return False

#==============================================================================
#==============================================================================
#Function which collects a data block from a .chk file and returns as a list.

def collect_simple(fin,n):

    data=[]

    #.chk files store 4 data per line; this syntax ensures the correct number
    #of lines is read
    for i in range(int(n/4)+(n%4>0)):
        for j in fin.readline().split():
            data.append(j)

    return data

#==============================================================================
#==============================================================================
#Function which reshapes a 1 by x list into a y by 3 list, needed for positions.

def reshape_x3(data):

    reshaped_data=[]

    #Data should always be a multiple of 3 by design.
    for i in range(int(len(data)/3)):
        reshaped_data.append(
            [data[i*3],data[i*3+1],data[i*3+2]])

    return reshaped_data

#==============================================================================
#==============================================================================
#Function which reads a .chk file and collects atom numbers, positions, charges,
#force field types, and cell vectors, formats and returns them.
#Returning to the start of the file, in principle, should avoid any problems
#that might arise from data blocks not being in the expected order.

def collect_data(filename):

    with open(filename,'r') as fin:

        #Move to charges in file and collect number of atoms and charges
        natoms=int(readuntil(fin,'charges').split()[-1])
        charges=collect_simple(fin,natoms)
        charges=[float(i) for i in charges]

        #Collect force field atom type ids
        fin.seek(0)
        readuntil(fin,'ffatype_ids')
        ffatype_ids=collect_simple(fin,natoms)
        ffatype_ids=[int(j) for j in ffatype_ids]

        #Collect force field atom types, generate conversion dictionary, and
        #turn them into their corresponding string labels
        fin.seek(0)
        types=int(readuntil(fin,'ffatypes').split()[-1])
        ffatypes=collect_simple(fin,types)
        ffadict=dict([(k,str(v)) for k,v in enumerate(ffatypes)])
        labels=[ffadict[l] for l in ffatype_ids]

        #Collect atomic numbers
        fin.seek(0)
        readuntil(fin,'numbers')
        anumbers=collect_simple(fin,natoms)
        anumbers=[int(m) for m in anumbers]

        #Collect atom coordinates and reshape list
        fin.seek(0)
        readuntil(fin,'pos')
        raw_coords=collect_simple(fin,natoms*3)
        #Convert coordinates from bohr to angstrom
        b2a=0.529177
        raw_coords=[float(c)*b2a for c in raw_coords]
        coords=reshape_x3(raw_coords)

        #Collect cell vectors and reshape list
        fin.seek(0)
        readuntil(fin,'rvecs')
        raw_vecs=collect_simple(fin,9)
        #Convert cell vectors from bohr to angstrom
        raw_vecs=[float(v)*b2a for v in raw_vecs]

    return anumbers,coords,raw_vecs,labels,charges

#==============================================================================
#==============================================================================
#Short function which takes a 1x9 list of cell vector elements and converts
#them to cell parameters. Angles are in radians.

def vecs_to_params(vecs):

    a=math.sqrt(vecs[0]*vecs[0]+vecs[1]*vecs[1]+vecs[2]*vecs[2])
    b=math.sqrt(vecs[3]*vecs[3]+vecs[4]*vecs[4]+vecs[5]*vecs[5])
    c=math.sqrt(vecs[6]*vecs[6]+vecs[7]*vecs[7]+vecs[8]*vecs[8])

    alpha=math.acos((vecs[3]*vecs[6]+vecs[4]*vecs[7]+vecs[5]*vecs[8])/(b*c))
    beta=math.acos((vecs[0]*vecs[6]+vecs[1]*vecs[7]+vecs[2]*vecs[8])/(b*c))
    gamma=math.acos((vecs[0]*vecs[3]+vecs[1]*vecs[4]+vecs[2]*vecs[5])/(b*c))

    return [a,b,c,alpha,beta,gamma]

#==============================================================================
#==============================================================================
#If the coordinates printed in the CP2K restart file are Cartesian, then they
#need to be converted to fractional using the corresponding transformation
#matrix. This function generates that matrix from the cell parameters.

def cart_to_frac_matrix(cell_params):

    a,b,c,alpha,beta,gamma=cell_params

    #Standard crystallographic transformation matrix involves 2 quantities,
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

def convert_coordinates(atoms,cell_params):
    
    t_matrix=cart_to_frac_matrix(cell_params)

    for i in range(len(atoms)):
        atoms[i]=matrix_multiply(t_matrix,atoms[i])

        for j in range(3):
            while atoms[i][j]<0.0:
                atoms[i][j]+=1

            while atoms[i][j]>1.0:
                atoms[i][j]-=1

    return atoms

#==============================================================================
#==============================================================================
#Function for writing out a .cif file. The space group information is printed
#at the top of the file, regardless of whether symmetry is requested.

def makecif(foutname,cell_params,labels,numbers,coords,charges):

    with open(foutname,'w') as fout:
        fout.write('#This file was created on {}\n'.format(date.today()))

        fout.write('data_image0\n')
        fout.write('_symmetry_space_group_name_H-M P 1\n')
        fout.write('_space_group_IT_number 1\n\n')

        fout.write('_cell_length_a {:.4f}\n'.format(cell_params[0]))
        fout.write('_cell_length_b {:.4f}\n'.format(cell_params[1]))
        fout.write('_cell_length_c {:.4f}\n'.format(cell_params[2]))
        fout.write('_cell_angle_alpha {:.4f}\n'.format(
            cell_params[3]*180/math.pi)
        )
        fout.write('_cell_angle_beta {:.4f}\n'.format(
            cell_params[4]*180/math.pi)
        )
        fout.write('_cell_angle_gamma {:.4f}\n\n'.format(
            cell_params[5]*180/math.pi)
        )

        fout.write('loop_\n')
        fout.write('_atom_site_label\n')
        fout.write('_atom_site_type_symbol\n')
        fout.write('_atom_site_fract_x\n')
        fout.write('_atom_site_fract_y\n')
        fout.write('_atom_site_fract_z\n')
        fout.write('_atom_site_charge\n')

        for i,a in enumerate(coords):
            fout.write('{} {} {:.5f} {:.5f} {:.5f} {:.5f}\n'.format(
                        labels[i],rdict[numbers[i]],a[0],
                        a[1],a[2],charges[i])
                    )

#==============================================================================
#==============================================================================

if len(sys.argv)>2:
    print('''Unexpected command line arguments given. Only the first argument
(a path+file name) will be processed, the rest will be ignored.\n''')

filename=sys.argv[1]

#Collect and format data.
numbers,coords,vecs,labels,charges=collect_data(filename)
cell_params=vecs_to_params(vecs)
frac_coords=convert_coordinates(coords,cell_params)

foutname=filename[:-4]+'.cif'
makecif(foutname,cell_params,labels,numbers,frac_coords,charges)

print('.cif file {} was generated.'.format(foutname))





