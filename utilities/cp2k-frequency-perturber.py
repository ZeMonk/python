import numpy as np
import sys

'''
Python script which takes a CP2K VIBRATIONAL_ANALYSIS output file (at MEDIUM
print level) to perturb the structure using some chosen eigenvectors. Useful
if some undesired negative frequencies have been found.

python3 cp2k-frequency-perturber.py inputfile scalefactor index1,index2,...

indices start at 0.
'''

###############################################################################

def read_until(file,target):

    while True:
        line=fin.readline()
        if target in line:
            return line
        elif not line:
            return False

###############################################################################

def read_coordinates(file):

    read_until(file,'TOTAL NUMBERS AND MAXIMUM NUMBERS')
    file.readline()
    file.readline()
    natoms=int(file.readline().strip().split()[-1])

    read_until(file,'Atom Kind Element')

    coords=np.zeros([natoms,3],dtype=np.float64)
    species=[]
    for i in range(natoms):
        line=file.readline().split()
        coords[i,:]=line[4:7]
        species.append(str(line[2]))

    return coords,species,natoms

###############################################################################
#Reads normal modes from CP2K-formatted output file for a vibrational analysis
#carried out at medium print level.

def read_nm(file,natoms):

    normal_modes=[]
    frequencies=[]

    #Read all normal modes
    while True:
        line=read_until(file,'VIB|Frequency')
        #Once all normal modes read, read_until evaluates to False upon EOF
        if not line:
            return normal_modes,frequencies
        freqs=line.strip().split()[2:]
        for f in freqs:
            frequencies.append(float(f))
        #Skip some lines
        for i in range(3):
            file.readline()

        nm=np.zeros([natoms,3,3],dtype=np.float64)

        line=file.readline().strip().split()
        while line:
            index=int(line[0])-1
            #Use try-except clause to avoid errors when the number of normal
            #modes after translation and rotation have been removed is not a
            #multiple of 3
            try:
                nm[index,0,:]=line[2:5]
                nm[index,1,:]=line[5:8]
                nm[index,2,:]=line[8:]
            except IndexError:
                pass
            line=file.readline().strip().split()
        for i in range(3):
            normal_modes.append(nm[:,i,:])


###############################################################################

def xyz_writer(foutname,coords,species):

    with open(foutname,'w') as fout:

        fout.write('{}\n\n'.format(len(species)))
        for i in range(len(species)):
            fout.write('{} {} {} {}\n'.format(species[i],*coords[i]))

###############################################################################

filename=sys.argv[1]
scaling_factor=float(sys.argv[2])
indices=[int(i) for i in sys.argv[3].split(',')]

with open(filename,'r') as fin:

    coords,species,natoms=read_coordinates(fin)
    normal_modes,frequencies=read_nm(fin,natoms)
    
    for index in indices:
        for i in range(len(coords)):
            for j in range(3):
                coords[i,j]+=scaling_factor*normal_modes[index][i][j]
    xyz_writer('perturbed.xyz',coords,species)
