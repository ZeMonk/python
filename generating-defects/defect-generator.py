###############################################################################
import numpy as np
import os
import sys
import networkx as nx
from datetime import date
###############################################################################

'''
Script for generating missing linker defects in Zr6 MOFs. Takes in some command
line arguments: an input .cif file (which should have P1 symmetry), a cap
coordinate .cif file (the cap attached to a bare Zr6 cluster), and some indices
for the linkers to be replaced.
The script may be run blank to generate a .cif file from which indices can
be chosen.

python3 defect-generator.py 

will generate skeleton.cif file.

python3 defect-generator.py file.cif coordinates.cif 1,2,3

will replace linkers indexed as 1, 2, and 3 in file.cif with the capping species
contained in coordinates.cif
'''

#==============================================================================
#==============================================================================
#raycov_dictionary is a dictionary containing the covalent radii of
#all elements. data taken from the CSD at 
#https://www.ccdc.cam.ac.uk/support-and-resources/
#ccdcresources/Elemental_Radii.xlsx
#at 18:10 on 27/11/2018
#Two further dummy radii, for linkers and clusters, are defined.
raycov_dict={
'H':0.23,'He':1.5,'Li':1.28,'Be':0.96,'B':0.83,'C':0.68,
'N':0.68,'O':0.68,'F':0.64,'Ne':1.5,'Na':1.66,'Mg':1.41,'Al':1.21,'Si':1.2,
'P':1.05,'S':1.02,'Cl':0.99,'Ar':1.51,'K':2.03,'Ca':1.76,'Sc':1.7,'Ti':1.6,
'Va':1.53,'Cr':1.39,'Mn':1.61,'Fe':1.52,'Co':1.26,'Ni':1.24,'Cu':1.32,
'Zn':1.22,'Ga':1.22,'Ge':1.17,'As':1.21,'Se':1.22,'Br':1.21,'Kr':1.5,'Rb':2.2,
'Sr':1.95,'Y':1.9,'Zr':1.75,'Nb':1.64,'Mo':1.54,'Tc':1.47,'Ru':1.46,'Rh':1.42,
'Pd':1.39,'Ag':1.45,'Cd':1.54,'In':1.42,'Sn':1.39,'Sb':1.39,'Te':1.47,'I':1.4,
'Xe':1.5,'Cs':2.44,'Ba':2.15,'La':2.07,'Ce':2.04,'Pr':2.03,'Nd':2.01,'Pm':1.99,
'Sm':1.98,'Eu':1.98,'Gd':1.96,'Tb':1.94,'Dy':1.92,'Ho':1.92,'Er':1.89,'Tm':1.9,
'Yb':1.87,'Lu':1.87,'Hf':1.75,'Ta':1.7,'W':1.62,'Re':1.51,'Os':1.44,'Ir':1.41,
'Pt':1.36,'Au':1.36,'Hg':1.32,'Tl':1.45,'Pb':1.46,'Bi':1.48,'Po':1.4,'At':1.21,
'Rn':1.5,'Fr':2.6,'Ra':2.21,'Ac':2.15,'Th':2.06,'Pa':2,'U':1.96,'Np':1.9,
'Pu':1.87,'Am':1.8,'Cm':1.69,'Bk':1.54,'Cf':1.83,'Es':1.5,'Fm':1.5,'Md':1.5,
'No':1.5,'Lr':1.5,'Rf':1.5,'Db':1.5,'Sg':1.5,'Bh':1.5,'Hs':1.5,
'Mt':1.5,'Ds':1.5,'Linker':0,'Cluster':0
}

#==============================================================================
#==============================================================================
#==============================================================================
#Class for storing data about atoms. Normally contains an integer label,
#the atomic symbol, the covalent radius, position coordinates, a boolean
#which describes whether the atom is real or a periodic image, and a set of
#bonded neighbours.

class Atom:

    def __init__(self,label,atom_type,position,image):

        self.L=label
        self.t=atom_type
        self.r=raycov_dict[self.t]
        self.p=np.asarray(position,dtype=np.float64)
        self.i=image
        self.n=set()

#==============================================================================
#==============================================================================
#Class for storing data and functions related to the unit cell. Initialised
#using 6 cell parameters (angles in degrees).

class Crystal:

    def __init__(self,cell_params,n_atoms):

        if len(cell_params)==6:
            self.cell_params=cell_params
        else:
            print('Cell parameters should contain 6 entries.')
            sys.exit()

        #Also store number of atoms and supercell status.
        self.n=int(n_atoms)
        self.supercell=(1,1,1)

#====
    #Method for building a supercell from the current unit cell. Adjusts
    #positions and cell parameters while adding new atoms.

    def make_supercell(self,supercell,atoms):

        #A supercell should really be generated only once, so attributes 
        #which would be needed to generate it a second time are not tracked.
        if self.supercell[0]*self.supercell[1]*self.supercell[2]!=1:
            print('Supercell should only be made once and will not be made.\n')
            return

        self.supercell=supercell
        #Edit the cell parameters.
        for i in range(3):
            self.cell_params[i]=self.cell_params[i]*self.supercell[i]

        #Variable N used to label the new atoms
        N=self.n
        #Add new atoms for each unit cell we are adding to the supercell.
        for i in range(self.supercell[0]):
            for j in range(self.supercell[1]):
                for k in range(self.supercell[2]):
                    #Skip the original unit cell.
                    if i+j+k==0:
                        continue
                    for x in range(self.n):
                        new_position=atoms[x].p+[i,j,k]
                        atoms.append(
                            Atom(x+N,atoms[x].t,new_position,False)
                            )
                        N+=1

        #Edit the coordinates of the original cell to reflect changes in cell
        #parameters.
        for q in range(self.n):
            for w in range(3):
                atoms[q].p[w]=atoms[q].p[w]/self.supercell[w]

#====
    #Builds and returns the metric tensor for the given unit cell from the 
    #cell parameters.

    def make_metric_tensor(self):

        a=self.cell_params[0]
        b=self.cell_params[1]
        c=self.cell_params[2]
        #Assumes the cell parameters are in degrees
        alpha=np.radians(self.cell_params[3])
        beta=np.radians(self.cell_params[4])
        gamma=np.radians(self.cell_params[5])

        #Standard crystallographic formula for an abritrary unit cell.
        m_tensor=np.array([
            [a*a,a*b*np.cos(gamma),a*c*np.cos(beta)],
            [a*b*np.cos(gamma),b*b,b*c*np.cos(alpha)],
            [a*c*np.cos(beta),b*c*np.cos(alpha),c*c]
            ]
            )

        return m_tensor

#====
    #Builds and returns the transformation matrix for going from fractional to
    #cartesian coordinates from the cell parameters. Optional flag allows for
    #the inverse (cartesian to fractional) to be made instead.

    def make_frac2cart(self,inverse=False):

        a=self.cell_params[0]
        b=self.cell_params[1]
        c=self.cell_params[2]
        #Assumes the cell parameters are in degrees
        alpha=np.radians(self.cell_params[3])
        beta=np.radians(self.cell_params[4])
        gamma=np.radians(self.cell_params[5])

        #Omega quantity is required to define the transformation matrix.
        omega=a*b*c*np.sqrt(
            1-np.power(np.cos(alpha),2) \
            -np.power(np.cos(beta),2) \
            -np.power(np.cos(gamma),2) \
            +2*np.cos(alpha)*np.cos(beta)*np.cos(gamma)
            )
        
        #Standard crystallographic formula for an arbitrary unit cell.
        t_matrix=np.array(
            [
            [a,b*np.cos(gamma),c*np.cos(beta)],
            [0,b*np.sin(gamma),c*(np.cos(alpha) \
                -np.cos(beta)*np.cos(gamma))/np.sin(gamma)],
            [0,0,(omega)/(a*b*np.sin(gamma))]
            ]
            )

        if inverse:
            return np.linalg.inv(t_matrix)
        else:
            return t_matrix

#====
    #Writes a cif file for the current structure given a set of atom objects.

    def makecif(self,outfile,atoms):

        with open(outfile,'w') as file:

            file.write('#File created on {}\n\n'.format(date.today()))

            #Pre-position block.
            file.write('data_image0\n'\
            +'_symmetry_space_group_name_H-M    "P 1"\n'\
            +'_symmetry_int_tables_number       1\n\n'\
            +'loop_\n'\
            +' _symmetry_equiv_pos_as_xyz\n'\
            +" 'x, y, z'\n\n")

            #Cell parameter block.
            file.write(
            '_cell_length_a {:.4f}\n'.format(self.cell_params[0])\
            +'_cell_length_b {:.4f}\n'.format(self.cell_params[1])\
            +'_cell_length_c {:.4f}\n'.format(self.cell_params[2])\
            +'_cell_angle_alpha {:.4f}\n'.format(self.cell_params[3])\
            +'_cell_angle_beta {:.4f}\n'.format(self.cell_params[4])\
            +'_cell_angle_gamma {:.4f}\n\n'.format(self.cell_params[5]))

            #Position block.
            file.write('loop_\n'\
            +' _atom_site_label\n'\
            +' _atom_site_type_symbol\n'\
            +' _atom_site_fract_x\n'\
            +' _atom_site_fract_y\n'\
            +' _atom_site_fract_z\n')

            for n,a in enumerate(atoms):
                file.write('{} {} {:.5f} {:.5f} {:.5f}\n'.format(a.t+str(n),a.t
                        ,a.p[0],a.p[1],a.p[2]))

        print('New .cif file {} was successfully generated.'.format(outfile))

#==============================================================================
#==============================================================================
#Function to read lines until the target string is found. Line is returned
#if found, else False returned if EOF is reached.

def readuntil(fin,target):

    while True:
        line=fin.readline()
        if target in line:
            return line
        #Evaluates to False for empty string which EOF returns
        elif not line:
            return False

#==============================================================================
#==============================================================================
#Function for reading a cif file and pulling cell and atom information. Works
#only for cifs of P1 symmetry, though this functionality could be added using
#pymatgen.

def cifreader(filename):

    with open(filename,'r') as fin:

        #Read cif file until cell parameters
        a=float(readuntil(fin,'_cell_length_a').split()[-1])
        b=float(fin.readline().split()[-1])
        c=float(fin.readline().split()[-1])
        alpha=float(fin.readline().split()[-1])
        beta=float(fin.readline().split()[-1])
        gamma=float(fin.readline().split()[-1])

        #Read until coordinate and atom type headers and collect them
        headers=[]
        headerstring='_atom_site_'
        headers.append(readuntil(fin,headerstring).strip())
        while True:
            line=fin.readline()
            if headerstring in line:
                headers.append(line.strip())
            else:
                break

        #Once headers have been collected, move back to start of atoms block
        fin.seek(0)
        readuntil(fin,headers[-1])

        atoms=[]
        i=0
        while True:
            line=fin.readline()
            if line:
                try:
                    line=line.split()
                    #Using the index of the needed headers in headers array
                    #ensures that, even if other fields are present, only the
                    #atom type, x, y, and z coordinates will be pulled.
                    atoms.append(Atom(
                        i,line[headers.index('_atom_site_type_symbol')],
                        [line[headers.index('_atom_site_fract_x')],
                        line[headers.index('_atom_site_fract_y')],
                        line[headers.index('_atom_site_fract_z')]],
                        False
                        ))
                    i+=1
                except IndexError:
                    break
            else:
                break

    crystal=Crystal(
        np.array([a,b,c,alpha,beta,gamma],dtype=np.float64),len(atoms)
            )

    return atoms,crystal

#==============================================================================
#==============================================================================
#This function finds all atoms too close to the far edge of the unit cell
#(as controlled by cutoff) and appends their periodic images at the 
#appropriate position (i.e. current position -1) to atoms list to mimic 
#translational symmetry.

def t_symm(atoms,crystal,cutoff=3):

    for tt in range(0,len(atoms)):
        #first, make sure all atoms are inside the unit cell
        for s in range(3):
            while atoms[tt].p[s]<0.0:
                atoms[tt].p[s]=atoms[tt].p[s]+1.0
            while atoms[tt].p[s]>1.0:
                atoms[tt].p[s]=atoms[tt].p[s]-1.0

        #track adjacency to a cell boundary with 3 flags (x,y,z)
        edge_flags=[False,False,False]
        for ttt in range(0,3):
            if (1-atoms[tt].p[ttt])*crystal.cell_params[ttt]<cutoff:
                edge_flags[ttt]=True
                #for each far edge that an atom is near, set flag to true

        #for each combination of True flags, add an image atom in the
        #corresponding position
        #e.g. [False,True,True] will add atoms at p-[0,0,1], p-[0,1,0], and
        #p-[0,1,1]
        for i in range(int(edge_flags[0])+1):
            for j in range(int(edge_flags[1])+1):
                for k in range(int(edge_flags[2])+1):
                    if i+j+k!=0:
                        vector=[i,j,k]
                        atoms.append(
                            Atom(
                                atoms[tt].L,atoms[tt].t,
                                atoms[tt].p-vector,True
                                )
                            )
                        atoms[tt].edge=edge_flags
                        #copy all attributes of the original atom, except
                        #the image label which is True

    return atoms

#==============================================================================
#==============================================================================
#Function for finding bonds between the given set of atoms. Requires a crystal
#object containing some cell parameters in order to generate the metric
#tensor, which is used to calculate distances valid for any crystal system.
#Default tolerance of 0.4 angstrom is a standard value.

def build_bonds(atoms,crystal,tol=0.4):

    #Metric tensor object is a 3x3 matrix.
    m_tensor=crystal.make_metric_tensor()

    for i in range(0,len(atoms)):
        for j in range(i+1,len(atoms)):
            vector=atoms[i].p-atoms[j].p
            #Generalised Pythagoras formula works as long as m_tensor is correct
            mag=np.sqrt(vector@m_tensor@vector)
            if atoms[i].r+atoms[j].r-tol < mag < atoms[i].r+atoms[j].r+tol:
                #Due to method used for translational symmetry, image atoms are
                #also present. Hence add the labels of the ith and jth atoms
                #instead (since these were copied from their counterparts).
                atoms[atoms[i].L].n.add(atoms[j].L)
                atoms[atoms[j].L].n.add(atoms[i].L)

#==============================================================================
#==============================================================================
#Function for finding the centroid of a molecular fragment. Takes into account
#possible disconnects across unit cell boundaries

def get_centroid(atoms,fragment):

    #Extract the atom positions to avoid editing atom positions.
    positions=[atoms[i].p for i in fragment]
    #Set one of the atoms as the initial geometric centre to compare to.
    gc=positions[0]


    #For each atom, if it is farther away than 0.75 for any coordinate from
    #the centroid, edit its position to place it in a closer image position
    for p in positions:
        for i in range(3):
            if gc[i]-p[i]>0.75:
                p[i]+=1
            elif gc[i]-p[i]<-0.75:
                p[i]-=1

    #Geometric centroid simply the average of all corrected coordinates.
    gc=np.mean(positions,axis=0)

    return gc

#==============================================================================
#==============================================================================
#Function for generating a local Cartesian basis set given a pair of oxygens
#corresponding to linker-SBU bonds

def basis_generator(atoms,ends,crystal):

    metals=['Zr']
    vectors=[]

    #Find O-Zr vectors
    for oxygen in ends:
        for j in atoms[oxygen].n:
            if atoms[j].t in metals:
                vectors.append(atoms[oxygen].p-atoms[j].p)

    #Check that vector isn't cut across a unit cell boundary
    for v in vectors:
        for i in range(3):
            if v[i]>0.5:
                v[i]-=1
            elif v[i]<-0.5:
                v[i]+=1

    #Transform vectors to Cartesian
    t_matrix=crystal.make_frac2cart()
    vectors=[t_matrix@v for v in vectors]

    #Basis vector 1 defined using mean of O-Zr bonds, is normalised.
    basis1=np.mean([vectors],axis=1)[0]
    basis1/=np.linalg.norm(basis1)

    #Basis vector 2 defined using cross product of O-Zr bonds, is normalised.
    basis2=np.cross(vectors[0],vectors[1])
    basis2/=np.linalg.norm(basis2)

    #Basis vector 3 defined using cross product of other 2 basis vectors and is
    #already normalised.
    basis3=np.cross(basis1,basis2)
    basis_set=np.array([basis1,basis2,basis3],dtype=np.float64)

    return basis_set

#==============================================================================
#==============================================================================
#Function for generating capping species coordinates given a coordinate system
#and an origin corresponding to the specific part of the linker being replaced.

def capping_generator(filename,mof_basis,mof_crystal,origin):

    #Read capping species coordinates and separate into linker+SBU
    cap_atoms,cap_crystal=cifreader(filename)
    build_bonds(cap_atoms,cap_crystal,tol=0.4)
    cap,cluster=find_fragments(cap_atoms)

    #Find the capping species' end oxygens in order to generate a local basis
    #set.
    ends=[]
    for a in cap_atoms:
        ntypes=set([cap_atoms[i].t for i in a.n])
        if ntypes.intersection(set(['O','C'])) and len(a.n)==3:
            for j in a.n:
                if cap_atoms[j].t=='O':
                        ends.append(j)

    #Generate cap basis set and fractional-to-cartesian matrix
    cap_basis=basis_generator(cap_atoms,ends,cap_crystal)
    cap_t_matrix=cap_crystal.make_frac2cart()

    #Generate a set of transform matrices from the MOF coordinate system.
    inverse=mof_basis.transpose()
    t_matrix=mof_crystal.make_frac2cart(inverse=True)


    cap_origin=np.mean([cap_atoms[k].p for k in ends],axis=0)

    '''In order:
    1. set the capping species' origin to be the centroid of the 2 end oxygens
    2. transform cap coordinates to Cartesian, then to local basis
    3. transform edited cap coordinates to MOF basis, then to MOF fractional,
    then add the origin of the linker end in the MOF'''
    for a in cap_atoms:
        a.p-=cap_origin
        a.p=cap_basis@cap_t_matrix@a.p
        a.p=t_matrix@inverse@a.p+origin

    #Keep only the atoms corresponding to the cap atoms
    cluster_atoms=[cap_atoms[l] for l in cluster[0].members]
    cap_atoms=[a for a in cap_atoms if a not in cluster_atoms]

    return cap_atoms

#==============================================================================
#==============================================================================
#Function for generating missing linker defects. Take a set of target indices
#from the array of identified linkers, as well as an external coordinate file
#which should contain the cap attached to the SBU.

def defect_maker(targets,total_array,atoms,crystal,coordfile):

    #Start with any possible clusters
    targets.sort(reverse=True)
    linker_ends=[]
    removal=[]
    removed_linkers=[]


    #Identify all sets of terminating oxygens
    for fragment in targets:
        ends=[]

        #For a cluster, identify all terminating oxygens on the far sides of
        #connected linkers
        if total_array[fragment].t=='Cluster':
            m_tensor=crystal.make_metric_tensor()
            for neighbour in total_array[fragment].n:
                removed_linkers.append(neighbour)
                linker_atoms=[atoms[i] for i in total_array[neighbour].members]
                for at in linker_atoms:
                    #Track which atoms to remove
                    removal.append(at)
                    #Check how far the target oxygens are and keep only the far
                    #ones
                    if at.t=='O':
                        vector=at.p-total_array[fragment].p
                        #Shortest image distance
                        for i in range(3):
                            if vector[i]>0.5:
                                vector[i]-=1
                            elif vector[i]<-0.5:
                                vector[i]+=1
                        #Oxygen on the far end of a linker should be around
                        #9-10 angstrom away, here the distance squared is
                        #checked. Would need adjustment for other MOFs
                        if vector@m_tensor@vector>80:
                            ends.append(at.L)
                linker_ends.append(ends)
                ends=[]

        for member in total_array[fragment].members:
            #Track all atoms to be removed
            removal.append(atoms[member])

            #Skip any linkers already removed
            if fragment in removed_linkers:
                continue
            elif total_array[fragment].t=='Linker':
            #Simple criteria for finding the ending oxygens
                ntypes=set([atoms[i].t for i in atoms[member].n])
                if ntypes.intersection(set(['O','C'])) and \
                len(atoms[member].n)==3:
                    for j in atoms[member].n:
                        if atoms[j].t=='O':
                            ends.append(j)

                    #Append ends only if it contains items; empty list evaluates
                    #to False.
                    if ends:
                        linker_ends.append(ends)
                        ends=[]

    #Generate a set of coordinates for the capping species for each end of each
    #linker that was given.
    new_atoms=[]
    for ends in linker_ends:
        origin=np.mean([atoms[k].p for k in ends],axis=0)
        mof_basis=basis_generator(atoms,ends,crystal)
        new_atoms.append(capping_generator(
                coordfile,mof_basis,
                crystal,origin
                ))
    #Remove original linker atoms
    atoms=[a for a in atoms if a not in removal]
    for new in new_atoms:
        for k in new:
            atoms.append(k)

    return atoms

#==============================================================================
#==============================================================================
#Function to isolate building blocks in the given MOF. Builds a graph object
#by initialising every bond as an edge. Then removes edges corresponding to
#linker-SBU bonds, and finally finds connected subgraphs.

def find_fragments(atoms):

    metals=set(['Zr'])
    removebonds=[]
    
    #For Zr6 MOFs, this set of criteria is sufficient for identifying the
    #bonds which should be removed from the graph.
    for a in atoms:
        if len(a.n)==2 and a.t=='O':
            for b in a.n:
                if atoms[b].t in metals:
                    removebonds.append((a.L,b))

    #Initialise the graph object, build the edges, then remove SBU-linker
    #bonds
    graph=nx.Graph()
    for a in atoms:
        for b in a.n:
            graph.add_edge(a.L,b)
    for rb in removebonds:
        graph.remove_edge(rb[0],rb[1])

    #Find molecular fragments by identifying connected subgraphs
    subgraphs=[]
    for g in nx.connected_components(graph):
        subgraphs.append(list(graph.subgraph(g).nodes))

    #Each subgraph is stored as an atom object with a label and position.
    clusters=[]
    linkers=[]
    for g in range(len(subgraphs)):
        fragment=[i for i in subgraphs[g]]
        atomtypes=set([atoms[j].t for j in fragment])
        #A set with elements in it evaluates to True, hence if there is an
        #overlap with the metals set, must be an SBU
        if atomtypes.intersection(metals):
            cluster=Atom(0,'Cluster',get_centroid(atoms,fragment),False)
            cluster.members=fragment #Members attribute not present by default
            #for atoms but useful to add to molecular fragments.
            clusters.append(cluster)
        else:
            linker=Atom(0,'Linker',get_centroid(atoms,fragment),False)
            linker.members=fragment
            linkers.append(linker)

    #Relabel fragments for convenience.
    for q in range(len(clusters)):
        clusters[q].L=q
    for q in range(len(linkers)):
        linkers[q].L=q

    return linkers,clusters

#==============================================================================
#==============================================================================
#Wrapper calls the necessary functions sequentially.

def _wrapper(filename,coordfile,indices,skeleton=False):

    #Read the input data, apply PBCs.
    atoms,crystal=cifreader(filename)
    t_symm(atoms,crystal)

    #Build bonds and then remove image atoms.
    build_bonds(atoms,crystal,tol=0.4)
    atoms=[gg for gg in atoms if gg.i==False]

    #Find linkers and clusters using graph algorithm
    linker_array,cluster_array=find_fragments(atoms)
    total_array=linker_array+cluster_array

    #Relabel fragments in total_array, find which linkers are connected to
    #each cluster
    for i,fragment in enumerate(total_array):
        fragment.L=i
    t_symm(total_array,crystal,cutoff=12)
    build_bonds(total_array,crystal,tol=12)
    total_array=[gg for gg in total_array if gg.i==False]

    #Make skeleton cif file if requested, which contains only linkers and
    #clusters as dummy Xe and La atoms.
    if skeleton:
        for c in cluster_array:
            c.t='Xe'
        for l in linker_array:
            l.t='La'
        crystal.makecif('skeleton.cif',total_array)
        sys.exit()


    foutname=filename[:-4]+'-{}-defects.cif'.format(len(indices))
    new_atoms=defect_maker(indices,total_array,atoms,crystal,coordfile)
    crystal.makecif(foutname,new_atoms)

#==============================================================================
#==============================================================================
#Run the code via the wrapper function. Receive and format some command-line
#arguments for file name, cap coordinates, and defect indices.

if __name__=="__main__":

    if 1 < len(sys.argv)<4:
        print('No coordinate file or indices given; skeleton.cif file will'+\
            ' be generated instead \n')
        filename=sys.argv[1]
        coordfile='dummy'
        indices=[]
        _wrapper(sys.argv[1],coordfile,indices,True)

    elif len(sys.argv)>3:
        filename=sys.argv[1]
        coordfile=sys.argv[2]
        indices=[int(i) for i in sys.argv[3].split(',')]

        print('Arguments received for input  file: {}'.format(filename))
        print('Argument received for coordinate file: {}'.format(coordfile))
        print('Arguments received for linker indices: {}\n'.format(sys.argv[3]))

        _wrapper(sys.argv[1],sys.argv[2],indices,False)
    else:
        print('Insufficient arguments received; no defects will be generated.')
        sys.exit()
