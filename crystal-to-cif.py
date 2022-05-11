import sys
from datetime import date

'''
This short script converts the output from a CRYSTAL structural optimisation
to cif format. Given a valid output file from a converged calculation, this
will generate 2 cif files, one with the same space group symmetry containing
only the asymmetric unit, and another with P1 symmetry containing all atoms.

Run using python3 crystal-to-cif.py some/example/path/filename
'''

#==============================================================================
#==============================================================================
#Dictionary of space group names and numbers, along with the reverse.

spacegroups={
'P 1':1,'P -1':2,'P 2':3,'P 21':4,'C 2':5,'P M':6,'P C':7,'C M':8,'C C':9,
'P 2/M':10,'P 21/M':11,'C 2/M':12,'P 2/C':13,'P 21/C':14,'C 2/C':15,
'P 2 2 2':16,'P 2 2 21':17,'P 21 21 2':18,'P 21 21 21':19,'C 2 2 21':20,
'C 2 2 2':21,'F 2 2 2':22,'I 2 2 2':23,'I 21 21 21':24,'P M M 2':25,
'P M C 21':26,'P C C 2':27,'P M A 2':28,'P C A 21':29,'P N C 2':30,
'P M N 21':31,'P B A 2':32,'P N A 21':33,'P N N 2':34,'C M M 2':35,
'C M C 21':36,'C C C 2':37,'A M M 2':38,'A B M 2':39,'A M A 2':40,'A B A 2':41,
'F M M 2':42,'F D D 2':43,'I M M 2':44,'I B A 2':45,'I M A 2':46,'P M M M':47,
'P N N N':48,'P C C M':49,'P B A N':50,'P M M A':51,'P N N A':52,'P M N A':53,
'P C C A':54,'P B A M':55,'P C C N':56,'P B C M':57,'P N N M':58,'P M M N':59,
'P B C N':60,'P B C A':61,'P N M A':62,'C M C M':63,'C M C A':64,'C M M M':65,
'C C C M':66,'C M M A':67,'C C C A':68,'F M M M':69,'F D D D':70,'I M M M':71,
'I B A M':72,'I B C A':73,'I M M A':74,'P 4':75,'P 41':76,'P 42':77,'P 43':78,
'I 4':79,'I 41':80,'P -4':81,'I -4':82,'P 4/M':83,'P 42/M':84,'P 4/N':85,
'P 42/N':86,'I 4/M':87,'I 41/A':88,'P 4 2 2':89,'P 4 21 2':90,'P 41 2 2':91,
'P 41 21 2':92,'P 42 2 2':93,'P 42 21 2':94,'P 43 2 2':95,'P 43 21 2':96,
'I 4 2 2':97,'I 41 2 2':98,'P 4 M M':99,'P 4 B M':100,'P 42 C M':101,
'P 42 N M':102,'P 4 C C':103,'P 4 N C':104,'P 42 M C':105,'P 42 B C':106,
'I 4 M M':107,'I 4 C M':108,'I 41 M D':109,'I 41 C D':110,'P -4 2 M':111,
'P -4 2 C':112,'P -4 21 M':113,'P -4 21 C':114,'P -4 M 2':115,'P -4 C 2':116,
'P -4 B 2':117,'P -4 N 2':118,'I -4 M 2':119,'I -4 C 2':120,'I -4 2 M':121,
'I -4 2 D':122,'P 4/M M M':123,'P 4/M C C':124,'P 4/N B M':125,'P 4/N N C':126,
'P 4/M B M':127,'P 4/M N C':128,'P 4/N M M':129,'P 4/N C C':130,
'P 42/M M C':131,'P 42/M C M':132,'P 42/N B C':133,'P 42/N N M':134,
'P 42/M B C':135,'P 42/M N M':136,'P 42/N M C':137,'P 42/N C M':138,
'I 4/M M M':139,'I 4/M C M':140,'I 41/A M D':141,'I 41/A C D':142,'P 3':143,
'P 31':144,'P 32':145,'R 3':146,'P -3':147,'R -3':148,'P 3 1 2':149,
'P 3 2 1':150,'P 31 1 2':151,'P 31 2 1':152,'P 32 1 2':153,'P 32 2 1':154,
'R 3 2':155,'P 3 M 1':156,'P 3 1 M':157,'P 3 C 1':158,'P 3 1 C':159,
'R 3 M':160,'R 3 C':161,'P -3 1 M':162,'P -3 1 C':163,'P -3 M 1':164,
'P -3 C 1':165,'R -3 M':166,'R -3 C':167,'P 6':168,'P 61':169,'P 65':170,
'P 62':171,'P 64':172,'P 63':173,'P -6':174,'P 6/M':175,'P 63/M':176,
'P 6 2 2':177,'P 61 2 2':178,'P 65 2 2':179,'P 62 2 2':180,'P 64 2 2':181,
'P 63 2 2':182,'P 6 M M':183,'P 6 C C':184,'P 63 C M':185,'P 63 M C':186,
'P -6 M 2':187,'P -6 C 2':188,'P -6 2 M':189,'P -6 2 C':190,'P 6/M M M':191,
'P 6/M C C':192,'P 63/M C M':193,'P 63/M M C':194,'P 2 3':195,'F 2 3':196,
'I 2 3':197,'P 21 3':198,'I 21 3':199,'P M 3':200,'P N 3':201,'F M 3':202,
'F D 3':203,'I M 3':204,'P A 3':205,'I A 3':206,'P 4 3 2':207,'P 42 3 2':208,
'F 4 3 2':209,'F 41 3 2':210,'I 4 3 2':211,'P 43 3 2':212,'P 41 3 2':213,
'I 41 3 2':214,'P -4 3 M':215,'F -4 3 M':216,'I -4 3 M':217,'P -4 3 N':218,
'F -4 3 C':219,'I -4 3 D':220,'P M 3 M':221,'P N 3 N':222,'P M 3 N':223,
'P N 3 M':224,'F M 3 M':225,'F M 3 C':226,'F D 3 M':227,'F D 3 C':228,
'I M 3 M':229,'I A 3 D':230,
}
rspacegroups=dict([(v,k) for k,v in spacegroups.items()])

#==============================================================================
#==============================================================================
#Function to read lines until the target string is found. Line is returned
#if found, else False returned if EOF is reached.

def readuntil(file,target):

    while True:
        line=fin.readline()
        if target in line:
            return line
        elif not line:
            return False

#==============================================================================
#==============================================================================
#Function for collecting atom types and coordinates from CRYSTAL output.
#Assumes that current file position is right above first coordinate line.

def coords(fin,natoms):

    atoms=[]

    for i in range(natoms):
        l=fin.readline().split()

        #CRYSTAL prints atom types in all-caps. If atom type symbol has more
        #than 2 letters,convert those beyond the first to lower case.
        if len(l[3])>1:
            atomtype=l[3][0]+l[3][1:].lower()
        else:
            atomtype=l[3]

        atoms.append([l[1],atomtype,float(l[4]),float(l[5]),float(l[6])])

    return atoms

#==============================================================================
#==============================================================================
#Function for writing out a .cif file. The space group information is printed
#at the top of the file, regardless of whether symmetry is requested.

def makecif(cell_params,atoms,spacegroup,foutname,symmetry):

    with open(foutname,'w') as fout:
        fout.write('#This file was created on {}\n'.format(date.today()))
        fout.write(
            '#The original space group number was {}\n\n'.format(spacegroup)
        )

        fout.write('data_image0\n')
        if symmetry:
            fout.write('_symmetry_space_group_name_H-M {}\n'.format(
                rspacegroups[spacegroup]))
            fout.write('_space_group_IT_number {}\n'.format(spacegroup))
        else:
            fout.write('_symmetry_space_group_name_H-M P 1\n')
            fout.write('_space_group_IT_number 1\n\n')

        fout.write('_cell_length_a {:.4f}\n'.format(cell_params[0]))
        fout.write('_cell_length_b {:.4f}\n'.format(cell_params[1]))
        fout.write('_cell_length_c {:.4f}\n'.format(cell_params[2]))
        fout.write('_cell_angle_alpha {:.4f}\n'.format(cell_params[3]))
        fout.write('_cell_angle_beta {:.4f}\n'.format(cell_params[4]))
        fout.write('_cell_angle_gamma {:.4f}\n\n'.format(cell_params[5]))

        fout.write('loop_\n')
        fout.write('_atom_site_label\n')
        fout.write('_atom_site_type_symbol\n')
        fout.write('_atom_site_fract_x\n')
        fout.write('_atom_site_fract_y\n')
        fout.write('_atom_site_fract_z\n')

        for i,a in enumerate(atoms):

            #If symmetry maintained,print out only the asymmetric unit.
            if symmetry:
                if a[0]=='T':
                    fout.write('{} {} {:.5f} {:.5f} {:.5f}\n'.format(
                        a[1]+str(i),a[1],a[2],a[3],a[4])
                    )

            #Else print all atoms.
            else:
                fout.write('{} {} {:.5f} {:.5f} {:.5f}\n'.format(
                        a[1]+str(i),a[1],a[2],a[3],a[4])
                )

#==============================================================================
#==============================================================================
#Small check on command line arguments.

if len(sys.argv)==2:
    filename=sys.argv[1]
else:
    try:
        filename=sys.argv[1]
    except IndexError:
        print('\nNo input file name was provided. Please try again.\n')

#Read input file, collect necessary information, and write output cif files.

with open(filename,'r') as fin:

    #Grab the space group number for output cif file.
    spacegroup=int(readuntil(fin,'SPACE  GROUP  N').split()[-1])

    #Check first whether the CRYSTAL optimisation actually converged.
    if readuntil(fin,'CONVERGENCE TESTS SATISFIED'):

        #Search for final set of coordinates
        readuntil(fin,'FINAL OPTIMIZED GEOMETRY')

        #Skip some blank lines
        for i in range(5):
            fin.readline()

        #Collect the cell parameters
        cell_params=fin.readline().split()
        cell_params=[float(j) for j in cell_params]

        #Skip some blank lines and collect the number of atoms
        fin.readline()
        natoms=int(fin.readline().split()[-1])
        fin.readline()
        fin.readline()

        #Collect atom types and coordinates, as well as whether they are part
        #of the asymmetric unit.
        atoms=coords(fin,natoms)

        #Generate cif file.
        foutname1=sys.argv[1]+'-symmetry-{}.cif'.format(spacegroup)
        makecif(cell_params,atoms,spacegroup,foutname1,True)

        foutname2=sys.argv[1]+'-P1.cif'
        makecif(cell_params,atoms,spacegroup,foutname2,False)

        print('\nTwo .cif files have been generated, with names ending in '+\
            '"-symmetry-{}.cif" and "-P1.cif"'.format(spacegroup))

    else:
        print('Optimisation does not appear to have successfully completed.'+\
            'Please double check output file. Not proceeding with conversion.')
        sys.exit()

