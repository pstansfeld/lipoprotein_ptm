import os, sys 
import numpy as np
import argparse
from subprocess import Popen, PIPE
import subprocess, shlex
from time import gmtime, strftime
from shutil import copyfile
from pathlib import Path
from collections import Counter  
import re
import time

parser = argparse.ArgumentParser()

parser.add_argument('-itp', help='cg itp file',metavar='protein.itp', type=str, required=True)
parser.add_argument('-cys', help='cysteines to acylate ',metavar='22',type=int, nargs='*')
parser.add_argument('-c', help='pdb file',metavar='protein.pdb', type=str, required=True)
parser.add_argument('-res', help='acyl type',metavar='CYSM, CYSD or CYST', type=str, nargs='*', required=True)
parser.add_argument('-z', help='z scaling factor ',metavar='between -1 and 1 scales z orientations',type=float, nargs='*')


args = parser.parse_args()
options = vars(args)



def acyl_group(cystype, initial_res, add_acyl):
	## acyl_type bead type and charge separated by : [bead type:charge]
	## atom_type is atom name 
	## bond is the bonded information separated by : 
	## bdt is the bonded info (func, length, fc)
	## angle is the angle information separated by :
	## ant is the bonded info (func, angle, fc)
	## order is 1st number is vector scaling factor from the CYS, 2nd value is the offset on the X axis in Angstrom, 3rd value if not 0 means add to BB bead instead of SC 
	## if terminal residue uses charged bead
	if initial_res == True:
		acyl_type = ['Qd:1']
	else:
		acyl_type = ['P5:0']

	if cystype == 'CYSD':
		acyl_type+=['C5:0',  'Na:0', 'Na:0', 'C1:0', 'C3:0', 'C1:0', 'C1:0', 'C1:0', 'C1:0', 'C1:0','C1:0']
		atom_type = ['BB', 'SC1','GL1','GL2','C1A','D2A','C3A','C4A','C1B','C2B','C3B','C4B']
		bdt=[':   1	0.31    7500',': 	1 	0.37 	1250',': 	1 	0.47 	1250']
		bond = ['1:2'+bdt[0], '2:3'+bdt[1], '3:4'+bdt[1], '3:5'+bdt[2],'5:6'+bdt[2],'6:7'+bdt[2],'7:8'+bdt[2],'4:9'+bdt[2],\
		'9:10'+bdt[2],'10:11'+bdt[2],'11:12'+bdt[2]]	
		ant=[': 	2 	180.0 	25.0', ': 	2 	120.0 	25.0', ': 	2 	180.0 	25.0', ': 	2 	120.0 	45.0']
		angle = ['1:2:3'+ant[0],'2:3:4'+ant[1],'2:3:5'+ant[2],'3:5:6'+ant[2],'5:6:7'+ant[3],\
		'6:7:8'+ant[2],'4:9:10'+ant[2],'9:10:11'+ant[2],'10:11:12'+ant[2]]
		order     = ['0:0:0','0:0:0','1:2:0','1:0:0','2:0:0','3:0:0','4:0:0','5:0:0','2:2:0','3:2:0','4:2:0','5:2:0']
	elif cystype == 'CYSP':
		acyl_type+=['Na:0','C1:0','C1:0','C1:0','C1:0']
		atom_type = ['BB','SC1','C1A','C2A','C3A','C4A']
		bdt=[':   1	0.31    7500',': 	1 	0.37 	1250',': 	1 	0.47 	1250']
		bond = ['1:2'+bdt[0],'2:3'+bdt[1],'3:4'+bdt[2],'4:5'+bdt[2],'5:6'+bdt[2]]
		ant=[': 	2 	180.0 	25.0']
		angle = ['1:2:3'+ant[0],'2:3:4'+ant[0],'3:4:5'+ant[0],'4:5:6'+ant[0]]
		order = ['0:0:0','0:0:0','1:0:0','2:0:0','3:0:0','4:0:0']
	elif cystype == 'CYST' and initial_res==True:
		acyl_type = ['P5:0', 'C5:0', 'Na:0', 'Na:0', 'C1:0', 'C3:0', 'C1:0', 'C1:0', 'C1:0', 'C1:0', 'C1:0', 'C1:0', 'C1:0', 'C1:0', 'C1:0', 'C1:0']
		atom_type = ['BB','SC1','GL1','GL2','C1A','D2A','C3A','C4A','C1B','C2B','C3B','C4B','C1C','C2C','C3C','C4C']
		order     = ['0:0:0','0:0:0','1:2:0','1:0:0','2:0:0','3:0:0','4:0:0','5:0:0','2:2:0','3:2:0','4:2:0','5:2:0','1:0:2','2:0:2','3:0:2','4:0:2']
		bdt=[':   1	0.31    7500',': 	1 	0.37 	1250',': 	1 	0.47 	1250']
		bond = ['1:2'+bdt[0], '2:3'+bdt[1], '3:4'+bdt[1], '3:5'+bdt[2],'5:6'+bdt[2],'6:7'+bdt[2],'7:8'+bdt[2],'4:9'+bdt[2],\
		'9:10'+bdt[2],'10:11'+bdt[2],'11:12'+bdt[2],'1:13'+bdt[1],'13:14'+bdt[2],'14:15'+bdt[2],'15:16'+bdt[2]]
		ant=[': 	2 	180.0 	25.0', ': 	2 	120.0 	25.0', ': 	2 	180.0 	25.0', ': 	2 	120.0 	45.0']
		angle = ['1:2:3'+ant[0],'2:3:4'+ant[1],'2:3:5'+ant[2],'3:5:6'+ant[2],'5:6:7'+ant[3],'6:7:8'+ant[2],'4:9:10'+ant[2],\
		'9:10:11'+ant[2],'10:11:12'+ant[2],'1:13:14'+ant[2],'13:14:15'+ant[2],'14:15:16'+ant[2]]

	else:
		sys.exit("acyl group is not possible for this residue: "+str(cystype)+'\t'+str(add_acyl))
	return acyl_type, atom_type, bond, angle, order


def read_topology(itp):
	## reads in topology seperated by '[', sections = [[moleculetype],[atoms],...] 
	section_number=0
	sections={}
	block='header'
	sections[block]=[]
	for line in open(itp, 'r').readlines():
		if len(line.split()) > 0: 
			if line.split()[0] == '[':
				block = line.split()[1]
				sections[block]=[]
			sections[block].append(line.rstrip())
	return sections

def check_cysteine(section, add_acyl,cystype):
	## checks if selected residue is a CYS
	is_cys=False
	initial_res=True
	for atoms in section:
		try:
			int(atoms.split()[0])
			if atoms.split()[2] == str(add_acyl) and atoms.split()[3] == 'CYS':
				is_cys = True
				return initial_res
			initial_res=False
		except:
			pass
	if not is_cys:
		sys.exit("acyl group is not possible for this residue: "+str(cystype)+'\t'+str(add_acyl))

def update_moleculetype(moleculetype_section):
	for moleculetype_val, moleculetype_line in enumerate(moleculetype_section):
		if moleculetype_line.split()[0] not in ['[', ';']:
			moleculetype_section[moleculetype_val]='acyl_'+moleculetype_line
	return moleculetype_section

def update_atoms(atoms_section, add_acyl,cystype):
	offset=0
	for at_val, atom in enumerate(atoms_section):
		try:
			int(atom.split()[0]) ## removes empty lines
			if int(atom.split()[2]) in [add_acyl]:  ## if correct CYS
				if atom.split()[4] =='BB':  ##  updates BB bead to correct bead type
					new = '    '+str(int(atom.split()[0]))+'\t'+acyl_type[0].split(':')[0]+'\t'+str(atom.split()[2])+\
					'\t'+cystype+'\t'+atom_type[0]+'\t'+str(int(atom.split()[5]))+'\t'+acyl_type[0].split(':')[1]+'.0000'
					atoms_section[at_val]=new
				elif atom.split()[4] =='SC1':
					offset_initial=int(atom.split()[0])
					for atoms in range(1,len(acyl_type)): ## runs through acyl group and adds to [ ATOMS ]
						new = '    '+str(int(atom.split()[0])+atoms-1)+'\t'+acyl_type[atoms].split(':')[0]+'\t'+str(atom.split()[2])\
						+'\t'+cystype+'\t'+atom_type[atoms]+'\t'+str(int(atom.split()[5])+atoms-1)+'\t'+acyl_type[atoms].split(':')[1]+'.0000'
						if atoms==1:
							atoms_section[at_val]=new
						else:
							atoms_section.insert(at_val+atoms-1, new)
					offset=atoms-1
			elif int(atom.split()[2]) not in [add_acyl]: ## adds other beads back into to the atoms sections with updated indexes
				atoms_section[at_val] = '    '+str(int(atom.split()[0])+offset)+'\t'+atom.split()[1]+'\t'+str(atom.split()[2])\
				+'\t'+atom.split()[3]+'\t'+atom.split()[4]+'\t'+str(int(atom.split()[5])+offset)+'\t'+str(atom.split()[6])
		except:
			pass
	return atoms_section, offset_initial, offset

def update_bonds(bonds_section, add_acyl,cystype):
	additional_bonds=[]
	sidechain=False
	for bond_val, bond in enumerate(bonds_section):
		if bond == '; Sidechain bonds':
			sidechain=True
		try:
			int(bond.split()[0]) ## removes empty lines
			if int(bond.split()[1]) != offset_initial:
				col1=greater_than(bond.split()[0], offset_initial, offset)
				col2=greater_than(bond.split()[1], offset_initial, offset)
				try:
					int(bond.split()[4])
					new = col1+'\t'+col2+'\t'+str(int(bond.split()[2]))+'\t'+bond.split()[3]+'\t'+bond.split()[4]+'\t'+bond.split()[5]+'\t'+bond.split()[6]
					
				except:	
					new = col1+'\t'+col2+'\t'+bond.split()[2]+'\t'+bond.split()[3]+'\t'+bond.split()[4]
				bonds_section[bond_val]=new
			elif sidechain:
				insertion=bond_val	
		except:
			pass
	del bonds_section[insertion]
	for offset_lines, bond_cys_line in enumerate(bond_cys):
		new = str(int(bond_cys_line.split(':')[0])+offset_initial-2)+'\t'+str(int(bond_cys_line.split(':')[1])+offset_initial-2)+bond_cys_line.split(':')[2]+'\t;\t'+cystype
		bonds_section.insert(insertion+offset_lines,new)
	return bonds_section

def update_contraints(contraints_section, add_acyl,cystype):
	for contraint_val, contraint in enumerate(contraints_section):
		try:
			int(contraint.split()[0])
			if int(contraint.split()[1]) != offset_initial:
				col1=greater_than(contraint.split()[0], offset_initial, offset)
				col2=greater_than(contraint.split()[1], offset_initial, offset)
			
				new = col1+'\t'+col2+'\t'+contraint.split()[2]+'\t'+contraint.split()[3]+'\t'+contraint.split()[4]+'\t'+contraint.split()[5]
				contraints_section[contraint_val]=new
		except:
			pass
	return contraints_section

def update_angle(angles_section, add_acyl,cystype):
	for angle_val, angle in enumerate(angles_section):
		try:
			int(angle.split()[0])
			if int(angle.split()[1]) != offset_initial:
				col1=greater_than(angle.split()[0], offset_initial, offset)
				col2=greater_than(angle.split()[1], offset_initial, offset)
				col3=greater_than(angle.split()[2], offset_initial, offset)
				angles_section[angle_val] = col1+'\t'+col2+'\t'+col3+'\t'+angle.split()[3]+'\t'+angle.split()[4]+'\t'+angle.split()[5]+'\t'+angle.split()[6]+'\t'+angle.split()[7]
		except:
			pass
	angles_section.append('; '+cystype)
	for offset_angle_val, angle_line in enumerate(angle_cys):
		new = str(int(angle_line.split(':')[0])+offset_initial-2)+'\t'+str(int(angle_line.split(':')[1])+offset_initial-2)+'\t'+str(int(angle_line.split(':')[2])+offset_initial-2)\
		+angle_line.split(':')[3]+'\t:\t'+cystype
		angles_section.append(new)

	return angles_section

def update_dihedrals(dihedrals_section, add_acyl,cystype):
	for dihedral_val, dihedral in enumerate(dihedrals_section):
		try:
			int(dihedral.split()[0])
			if int(dihedral.split()[1]) != offset_initial:
				col1=greater_than(dihedral.split()[0], offset_initial, offset)
				col2=greater_than(dihedral.split()[1], offset_initial, offset)
				col3=greater_than(dihedral.split()[2], offset_initial, offset)
				col4=greater_than(dihedral.split()[3], offset_initial, offset)
				try:
					int(dihedral.split()[7])
					new = col1+'\t'+col2+'\t'+col3+'\t'+col4+'\t'+line.split()[4]+'\t'+line.split()[5]+'\t'+line.split()[6]\
					+'\t'+line.split()[7]+'\t'+line.split()[8]+'\t'+line.split()[9]
				except:
					new = col1+'\t'+col2+'\t'+col3+'\t'+col4+'\t'+line.split()[4]+'\t'+line.split()[5]+'\t'+line.split()[6]\
					+'\t'+line.split()[7]+'\t'+line.split()[8]
				dihedrals_section[dihedral_val]=new
		except:
			pass
	return dihedrals_section

def greater_than(value, offset_initial, offset):
	if int(value) > offset_initial:
		new_number=str(int(value)+offset)
	else:
		new_number=str(value)
	return new_number

def read_pdb(pdb_file):
	pdb=[]
	box_line=''
	for line in open(pdb_file, 'r').readlines():
		line = re.sub('\n','',line)
		if not line[0] in ['#', '@']:
			if line.split()[0] == 'ATOM':
				pdb.append(line)
	return pdb

def pdbatom(line):
### get information from pdb file
### atom number, atom name, residue name,chain, resid,  x, y, z, backbone (for fragment), connect(for fragment)
	try:
		return dict([('atom_number',str(line[7:11]).replace(" ", "")),('atom_name',str(line[12:16]).replace(" ", "")),('residue_name',str(line[17:21]).replace(" ", "")),('chain',line[21]),('residue_id',int(line[22:26])), ('x',float(line[30:38])),('y',float(line[38:46])),('z',float(line[46:54]))])
	except:
		sys.exit('\npdb line is wrong:\t'+line) 

def add_acyl_tails(pdb, add_acyl,cystype,zscale):
	pdbline = "ATOM  %5d %4s %4s%1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f"

	pdb_new={}
	count=0

	for val, line in enumerate(pdb):
		line_sep = pdbatom(line)
		if line_sep['residue_id']==add_acyl:
			if line_sep['atom_name'] == 'BB':
				bb=[line_sep['x'], line_sep['y'], line_sep['z']]
				pdb_new[val]={'atom_name':atom_type[0],'residue_id':line_sep['residue_id'], 'residue_name':cystype, 'x':line_sep['x'],'y':line_sep['y'],'z':line_sep['z']}

			elif line_sep['atom_name'] == 'SC1':
				pdb_new[val]={'atom_name':atom_type[1],'residue_id':line_sep['residue_id'], 'residue_name':cystype, 'x':line_sep['x'],'y':line_sep['y'],'z':line_sep['z']}
				sc=[line_sep['x'], line_sep['y'], line_sep['z']]
				vxsc=sc[0]-bb[0]
				vysc=sc[1]-bb[1]
				vzsc=sc[2]-bb[2]
				for atom in range(2, len(atom_type)):
					count+=1
					if int(order[atom].split(':')[2]) == 0:
						x=sc[0]+(vxsc*int(order[atom].split(':')[0]))+int(order[atom].split(':')[1])
						y=sc[1]+(vysc*int(order[atom].split(':')[0]))+int(order[atom].split(':')[1])
						z=sc[2]+(vzsc*int(order[atom].split(':')[0]))*zscale
						pdb_new[val+count]={'atom_name':atom_type[atom],'residue_id':line_sep['residue_id'], 'residue_name':cystype, 'x':x,'y':y,'z':z}
					else:
						x=bb[0]+(vxsc*int(order[atom].split(':')[0]))+int(order[atom].split(':')[1])
						y=bb[1]+(vysc*int(order[atom].split(':')[0]))+int(order[atom].split(':')[1])
						z=bb[2]+(vzsc*int(order[atom].split(':')[0]))*zscale
						pdb_new[val+count]={'atom_name':atom_type[atom],'residue_id':line_sep['residue_id'], 'residue_name':cystype, 'x':x,'y':y,'z':z}

					
		else:
			pdb_new[val+count]={'atom_name':line_sep['atom_name'],'residue_id':line_sep['residue_id'], 'residue_name':line_sep['residue_name'], 'x':line_sep['x'],'y':line_sep['y'],'z':line_sep['z']}
	pdb=[]
	for line_index in pdb_new:
		pdb.append(pdbline%((int(line_index), pdb_new[line_index]['atom_name'],pdb_new[line_index]['residue_name'],' ', pdb_new[line_index]['residue_id'],pdb_new[line_index]['x'],pdb_new[line_index]['y'],pdb_new[line_index]['z'],1,0)))
	return pdb

def collect_input():
#### collates all input files in input directory
	copyfile(args.c, working_dir+args.c.split('/')[-1])
	copyfile(args.itp, working_dir+args.itp.split('/')[-1])
	os.chdir(working_dir)
#### converts input files into pdb files 
	gromacs('gmx editconf -f '+args.c.split('/')[-1]+' -o CG_input.pdb')

def mkdir_directory(directory):
#### checks if folder exists, if not makes folder
	if not os.path.exists(directory):
		os.mkdir(directory)

### runs gromacs commands
def gromacs(cmd):
#### if the flag gromacs is used every gromacs command will be printed to the terminal 
	output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	err, out = output.communicate()
	exitcode = output.returncode
	out=out.decode("utf-8")
#### all gromacs outputs will be saved into gromacs_outputs within the folder it is run
	with open('gromacs_outputs', 'a') as checks:
		checks.write(out)
#### standard catch for failed gromacs commands
		if 'File input/output error:' in out:
			sys.exit('\n'+out)
		elif 'Error in user input:' in out:
			sys.exit('\n'+out)
		elif 'did not converge to Fmax ' in out:
			sys.exit('\n'+out)
		elif 'Segmentation fault (core dumped):' in out:
			sys.exit('\n'+out)
		elif 'Fatal error:' in out:
			sys.exit('\n'+out)
		elif 'but did not reach the requested Fmax' in out:
			sys.exit('\n'+out)

gmx='gmx'
pdbline = "ATOM  %5d %4s %4s%1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f"

### sets up file locations

timestamp =  strftime("%Y-%m-%d_%H-%M-%S", gmtime())

start_dir=os.getcwd()+'/'  ### initial working directory
working_dir=os.getcwd()+'/add_acyl_'+timestamp+'/'   ### working directory 
script_dir=os.path.dirname(os.path.realpath(__file__))+'/'



mkdir_directory(working_dir)
script_dir=os.path.dirname(os.path.realpath(__file__))+'/'

user_at_input = collect_input()


itp_file=args.itp.split('/')[-1]

sections = read_topology(itp_file)
pdb=read_pdb('CG_input.pdb')

for val, cysteines in enumerate(args.cys):
	initial_res = check_cysteine(sections['atoms'], cysteines, args.res[val])
	acyl_type, atom_type, bond_cys, angle_cys,order = acyl_group(args.res[val], initial_res, cysteines)
	sections['moleculetype']=update_moleculetype(sections['moleculetype'])
	sections['atoms'], offset_initial, offset = update_atoms(sections['atoms'], cysteines, args.res[val])
	sections['bonds'] = update_bonds(sections['bonds'], cysteines, args.res[val])
	sections['constraints'] = update_contraints(sections['constraints'], cysteines, args.res[val])
	sections['angles'] = update_angle(sections['angles'], cysteines, args.res[val])
	sections['dihedrals'] = update_dihedrals(sections['dihedrals'], cysteines, args.res[val])
	pdb = add_acyl_tails(pdb, cysteines, args.res[val], args.z[val])

with open('acyl_'+itp_file, 'w') as file_write:
	for segment in sections:
		for line in sections[segment]:
			file_write.write(line+'\n')
		file_write.write('\n')
	
with open('acyl_CG_input.pdb', 'w') as file_write:
	for line_val,line in enumerate(pdb):
			file_write.write(line+'\n')
