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

def acyl_group(cystype, final_res, add_acyl):
	## acyl_type bead type and charge separated by : [bead type:charge]
	## atom_type is atom name 
	## bond is the bonded information separated by : 
	## bdt is the bonded info (func, length, fc)
	## angle is the angle information separated by :
	## ant is the bonded info (func, angle, fc)
	## order is 1st number is vector scaling factor from the CYS, 2nd value is the offset on the X axis in Angstrom, 3rd value if not 0 means add to BB bead instead of SC 
	## if terminal residue uses charged bead
	if final_res == True:
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
	elif cystype == 'CYST' and final_res==True:
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
	## reads in topology seperated by '[' 
	section_number=0
	sections=[[]]
	for line in open(itp, 'r').readlines():
		if len(line.split()) > 0: 
			if line.split()[0] == '[':
				section_number+=1
				sections.append([])		
			sections[section_number].append(line)
	return sections

def update(sections, add_acyl,cystype):
	## checks if selected residue is a CYS
	is_cys=False
	for residues in sections[2]:
		if residues.split()[2] == str(add_acyl) and residues.split()[3] == 'CYS':
			is_cys = True
	if not is_cys:
		sys.exit("acyl group is not possible for this residue: "+str(cystype)+'\t'+str(add_acyl))
	##checks if selected residue is the final residue
	final_res=False
	if sections[2][-1].split()[2] == str(add_acyl):
		final_res=True

	acyl_type, atom_type, bond, angle,order = acyl_group(cystype, final_res, add_acyl) ## fetches acyl group info
	offset=0
	offset_ini=0
	bond_offset=0
	angle_offset=0
	additional_bonds, additional_angles=[],[]
	for seg_val,segment in enumerate(sections):
		cont=False
		for line_val,line in enumerate(segment):
			if seg_val < 2:  ### gumpf at start of itp file 
				new = re.sub('\n','',line)
				sections[seg_val][line_val]=new
			if seg_val == 2: ### [ ATOMS ] 
				if len(line.split()) > 5: ## removes empty lines
					if int(line.split()[2]) in [add_acyl]:  ## if correct CYS
						if line.split()[4] =='BB':  ##  updates BB bead to correct bead type
							new = '\t'+str(int(line.split()[0]))+'\t'+acyl_type[0].split(':')[0]+'\t'+str(line.split()[2])+\
							'\t'+cystype+'\t'+atom_type[0]+'\t'+str(int(line.split()[5]))+'\t'+acyl_type[0].split(':')[1]+'.0000'
							sections[seg_val][line_val]=new
						elif line.split()[4] =='SC1':
							offset_ini=int(line.split()[0])
							for atoms in range(1,len(acyl_type)): ## runs through acyl group and adds to [ ATOMS ]
								new = '\t'+str(int(line.split()[0])+atoms-1)+'\t'+acyl_type[atoms].split(':')[0]+'\t'+str(line.split()[2])\
								+'\t'+cystype+'\t'+atom_type[atoms]+'\t'+str(int(line.split()[5])+atoms-1)+'\t'+acyl_type[atoms].split(':')[1]+'.0000'
								if atoms==1:
									sections[seg_val][line_val]=new
								else:
									sections[seg_val].insert(line_val+atoms-1, new)
							offset=atoms-1
					elif int(line.split()[2]) not in [add_acyl]: ## adds other beads back into to the atoms sections with updated indexes
						new = '\t'+str(int(line.split()[0])+offset)+'\t'+line.split()[1]+'\t'+str(line.split()[2])\
						+'\t'+line.split()[3]+'\t'+line.split()[4]+'\t'+str(int(line.split()[5])+offset)+'\t'+str(line.split()[6])
						sections[seg_val][line_val]=new
			if seg_val == 3:
				if len(line.split()) >0:
					try:
						float(line.split()[0])
						if int(line.split()[1]) != offset_ini:
							col1=greater_than(line.split()[0], offset_ini, offset)
							col2=greater_than(line.split()[1], offset_ini, offset)
							if len(line.split()) == 5:	
								new = col1+'\t'+col2+'\t'+line.split()[2]+'\t'+line.split()[3]+'\t'+line.split()[4]
							if len(line.split())==7:	
								new = col1+'\t'+col2+'\t'+str(int(line.split()[2]))+'\t'+line.split()[3]+'\t'+line.split()[4]+'\t'+line.split()[5]+'\t'+line.split()[6]
							sections[seg_val][line_val]=new
						elif cont==True and int(line.split()[1]) == offset_ini:
							insertion=line_val
							for offset_lines, bond_line in enumerate(bond):
								new = str(int(bond_line.split(':')[0])+offset_ini-2)+'\t'+str(int(bond_line.split(':')[1])+offset_ini-2)\
								+bond_line.split(':')[2]+'\t;\t'+cystype
								additional_bonds.append(new)
							bond_offset=offset_lines
						if bond_offset > 0 and int(line_val)==len(sections[seg_val])-1:
								del sections[seg_val][insertion]
								for val,bond in enumerate(additional_bonds):
									sections[seg_val].insert(insertion+val,bond)
								break							
					except:
						new = re.sub('\n','',line)
						sections[seg_val][line_val]=new
						if len(line.split()) > 1:
							if line.split()[1]=='Sidechain':
								cont=True
			if seg_val == 4:
				if len(line.split()) >0:
					try:
						float(line.split()[0])
						if int(line.split()[1]) != offset_ini:
							col1=greater_than(line.split()[0], offset_ini, offset)
							col2=greater_than(line.split()[1], offset_ini, offset)
						
							new = col1+'\t'+col2+'\t'+line.split()[2]+'\t'+line.split()[3]+'\t'+line.split()[4]+'\t'+line.split()[5]
							sections[seg_val][line_val]=new
					except:
						new = re.sub('\n','',line)
						sections[seg_val][line_val]=new
			if seg_val == 5:
				if len(line.split()) >0:
					try:
						float(line.split()[0])
						if int(line.split()[1]) != offset_ini:
							col1=greater_than(line.split()[0], offset_ini, offset)
							col2=greater_than(line.split()[1], offset_ini, offset)
							col3=greater_than(line.split()[2], offset_ini, offset)

							new = col1+'\t'+col2+'\t'+col3+'\t'+line.split()[3]+'\t'+line.split()[4]+'\t'+line.split()[5]+'\t'+line.split()[6]+'\t'+line.split()[7]
							sections[seg_val][line_val]=new
						if cont==True:
							if int(line.split()[0]) > offset_ini or line_val==len(sections[seg_val])-1:
								if line_val==len(sections[seg_val])-1:
									insertion=line_val+1
								else:
									insertion=line_val
								sections[seg_val][line_val]=new
								cont=False
								for offset_angle_val, angle_line in enumerate(angle):
									new = str(int(angle_line.split(':')[0])+offset_ini-2)\
										+'\t'+str(int(angle_line.split(':')[1])+offset_ini-2)\
										+'\t'+str(int(angle_line.split(':')[2])+offset_ini-2)\
										+angle_line.split(':')[3]+'\t:\t'+cystype	
									additional_angles.append(new)
								angle_offset=offset_angle_val
						if angle_offset > 0 and line_val==len(sections[seg_val])-1:
								for val,angle in enumerate(additional_angles):
									sections[seg_val].insert(insertion+val,angle)
								break
					except:
						new = re.sub('\n','',line)
						sections[seg_val][line_val]=new
						if len(line.split()) > 1:
							if line.split()[1]=='Sidechain':
								cont=True
			if seg_val == 6:
				if len(line.split()) >0:
					try:
						float(line.split()[0])
						if int(line.split()[1]) != offset_ini:
							col1=greater_than(line.split()[0], offset_ini, offset)
							col2=greater_than(line.split()[1], offset_ini, offset)
							col3=greater_than(line.split()[2], offset_ini, offset)
							col4=greater_than(line.split()[3], offset_ini, offset)
							
							if len(line.split()) == 10:	
								new = col1+'\t'+col2+'\t'+col3+'\t'+col4+'\t'+line.split()[4]+'\t'+line.split()[5]+'\t'+line.split()[6]\
								+'\t'+line.split()[7]+'\t'+line.split()[8]+'\t'+line.split()[9]
							if len(line.split())==9:	
								new = col1+'\t'+col2+'\t'+col3+'\t'+col4+'\t'+line.split()[4]+'\t'+line.split()[5]+'\t'+line.split()[6]\
								+'\t'+line.split()[7]+'\t'+line.split()[8]					
							sections[seg_val][line_val]=new
					except:
						new = re.sub('\n','',line)
						sections[seg_val][line_val]=new
	return sections

def greater_than(value, offset_ini, offset):
	if int(value) > offset_ini:
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
	# strip information from pdb line
	try:
		# return [str(line[12:16]),str(line[17:21]),line[21],int(line[22:26]),float(line[30:38]),float(line[38:46]),float(line[46:54])]
		return dict([('atn',str(line[12:16])),('ren',str(line[17:21])),('ch',line[21]),('rid',int(line[22:26])), ('x',float(line[30:38])),('y',float(line[38:46])),('z',float(line[46:54]))])
	except:
		print('pdb line is wrong:\n'+line)
		return False

def add_acyl_tails(pdb, add_acyl,cystype,zscale):
	pdbline = "ATOM  %5d %4s %4s%1s %3d    %8.3f%8.3f%8.3f%6.2f%6.2f"
	final_res=False
	if pdbatom(pdb[-1])['rid'] == add_acyl:
		final_res=True
	is_cys=False
	for residues in pdb:
		if pdbatom(residues)['rid'] == add_acyl and pdbatom(residues)['ren'] == 'CYS ':
			is_cys = True
	if not is_cys:
		sys.exit('it is not possible to add '+str(cystype)+'\t'+str(add_acyl)+" to this residue: ")

	acyl_type, atom_type, bond, angle,order = acyl_group(cystype, final_res, add_acyl)
	CYS_new=[]
	xsc, ysc, zsc=[],[],[]
	for val, line in enumerate(pdb):
		offset=0
		try:
			int(line.split()[4])
		except:
			offset=1
		line_prev=line
		if int(line.split()[4+offset])==add_acyl:
			xsc.append(float(line.split()[5+offset]))
			ysc.append(float(line.split()[6+offset]))
			zsc.append(float(line.split()[7+offset]))
			if line.split()[2] == 'BB':
				rm=val
				line_sep_bb = pdbatom(line)
				CYS_new.append(pdbline%((add_acyl,atom_type[0],cystype,line_sep_bb['ch'],line_sep_bb['rid'],line_sep_bb['x'],line_sep_bb['y'],line_sep_bb['z'],1,0)))
				vxbb=float(line.split()[5+offset])-xbb
				vybb=float(line.split()[6+offset])-ybb
				vzbb=float(line.split()[7+offset])-zbb
		if int(line.split()[4+offset])==add_acyl and line.split()[2] == 'SC1':
			vxsc=xsc[1]-xsc[0]
			vysc=ysc[1]-ysc[0]
			vzsc=zsc[1]-zsc[0]
			line_sep = pdbatom(line)
			CYS_new.append(pdbline%((add_acyl,atom_type[1],cystype,line_sep['ch'],line_sep['rid'],line_sep['x'],line_sep['y'],line_sep['z'],1,0)))
			for atom in range(2, len(atom_type)):
				if int(order[atom].split(':')[2]) == 0:
					CYS_new.append(pdbline%((add_acyl,atom_type[atom],cystype,line_sep['ch'],line_sep['rid'],line_sep['x']\
					+(vxsc*int(order[atom].split(':')[0]))+int(order[atom].split(':')[1]),line_sep['y']+(vysc*int(order[atom].split(':')[0])),\
					line_sep['z']+float(vzsc*int(order[atom].split(':')[0]))*zscale,1,0)))
					print(float(vzsc*int(order[atom].split(':')[0]))*zscale)
				else:
					CYS_new.append(pdbline%((add_acyl,atom_type[atom],cystype,line_sep_bb['ch'],line_sep_bb['rid'],line_sep_bb['x']\
					+(vxsc*int(order[atom].split(':')[0]))+vxbb,line_sep_bb['y']+(vysc*int(order[atom].split(':')[0]))+vybb,\
					line_sep_bb['z']+float((vzsc*int(order[atom].split(':')[0]))+vzbb)*zscale,1,0)))
					print(float(vzsc*int(order[atom].split(':')[0]))*zscale)
			for i in CYS_new:
				print(i)
			break
		xbb=float(line.split()[5+offset])
		ybb=float(line.split()[6+offset])
		zbb=float(line.split()[7+offset])	

	for val,atom in enumerate(CYS_new):
		pdb.insert(rm+val,atom)	
	del pdb[rm+val+1:rm+val+3]
	return pdb

parser = argparse.ArgumentParser()

parser.add_argument('-itp', help='cg itp file',metavar='protein.itp', type=str)
parser.add_argument('-cys', help='cysteines to acylate ',metavar='22',type=int, nargs='*')
parser.add_argument('-pdb', help='pdb file',metavar='protein.pdb', type=str)
parser.add_argument('-res', help='acyl type',metavar='CYSM, CYSD or CYST', type=str, nargs='*')
parser.add_argument('-z', help='z scaling factor ',metavar='between -1 and 1 scales z orientations',type=float, nargs='*')


args = parser.parse_args()
options = vars(args)
timestamp =  strftime("%Y-%m-%d_%H-%M-%S", gmtime())

if args.itp is not None:
	sections = read_topology(args.itp)
	for val, cysteines in enumerate(args.cys):
		sections = update(sections, cysteines, args.res[val])
	with open('acyl_'+args.itp.split('/')[-1], 'w') as file_write:
		for seg_val,segment in enumerate(sections):
			for line_val,line in enumerate(segment):
				file_write.write(line+'\n')
			file_write.write('\n')
if args.pdb is not None:
	pdb=read_pdb(args.pdb)
	for val, cysteines in enumerate(args.cys):
		pdb = add_acyl_tails(pdb, cysteines, args.res[val], args.z[val])
	with open('acyl_'+args.pdb.split('/')[-1], 'w') as file_write:
		for line_val,line in enumerate(pdb):
				file_write.write(line+'\n')
