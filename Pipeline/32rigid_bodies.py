import sys
# from bfuncs import *
from math import *
import copy
import re
import numpy

from Scientific.Geometry import Vector, isVector, Tensor, isTensor
from Scientific.indexing import index_expression
from Scientific import N

from MMTK import *
from MMTK.MoleculeFactory import *
from MMTK.ForceFields import Amber99ForceField
from MMTK import Universe
from MMTK.InternalCoordinates import *
from MMTK.ChemicalObjects import *


########################################
########## AUX FUNCTIONS ###############
########################################

# START FUNC bsubstr @@@@@
def bsubstr(astring, pos, nchar):
	"""Awk substr equivalent""" 
	result = 0;
	l = len(astring);
	if((pos > 0) & (nchar <= l)):
		f = pos+nchar-1;
		if(f <= l):
			return astring[(pos-1):f];
	else:
		return result;

# END FUNC   bsubstr @@@@@


# START FUNC bsubstr @@@@@
def bcontains(basket, apple):
	"""
	Checks if apple is in the basket
	@basket: container
	@basket type: list
	@apple: subset
	@apple type type: list
	@returns 1 if apple is contained in basket
	"""
	found = 0
	lapple = len(apple)
	lbasket = len(basket)

	if lapple >= lbasket : return 0

	for i in range(lapple) :
		found = 0
		for j in range(lbasket) :
			if apple[i] == basket[j]:
				found = 1
				break
		if found == 0 :
			return 0

	return 1

# START FUNC bdiff @@@@@
def bPDBprint_atom(record, name, position):
	serial = 1
	altLoc = ' '
	resName = 'LIG'
	chainID = 'X'
	resSeq = 1
	occupancy = 1.0
	tempFactor = 0.0
	print '{0:<6}{1:>5}{2:>4}{3} {4} {5}{6:>4}    {7:8.3f}{8:8.3f}{9:8.3f}  {10:4.2f}{11:6.2f}'.format\
		(record, serial, name, altLoc, resName, chainID, resSeq, \
		position.x(), position.y(), position.z(), occupancy, tempFactor)
# END FUNC   bdiff  @@@@@


# START FUNC bisSaturated @@@@@
def bisSaturated(element, bonds):
	"""
	Checks if an atom is involved in multiple bonds
	@element: the atom
	@element type: string
	@bonds: how many bonds is the atom involved in
	@bonds type: integer
	@rvalue: difference between the valence and bonds (unsaturated value)
	@rtype: integer
	"""
	if element == 'C' :
		return 4 - bonds
	elif element == 'Si' :
		return 4 - bonds
	elif element == 'N' :
		return 3 - bonds
	elif element == 'O' :
		return 2 - bonds
	elif element == 'S' :
		return 2 - bonds
	elif element == 'P' :
		return 5 - bonds
	elif element == 'Se' :
		return 2 - bonds
	else :
		return 0

####################################
########## FUNCTIONS ###############
####################################

def anchor(universe):
	"""Anchor a universe"""
	if(universe.numberOfAtoms() >=3):
		#set anchor atoms indexes	
		a1i	= 0
		for i in range(len(universe.atomList()[a1i].bondedTo())):
			for j in range(len(universe.atomList()[a1i].bondedTo()[i].bondedTo())):
				if universe.atomList()[a1i] is not \
					universe.atomList()[a1i].bondedTo()[i].bondedTo()[j] :
					a2i = i;
					a3i = j;
					break
		print 'anchor using: ', a1i, a2i, a3i
		
		v1 = universe.atomList()[a1i].position();
		v2 = universe.atomList()[a1i].bondedTo()[a2i].position();
		v3 = universe.atomList()[a1i].bondedTo()[a2i].bondedTo()[a3i].position();
	
		vT01 = -v1;
		tra01 = Translation(vT01);
		universe.applyTransformation(tra01);
	
	
		x = universe.atomList()[a1i].bondedTo()[a2i].position().x();
		y = 0;
		z = universe.atomList()[a1i].bondedTo()[a2i].position().z();
		vR01 = Vector(x,y,z);
		vR02 = Vector(1,0,0);
	
		sign = 1.0
		if ((z > 0) and (x > 0)) or ((z < 0) and (x < 0)):
			sign = -1.0
		
		rot01 = Rotation(Vector(0,1,0), sign*vR01.angle(vR02));
		universe.applyTransformation(rot01);
	
		v4 = universe.atomList()[a1i].bondedTo()[a2i].position();
	
		sign = -1.0
		if ((x > 0) and (y > 0)) or ((x < 0) and (y < 0)):
			sign = 1.0
		rot02 = Rotation(Vector(0,0,1), sign*v4.angle(Vector(1,0,0)));
		universe.applyTransformation(rot02);
	
		x = 0;
		y = universe.atomList()[a1i].bondedTo()[a2i].bondedTo()[a3i].position().y();
		z = universe.atomList()[a1i].bondedTo()[a2i].bondedTo()[a3i].position().z();
	
		sign = 1.0
		if ((z > 0) and (y > 0)) or ((z < 0) and (y < 0)):
			sign = -1.0
	
		vR03 = Vector(x,y,z);
		vR04 = Vector(0,1,0);
		rot03 = Rotation(Vector(1,0,0), sign*vR03.angle(vR04));
		universe.applyTransformation(rot03);
	
	else:
		print "Molecule too little to be anchored";
#--------------------------------------------------



####################################
############ CLASSES ###############
####################################

class bConvertor():
	"""
	Converts from cartesian coordinates
	found in universe to bonds, angles, dihedrals
	lists found within self
	"""
	def __init__(self, universe):
		"""
		:param universe: an MMTK universe
		:type universe: MMTK universe
		"""
		#Variables
		self.bonds = []		#list of bonds
		self.angles = []	#list of angles
		self.dihedrals = [] #list of dihedrals

	def cart2all_internals(self):
		k = -1
		# Bonds
		for i in universe[0].bonds :
			k += 1;
			x = universe.distance(i.a1, i.a2);
			self.bonds.append([i.a1, i.a2, x]);
		
		# Angles
		k = -1;
		for i in self.bonds :
			for j in self.bonds:
				if j != i :
					if (i[0] == j[0]):
						k += 1;
						x = universe.angle(i[1], i[0], j[1]);
						self.angles.append([[i[1], i[0], i[2]],\
											[i[0], j[1], j[2]],\
										 	x]);
					elif (i[0] == j[1]):
						k += 1;
						x = universe.angle(i[1], i[0], j[0]);
						self.angles.append([[i[1], i[0], i[2]],\
											[i[0], j[0], j[2]],\
										 	x]);
					elif (i[1] == j[0]):
						k += 1;
						x = universe.angle(i[0], i[1], j[1]);
						self.angles.append([[i[0], i[1], i[2]],\
											[i[1], j[1], j[2]],\
										 	x]);
					elif (i[1] == j[1]):
						k += 1;
						x = universe.angle(i[0], i[1], j[0]);
						self.angles.append([[i[0], i[1], i[2]],\
											[i[1], j[0], j[2]],\
										 	x]);
		
		# Dihedrals
		k = -1;
		for i in self.angles :
			for j in self.angles:
				if j < i :
					if (i[0] == j[1]) :
						k += 1;
						x = universe.dihedral(j[0][0], j[0][1], j[1][1], i[1][1]);
						self.dihedrals.append([j,i,x])
					if (i[1] == j[0]) :
						k += 1;
						x = universe.dihedral(i[0][0], i[0][1], i[1][1], j[1][1]);
						self.dihedrals.append([i,j,x])
	#---------------------------------------------
	
	def print_internals(self):
		"""Writes Bonds, Angles & Dihedrals to stdout"""
		print "Bonds";
		for i in self.bonds :
			print i;
		print "Angles";
		for i in self.angles :
			print i
		print "Dihedrals"
		for i in self.dihedrals :
			print i;
	#---------------------------------------------
	
	def swap_bond(self, bond):
		"""Reverses the atoms order"""
		"""in a [a1,a2,val] type bond"""
		rbond = copy.deepcopy(bond)
		x = rbond[0]
		rbond[0] = rbond[1]
		rbond[1] = x
		return rbond
	#---------------------------------------------
	
	def swap_angle(self, angle):
		"""Reverses the atoms order in a"""
		"""[[a1,a2,bval], [a2,a3,bval], aval] type angle"""
		rangle = copy.deepcopy(angle)
		a1 =  angle[0][0]
		a2l = angle[0][1]
		bvall = angle[0][2]
		a2r = angle[1][0]
		a3  = angle[1][1]
		bvalr = angle[1][2]
		
		rangle[0][0] = a3
		rangle[0][1] = a2r
		rangle[0][2] = bvalr
		rangle[1][0] = a2l
		rangle[1][1] = a1
		rangle[1][2] = bvall
		return rangle
	#---------------------------------------------
	
	def swap_dihe(self, dihe):
		"""Reverses the atoms order in a"""
		"""[ [[a1,a2,bv1],[a2l,a3l,bv2l], av1],"""
		"""  [[a2r,a3r,bv2r],[a3,a4,bv3], av2],"""
		"""  dv ] type dihedral angle"""
		rdihe = copy.deepcopy(dihe)
	
		a1   = dihe[0][0][0]
		a2   = dihe[0][0][1]
		bv1  = dihe[0][0][2]
		a2l  = dihe[0][1][0]
		a3l  = dihe[0][1][1]
		bv2l = dihe[0][1][2]
		av1  = dihe[0][2]
		a2r  = dihe[1][0][0]
		a3r  = dihe[1][0][1]
		bv2r = dihe[1][0][2]
		a3   = dihe[1][1][0]
		a4   = dihe[1][1][1]
		bv3  = dihe[1][1][2]
		av2  = dihe[1][2]
	
		rdihe[0][0][0] = a4
		rdihe[0][0][1] = a3
		rdihe[0][0][2] = bv3
		rdihe[0][1][0] = a3r
		rdihe[0][1][1] = a2r
		rdihe[0][1][2] = bv2r
		rdihe[0][2] = av2
		rdihe[1][0][0] = a3l
		rdihe[1][0][1] = a2l
		rdihe[1][0][2] = bv2l
		rdihe[1][1][0] = a2
		rdihe[1][1][1] = a1
		rdihe[1][1][2] = bv1
		rdihe[1][2] = av1
		return rdihe
	#---------------------------------------------

	def sprint_zmat_dihe(self, dihe):
		"""Write dihedral to stdout"""
		""" in Z-matrix format"""
		a1   = dihe[0][0][0]
		a2   = dihe[0][0][1]
		bv1  = dihe[0][0][2]
		a2l  = dihe[0][1][0]
		a3l  = dihe[0][1][1]
		bv2l = dihe[0][1][2]
		av1  = dihe[0][2]
		a2r  = dihe[1][0][0]
		a3r  = dihe[1][0][1]
		bv2r = dihe[1][0][2]
		a3   = dihe[1][1][0]
		a4   = dihe[1][1][1]
		bv3  = dihe[1][1][2]
		av2  = dihe[1][2]
		return '{0:>4} {1:>4} {2:>4} {3:>4} {4:>10.6f} {5:>4} {6:>10.6f} {7:>4} {8:>10.6f}'.format(a1.index+1, a1.name, a1.name[0], a2.index+1, bv1*10.0, a3.index+1, av1, a4.index+1, dihe[2])
	#---------------------------------------------
	
	
	def print_zmat_dihe(self, dihe):
		"""Write dihedral to stdout"""
		""" in Z-matrix format"""
		a1   = dihe[0][0][0]
		a2   = dihe[0][0][1]
		bv1  = dihe[0][0][2]
		a2l  = dihe[0][1][0]
		a3l  = dihe[0][1][1]
		bv2l = dihe[0][1][2]
		av1  = dihe[0][2]
		a2r  = dihe[1][0][0]
		a3r  = dihe[1][0][1]
		bv2r = dihe[1][0][2]
		a3   = dihe[1][1][0]
		a4   = dihe[1][1][1]
		bv3  = dihe[1][1][2]
		av2  = dihe[1][2]
		print '{0:>4} {1:>4} {2:>4} {3:>4} {4:>10.6f} {5:>4} {6:>10.6f} {7:>4} {8:>10.6f}'.format(a1.index+1, a1.name, a1.name[0], a2.index+1, bv1*10.0, a3.index+1, av1, a4.index+1, dihe[2])
	#---------------------------------------------


	def NeRFinte2cart(self, dihe):	# Needs debug
		"""
		Reconstructs atomlist positions
		2005 Parsons et al - NeRF algorithm
		:param universe: MMTK universe
		:type universe: MMTK universe
		:param dihe: list of dihedrals
		:type dihe: [[[],[],aval] , [[],[],aval] , dval] where [] is a bond
		"""
		known_atoms = []
		unknown_atoms = copy.deepcopy(universe.atomList())
		C = Vector(0,0,0)
		D = Vector(0,0,0)
		D2 = Vector(0,0,0)
		l = len(dihe)
		for di in range(l):	# Create the START stack
			dihe.append(self.swap_dihe(dihe[di]))
	
		#for i in dihe:
		#	print_zmat_dihe(i)
	
		# Assign initial values
		i = dihe[0]
		a1 = i[0][0][0]
		a2 = i[0][0][1]
		a3 = i[1][1][0]
		a4 = i[1][1][1]
		R1 = i[0][0][2]
		R2 = i[1][1][2]
		theta1 = i[0][2]
		phi = i[2]
		universe.atomList()[a1.index].setPosition(Vector(0,0,0))
		universe.atomList()[a2.index].setPosition(Vector(R1,0,0))
		universe.atomList()[a3.index].setPosition(Vector(R1-R2*cos(theta1), R2*sin(theta1), 0))
		known_atoms.append(a1)
		known_atoms.append(a2)
		known_atoms.append(a3)
		#print 'initialize at: ', a1, a2, a3
		#bPDBprint_atom('HETATM', a1.name, 10*a1.position())
		#bPDBprint_atom('HETATM', a2.name, 10*a2.position())
		#bPDBprint_atom('HETATM', a3.name, 10*a3.position())
		# Main loop
		cnt=0;
		while len(known_atoms) < len(universe.atomList()):
			cnt += 1
			if cnt > 9999:
				print 'NeRFinte2cart() cnt EXCEEDED: 9999'
				break
			ai = -1
			for a in unknown_atoms:
				ai += 1
				for i in dihe:				# Search if atom a is a4 in any dihe
					a1 = i[0][0][0]
					a2 = i[0][0][1]
					a3 = i[1][0][1]
					a4 = i[1][1][1]
					a1known = 0
					a2known = 0
					a3known = 0
					a4known = 0
					if a4.index == a.index:
						for k in known_atoms:	# Check if atom is in KNOWN stack
							if a4.index == k.index:
								a4known = 1
								break
							if a1.index == k.index :
								a1known = 1
							if a2.index == k.index :
								a2known = 1
							if a3.index == k.index :
								a3known = 1
						if ((a4known == 0) and \
							(a1known == 1) and (a2known == 1) and (a3known == 1)):
								# Move everything so that a3 is in the center
								M = Translation(-1*a3.position())
								universe.applyTransformation(M)
								for k in known_atoms:
									M(k.position())
								# Solve a4
								#print 'SOLVE: ', a4
								print_zmat_dihe(i)
								A = a1.position()
								B = a2.position()
								C = a3.position()
								#print 'A = ', A, '(', a1.name , ')'
								#print 'B = ', B, '(', a2.name , ')'
								#print 'C = ', C, '(', a3.name , ')'
								R = i[1][1][2]
								#print 'R = ', R
								theta = i[1][2]
								if theta > pi/2:
									theta = pi - theta
								#print 'theta = ', theta, numpy.rad2deg(theta)
								phi = i[2]
								if phi < 0:
									phi = pi - phi
								#print 'phi = ', phi, numpy.rad2deg(phi)
								AB = B - A	# B - A
								#print 'AB = ', AB
								BC = C - B	# C - B
								#print 'BC = ', BC
								bc = BC.normal()
								#print 'bc = ', bc
								n = AB.cross(bc)
								n = n.normal()
								#print 'n = ', n
								Mx = bc
								#print 'Mx = ', Mx
								My = n.cross(bc)
								My = My.normal()
								#print 'My = ', My
								Mz = n
								#print 'Mz = ', Mz
								M = Tensor([ [Mx.x(), Mx.y(), Mx.z()],\
											 [My.x(), My.y(), My.z()],\
											 [Mz.x(), Mz.y(), Mz.z()] ])
								#print 'M: '
								#print M
								D2 = Vector(R*cos(theta), R*cos(phi)*sin(theta), R*sin(phi)*sin(theta))
								#print 'D2: '
								#print D2
								D = M*D2
								D = D + C
								#print 'D: '
								#print D
								a.setPosition(D)
								universe.atomList()[a4.index].setPosition(D)
								# Add a4 to KNOWN
								known_atoms.append(a4)
								unknown_atoms.pop(ai)
								#print '========================'
								break
	
		for k in known_atoms:
			bPDBprint_atom('HETATM', k.name, 10*k.position())
	
	
	# -----------------------------------
	
	
	def inte2cart(self, dihe):
		"""
		Reconstructs atomlist positions
		Classic
		:param universe: MMTK universe
		:type universe: MMTK universe
		:param dihe: list of dihedrals
		:type dihe: [[[],[],aval] , [[],[],aval] , dval] where [] is a bond
		"""
		known_atoms = []
		unknown_atoms = copy.deepcopy(universe.atomList())
		C = Vector(0,0,0)
		D0 = Vector(0,0,0)
		D1 = Vector(0,0,0)
		D = Vector(0,0,0)
		D2 = Vector(0,0,0)
		l = len(dihe)
		for di in range(l):	# Create the START stack
			dihe.append(self.swap_dihe(dihe[di]))
	
	
		# Assign initial values
		i = dihe[0]
		a1 = i[0][0][0]
		a2 = i[0][0][1]
		a3 = i[1][1][0]
		a4 = i[1][1][1]
		R1 = i[0][0][2]
		R2 = i[1][1][2]
		theta1 = i[0][2]
		phi = i[2]
		universe.atomList()[a1.index].setPosition(Vector(0,0,0))
		universe.atomList()[a2.index].setPosition(Vector(R1,0,0))
		universe.atomList()[a3.index].setPosition(Vector(R1-R2*cos(theta1), R2*sin(theta1), 0))
		known_atoms.append(a1)
		known_atoms.append(a2)
		known_atoms.append(a3)
	
		# Main loop
		cnt=0;
		while len(known_atoms) < len(universe.atomList()):
			cnt += 1
			if cnt > 9999:
				print 'inte2cart() cnt EXCEEDED: 9999'
				break
			ai = -1
			for a in unknown_atoms:
				ai += 1
				for i in dihe:				# Search if atom a is a4 in any dihe
					a1 = i[0][0][0]
					a2 = i[0][0][1]
					a3 = i[1][0][1]
					a4 = i[1][1][1]
					a1known = 0
					a2known = 0
					a3known = 0
					a4known = 0
					if a4.index == a.index:
						for k in known_atoms:	# Check if atom is in KNOWN stack
							if a4.index == k.index:
								a4known = 1
								break
							if a1.index == k.index :
								a1known = 1
							if a2.index == k.index :
								a2known = 1
							if a3.index == k.index :
								a3known = 1
						if ((a4known == 0) and \
							(a1known == 1) and (a2known == 1) and (a3known == 1)):
								# Move everything so that a3 is in the center
								M = Translation(-1*a3.position())
								universe.applyTransformation(M)
								for k in known_atoms:
									M(k.position())
								# Solve a4
								A = a1.position()
								B = a2.position()
								C = a3.position()
								#print 'SOLVE: ', a4
								R = i[1][1][2]
								theta = pi-i[1][2]
								phi = i[2]
								AB = B - A	# B - A
								BC = C - B	# C - B
								bc = BC.normal()
								n = AB.cross(bc)
								n = n.normal()
								D0 = C + (bc*R)
								M = Rotation(n, theta)
								D1 = M(D0)
								M = Rotation(bc, phi)
								D = M(D1)
								a.setPosition(D)
								universe.atomList()[a4.index].setPosition(D)
								
								# Add a4 to KNOWN
								known_atoms.append(a4)
								unknown_atoms.pop(ai)
								break
	
		for k in known_atoms:
			bPDBprint_atom('HETATM', k.name, 10*k.position())
	
	# -----------------------------------


class bRBSelector():
	"""
	Finds and writes rigid bodies as lists
	"""
	def __init__(self, universe, convertor):
		"""
		:param universe: an MMTK universe
		:type universe: MMTK universe
		:param bonds: lists of precomputed bonds
		:type bonds: a list of 3 element lists: [a1,a2,bval]
		"""
		#Variables
		self.bonds = convertor.bonds	#
		self.dihedrals = convertor.dihedrals
		self.convertor = convertor #
		self.cnt = 0		#checks if infinite loop
		self.bibonds = []	#an element is a bond with an index attached
		self.bibonds_dict = {}
		self.wostack = []	#stack which uses bibonds element types (work)
		self.rings = []		#stack which uses bibonds element types	(rings deposit)
		self.trash = []		#stack which uses bibonds element types (recovery)
		self.index = 0		#used to identify bonds !not the same as Atom.index
		self.inbuff = []	#indexes buffer used in inrings
		self.inrings = []	#list of sorted lists of indexes
		self.buff = []
		self.conv = 0
		self.atring = []	#list of lists of atom indexes in rings
		self.atbond = []	#list of 2-tuples atom indexes in rigid bonds
		self.atatom = []	#list of 2 tuples [atom, valence]
		self.atconj = []	#list of conjugate multiple bonds
		self.atrb = []
		self.stickbonds = []

	def __getstate__(self):
		return self.atring

	def findPaths(self):
		# Fill the START stack (bibonds)
		self.index = 0
		for i in self.bonds :
			self.index += 1
			self.bibonds.append([i, self.index])
			self.bibonds.append([convertor.swap_bond(i), self.index])
			self.bibonds_dict[self.index] = i
		# Kickstart - put the first element in WORK stack (wostack)
		self.wostack.append(self.bibonds[-1])
		self.bibonds.pop()
		# Main
		while True:
			if not self.bibonds:
				break
			elif not self.wostack:
				self.wostack.append(self.bibonds[-1])
				self.bibonds.pop()
			else:
				self.cnt += 1
				if self.cnt > (universe.numberOfAtoms()*universe.numberOfAtoms()):
					print 'findPaths(): cnt EXCEEDED!'
					break
				if not self.wostack:
					break
				else:
					found_link = 0
					for i in range(len(self.bibonds)) :
						if ( (self.wostack[-1][0][1].index == self.bibonds[i][0][0].index) and \
							 (self.wostack[-1][1] != self.bibonds[i][1])):
							found_link = 1
							self.wostack.append(self.bibonds[i])
							self.bibonds.pop(i)
							# Check for rings
							wolen = len(self.wostack)
							if(wolen >=3):
								for j in range(-3, -1 * wolen, -1) :
									#print 'wostack[',j,']:',wostack[j]
									#Check for circular path
									if ((self.wostack[-1][0][1].index == self.wostack[j][0][0].index) and\
										(self.wostack[-1][1] != self.wostack[j][1])):
										self.inbuff = []
										for k in range(wolen+j, wolen):
											self.inbuff.append(self.wostack[k][1])
										#Check if path has two way bonds
										bonds_doubled = 0
										if len(self.inbuff) > len(set(self.inbuff)) :
											bonds_doubled = 1
										#Check if ring is already in the rings stack
										self.inbuff.sort()
										ring_already_in = 0
										for t in self.inrings :
											if self.inbuff == t:
												ring_already_in = 1
												break
										#Copy a new ring to rings
										if (ring_already_in == 0) and (bonds_doubled == 0):
											for k in range(wolen+j, wolen):
												self.rings.append(self.wostack[k])
											self.inrings.append(self.inbuff)
										break
		
							break
					if found_link == 0:
						self.trash.append(self.wostack[-1])
						self.wostack.pop()

		# Prepare a list of names
		i = 0
		for m in self.inrings :
			self.buff = []
			for n in m :
				self.buff.append(self.bibonds_dict[n][0].index)
				self.buff.append(self.bibonds_dict[n][1].index)
			self.buff = set(self.buff)
			self.atring.append(list(self.buff))

		# ---------------------------------------------------------------

	def check_for_cond_rings(self, bonds, conv) :
		for mi in range(len(self.atring)) :
			for ni in range(len(self.atring)) :
				self.buff = []
				if bcontains(self.atring[mi], self.atring[ni]) == 1 :
					# Extract the hidden ring
					self.buff = list(set(self.atring[mi])-set(self.atring[ni]))
					buffbuff = []
					for o in self.atring[ni] :	# closed ring
						for p in self.buff :		# diff
							found_bond = 0
							if (p != o) :
								for r in bonds :
									if ( ((r[0].index == o) and (r[1].index == p)) or\
										 ((r[1].index == o) and (r[0].index == p)) ):
										buffbuff.append(o)
										found_bond = 1
										break
							if found_bond == 1 :
								break
					for s in list(set(buffbuff)) :
						self.buff.append(s)
					
					self.atring.append(self.buff)
					self.atring.pop(mi)
	
					c = 0
					for zi in (range(len(self.atring))) :
						c += len(self.atring[zi])
					if c == conv : # check for convergence
						return
	
					self.check_for_cond_rings( bonds, c)	#recursivity
					#break
	
	def untangleCondensed(self):
		self.conv = 0	#convergence
		self.check_for_cond_rings(self.bonds, self.conv)
		# Sort atring
		for i in self.atring:
			i.sort()
		self.atring.sort()

	def find_rings(self):
		self.findPaths()
		self.untangleCondensed()

	def print_atring(self):
		"""Prints atring list to stdout"""
		if(len(self.atring) > 1) :
			if(self.atring[0] != self.atring[1]) :
				print self.atring[0]
			for i in range(1,len(self.atring)):
				if self.atring[i] != self.atring[i-1] :
					print self.atring[i]
		elif (len(self.atring) == 1):
			print self.atring[0]
		else:
			print 'No rings found'
	#---------------------------------------------

	def print_rings(self):
		"""Prints ring atom names to stdout"""
		fin_string = ""
		fin_string += 'rings = ['
		t = 0
		if(len(self.atring) > 1) :
			if(self.atring[0] != self.atring[1]) :
				#print self.atring[0]
				if t != 0: fin_string += ', '
				t += 1
				fin_string += '['
				k = 0
				for z in self.atring[0]:
					if k != 0: fin_string += ', '
					fin_string += universe.atomList()[z].name
					k += 1
				fin_string += ']'
			for i in range(1,len(self.atring)):
				if self.atring[i] != self.atring[i-1] :
					#print self.atring[i]
					if t != 0: fin_string += ', '
					t += 1
					fin_string += '['
					k = 0
					for z in self.atring[i]:
						if k != 0: fin_string += ', '
						fin_string += universe.atomList()[z].name
						k += 1
					fin_string += ']'
		elif (len(self.atring) == 1):
			#print self.atring[0]
			if t != 0: fin_string += ', '
			t += 1
			fin_string += '['
			k = 0
			for z in self.atring[0]:
				if k != 0: fin_string += ', '
				fin_string += universe.atomList()[z].name
				k += 1
			fin_string += ']'
		else:
			pass
			#print 'No rings found'
		fin_string += ']'
		print fin_string
	#---------------------------------------------


	def growRB(self):
		"""
		Grows rigid bodies starting from
		rigid bonds (C-H, X=X) and rings
		"""
#		print 'BEGIN growRB ============='
		# Build atatom list ([atom, unsaturated rank])
		for ai in universe.atomList():
			unsat = bisSaturated(ai.type.symbol, len(ai.bondedTo()))
			if unsat > 0:
				self.atatom.append([ai.index, unsat])

		# Build atbond list (bonds with sp2 or sp atoms)
		self.atbond = []
		for ata in self.atatom:
			ring_found1 = 0
			for ri in self.atring:		# Check if ata is in any ring
				for rj in ri:
					if ata[0] == rj:
						ring_found1 = 1
						#break			# Atom 1 touches a ring
				if ring_found1 == 1:
					pass
					#break				# Checked
			for atb in self.atatom:	# Check if atb is in any ring
				ring_found2 = 0
				for ri in self.atring:
					for rj in ri:
						if atb[0] == rj:
							ring_found2 = 1
							#break	# Atom 2 touches a ring
				if (ring_found2 == 0) or (ring_found1 == 0):
					if ata.index < atb.index:
						for bo in universe.atomList()[ata[0]].bondedTo():
							if bo.index == atb[0]:
								self.atbond.append([ata[0],atb[0]])
								break


		fin_str = ""
		fin_str += 'non_ring_pi_bonds = [' # len(self.atbond)
		k = 0
		for bo in self.atbond:
			#print universe.atomList()[bo[0]].name, universe.atomList()[bo[1]].name
			if k != 0: fin_str += ', '
			fin_str += 'Bond('
			fin_str += universe.atomList()[bo[0]].name
			fin_str += ', '
			fin_str += universe.atomList()[bo[1]].name
			fin_str += ')'
			k += 1
		fin_str += ']'
		print fin_str
#			print bo
#		print '------------------'


# THE FOLLOWING CODE IS NOT NECESSARY
#		for ri in self.atring:	# Eliminates bonds sticked into rings BUG !!
#			for rj in ri:
#				le = len(self.atbond)
#				for bi in range(le):
#					if (self.atbond[bi][0] == rj) or (self.atbond[bi][1] == rj):
#						self.atbond.pop(bi)
#						break
#
#		print len(self.atbond), ': multiple bonds without sticked: '
#		for bo in self.atbond:
#			print universe.atomList()[bo[0]].name, universe.atomList()[bo[1]].name
#			#print bo
#		print '------------------'

		# Build atconj - merge multiple bonds in conjugated systems
		self.atconj = copy.deepcopy(self.atbond)
		buff = []
		cnt = 0
		merge_found = 1
		while (merge_found == 1) and (cnt < 100):
			cnt += 1
			merge_found = 0
			for coi in range(len(self.atconj)):
				for coj in range(len(self.atconj)):
					if coj != coi:
						# Check if they have a common atom (to merge them)
						for i in self.atconj[coi]:
							for j in self.atconj[coj]:
								if i == j:
									merge_found = 1
									break
							if merge_found == 1:
								break
						if merge_found == 1:	# Add coj to coi (merge)
							buff = []
							for j in self.atconj[coj]:
								common = 0
								for i in self.atconj[coi]:
									if i == j:
										common = 1
										break
								if common == 0:
									buff.append(j)
							for t in buff:
								self.atconj[coi].append(t)
							# Print the new conj
							# Remove bond coj 
							self.atconj.pop(coj) # Eliminate coj from atconj
							break
					if merge_found == 1:
						break
				if merge_found == 1:
					break


#		print 'Conjugated systems: '
#		for coi in self.atconj:
#				for k in coi:
#					print universe.atomList()[k].name, ' ',
#				print
			
#		print 'END growRB ================'


	def diffPseudoDihe(self):
		"""
		Eliminates correlated dihedrals by 
		replacing them with differences
		"""
		print 'BEGIN diffPseudoDihe() =========='
		#for ai in universe.atomList():
		#	print ai.position()

		#for bi in self.atbond:
		#	print 'atbond', universe.atomList()[bi[0]].name, universe.atomList()[bi[1]].name

		#for di in self.dihedrals:
		#	convertor.print_zmat_dihe(di)

		#Stickbonds
		self.stickbonds = []
		for bi in self.bonds:
			at1found =  -1
			at2found =  -1
			for ri in self.atring: # Search the rings
				for rj in ri:
					if (bi[0].index == rj) and \
						(universe.atomList()[bi[1].index].type.symbol != 'H'):
						at1found = rj
					elif (bi[1].index == rj)  and \
						(universe.atomList()[bi[0].index].type.symbol != 'H'):
						at2found = rj
			if ((at1found >= 0) and (at2found == -1) or\
				(at1found == -1) and (at2found >= 0)):
				self.stickbonds.append([bi[0].index, bi[1].index])
		
#		print 'Sticked bonds:'				
#		for bi in self.stickbonds:
#			print universe.atomList()[bi[0]].name, universe.atomList()[bi[1]].name


		for ri in self.atring:
			print ri
			for rj in ri:
				stick = -1
				for bi in self.stickbonds:
					if bi[0] == rj:
						stick = 1
					elif bi[1] == rj:
						stick = 0
					if stick >= 0:
						print '[',universe.atomList()[bi[0]].name,\
								universe.atomList()[bi[1]].name, '] :'
						# Search dihedrals
						for di in self.dihedrals:
							a2 = di[0][0][1]
							a3 = di[1][0][1]
							if ((a2.index == bi[0] and a3.index == bi[1]) or\
								(a2.index == bi[1] and a3.index == bi[0])):
									convertor.print_zmat_dihe(di)
						print '+++++++++++++++'
						break

		print 'END   diffPseudoDihe() =========='

#-------------------------------------------------------------------------


########
# MAIN #
########

# BUILD MOLECULE #

# "Declare" some variables
natm = 0;
atmi = 0;
x = .0;
y = .0;
z = .0;
elem = '--';
name = [];
##################

fpo1 = open(sys.argv[1], 'r');

factory = MoleculeFactory();
factory.createGroup('main');

linecnt = 0;
while True:
	line = fpo1.readline();
	if not line: break;
	linecnt += 1;
	if linecnt == 4:
		natm = int(bsubstr(line, 1, 3));
	if (linecnt >= 5) and (linecnt < (5 + natm)):   # Read cart coords
		atmi += 1;
		x = float(bsubstr(line,  2, 9));
		y = float(bsubstr(line, 12, 9));
		z = float(bsubstr(line, 22, 9));
		elem = bsubstr(line, 32, 2);
		if elem[1] == ' ':
			elem = elem[0];
		name.append(elem + str(atmi));
		factory.addAtom('main',name[atmi-1],elem);
		factory.setPosition('main',name[atmi-1],Vector\
			(x*Units.Ang, y*Units.Ang, z*Units.Ang));
	if linecnt >= (5 + natm):
		if line[0] == 'M': break
		factory.addBond('main',\
			name[int(bsubstr(line,1,3))-1], name[int(bsubstr(line,4,3))-1]);

fpo1.close();

newmol = factory.retrieveMolecule('main');
universe = InfiniteUniverse(Amber99ForceField(mod_files=['frcmod.ff99SB']));
universe.addObject(newmol);
universe.configuration();

# -------------------------------------------

#anchor(universe)	# BUG pseudo01.sdf when first atom has only one bond

#universe.writeXML(file('out.xml', 'w'));
#universe.writeToFile('out.pdb');

convertor = bConvertor(universe)
convertor.cart2all_internals()
#convertor.print_internals()

#print '\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n'
#for i in convertor.dihedrals:
#	convertor.print_zmat_dihe(i)
#convertor.inte2cart(convertor.dihedrals)

# Print rings
rb_selector = bRBSelector(universe, convertor)
rb_selector.find_rings()
#rb_selector.print_atring()
rb_selector.print_rings()
rb_selector.growRB()
#rb_selector.diffPseudoDihe()







