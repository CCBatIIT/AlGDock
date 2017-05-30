# Converts AMBER prmtop (and inpcrd) files to a MMTK database,
# which is basically a python script defining a molecule's
# atoms, bonds, and default coordinates

import os, sys, time, numpy

#################
### CONSTANTS ###
#################

# Map flag to variable names
varnames = ['POINTERS','TITLE','ATOM_NAME','AMBER_ATOM_TYPE','CHARGE','MASS',\
            'NONBONDED_PARM_INDEX','LENNARD_JONES_ACOEF','LENNARD_JONES_BCOEF',\
            'ATOM_TYPE_INDEX','BONDS_INC_HYDROGEN','BONDS_WITHOUT_HYDROGEN',\
            'RADII','SCREEN']
# Masses for atoms defined in the GAFF force field
mass2symbols = {12.01:'C',1.008:'H',19.00:'F',35.45:'Cl',79.90:'Br',126.9:'I',\
                14.01:'N',16.00:'O',30.97:'P',32.06:'S', 65.4:'Zn'}

############
### MAIN ###
############

import argparse
parser = argparse.ArgumentParser(
  description='Convert AMBER prmtop and inpcrd files to a MMTK database file')
parser.add_argument('prmtop_FN', help='AMBER prmtop file')
parser.add_argument('db_FN', help='MMTK Database file')
parser.add_argument('--inpcrd_FN', help='AMBER inpcrd file')
parser.add_argument('-f', help='Does nothing')
args = parser.parse_args()

print "Creating database "+args.db_FN

### Loads AMBER parameter file
import AlGDock
import AlGDock.IO
prmtop_IO = AlGDock.IO.prmtop()
prmtop = prmtop_IO.read(args.prmtop_FN, varnames)

# Modify atom names to be acceptable python variables and include the atom index
prmtop['ATOM_NAME'] = ['%si%s'%(a_n.replace("'","p").strip(),ind) \
  for (a_n,ind) in zip(prmtop['ATOM_NAME'],range(len(prmtop['ATOM_NAME'])))]

NATOM = prmtop['POINTERS'][0]
NTYPES = prmtop['POINTERS'][1]

### Extract Lennard-Jones well depth and radii for each atom
LJ_radius = numpy.ndarray(shape=(NTYPES), dtype=float)
LJ_depth = numpy.ndarray(shape=(NTYPES), dtype=float)
for i in range(NTYPES):
  LJ_index = prmtop['NONBONDED_PARM_INDEX'][NTYPES*i+i]-1
  if prmtop['LENNARD_JONES_ACOEF'][LJ_index]<1.0e-6:
    LJ_radius[i] = 0
    LJ_depth[i] = 0
  else:
    factor = 2 * prmtop['LENNARD_JONES_ACOEF'][LJ_index] / \
      prmtop['LENNARD_JONES_BCOEF'][LJ_index]
    LJ_radius[i] = pow(factor, 1.0/6.0) * 0.5
    LJ_depth[i] = prmtop['LENNARD_JONES_BCOEF'][LJ_index] / 2 / factor
# More useful for later calculations
root_LJ_depth = numpy.sqrt(LJ_depth)
LJ_diameter = LJ_radius*2
del i, LJ_index, factor

### Loads AMBER inpcrd file
if (args.inpcrd_FN is not None) and os.path.isfile(args.inpcrd_FN):
  inpcrdF = open(args.inpcrd_FN,'r')
  inpcrd = inpcrdF.read().split('\n')
  inpcrdF.close()
  inpcrd.pop(0) # Title
  NATOM = int(inpcrd.pop(0)) # Number of atoms
  w = 12 # Width of field
  coords = []
  for line in inpcrd:
        coords = coords + [float(line[x:x+w]) for x in range(0,len(line),w)]
  del w, inpcrdF, inpcrd
else:
  coords = None

### Writes database
db_dir = os.path.dirname(args.db_FN)
if not (db_dir=='' or os.path.exists(db_dir)):
  os.system('mkdir -p '+db_dir)

db = open(args.db_FN,'w')

db.write("name='%s'\n"%os.path.basename(args.db_FN))

for [name,mass] in zip(prmtop['ATOM_NAME'],prmtop['MASS']):
  if mass in mass2symbols.keys():
    db.write(name.strip()+" = Atom('"+mass2symbols[mass]+"')\n")
  else:
    raise Exception('Unknown atom with mass: %f!'%mass)
    sys.exit()

# The order of atoms in the prmtop file (not used by MMTK)
db.write("prmtop_order = [" + ", ".join(["%s"%name.strip() \
  for name in prmtop['ATOM_NAME']]) + "]\n")

bondList = list(prmtop['BONDS_INC_HYDROGEN']) + \
  list(prmtop['BONDS_WITHOUT_HYDROGEN'])
db.write("bonds = [" + ", ".join(["Bond(%s, %s)"%(\
  prmtop['ATOM_NAME'][bondList[i]/3].strip(),prmtop['ATOM_NAME'][bondList[i+1]/3].strip()) \
  for i in range(0,len(bondList),3)]) + "]\n")

db.write("amber12_atom_type = {" + ", ".join(["%s: '%s'"%(name.strip(),type.strip()) \
  for (name,type) in zip(prmtop['ATOM_NAME'],prmtop['AMBER_ATOM_TYPE'])]) + "}\n")

# Write the charge, converted to units of electric charge
# AMBER prmtop files multiply the actual charge by 18.2223, hence the division
db.write("amber_charge = {" + ", ".join(\
  ["%s: '%f'"%(name.strip(),charge/18.2223) \
   for (name,charge) in zip(prmtop['ATOM_NAME'],prmtop['CHARGE'])]) + "}\n")

# Write the grid scaling factors
# Because the grids are in units of kcal/mol, the scaling factors are multiplied by 4.184 to convert to kJ/mol
db.write("scaling_factor_electrostatic = {" + \
  ", ".join(["%s: %f"%(name.strip(),4.184*charge/18.2223) 
    for (name,charge) in zip(prmtop['ATOM_NAME'],prmtop['CHARGE'])]) \
  + "}\n")
atom_type_indicies = [prmtop['ATOM_TYPE_INDEX'][atom_index]-1 \
  for atom_index in range(NATOM)]
db.write("scaling_factor_LJr = {" + \
  ", ".join(["%s: %f"%(name.strip(),4.184*root_LJ_depth[type_index]*(LJ_diameter[type_index]**6)) \
             for (name,type_index) in zip(prmtop['ATOM_NAME'],atom_type_indicies)]) + "}\n")
db.write("scaling_factor_LJa = {" + \
  ", ".join(["%s: %f"%(name.strip(),4.184*root_LJ_depth[type_index]*(LJ_diameter[type_index]**3)) \
             for (name,type_index) in zip(prmtop['ATOM_NAME'],atom_type_indicies)]) + "}\n")

# Write the generalized Born implicit solvent parameters.
# The radii will be converted from Angstroms to nanometers.
db.write("scaling_factor_BornRadii = {" + \
  ", ".join(["%s: %f"%(name.strip(),radius/10.) \
             for (name,radius) in zip(prmtop['ATOM_NAME'],prmtop['RADII'])]) + "}\n")
db.write("scaling_factor_BornScreening = {" + \
  ", ".join(["%s: %f"%(name.strip(),radius) \
             for (name,radius) in zip(prmtop['ATOM_NAME'],prmtop['SCREEN'])]) + "}\n")

# Write the coordinates, converted from Angstroms to nanometers
if coords is not None:
  db.write("configurations = {\n")
  db.write("'default': Cartesian({" + ", ".join(["%s: (%f, %f, %f)"%(name.strip(), coords[i]/10.0, coords[i+1]/10.0, coords[i+2]/10.0) for (name,i) in zip(prmtop['ATOM_NAME'],range(0,NATOM*3,3))]) + "})}\n")

db.close()
