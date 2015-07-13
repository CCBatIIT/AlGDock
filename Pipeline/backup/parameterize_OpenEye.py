# Generate AMBER parameters for an input file

import os, sys

startdir = os.getcwd()

import openeye as oe

genConformers = oe.oeomega.OEOmega()
genConformers.SetMaxConfs(20)
genConformers.SetStrictStereo(False)

# For testing
sys.argv.append('/home/pgupta13/prep_ligand/protomers_2J79_GTL.smi')
sys.argv.append('protomers_2J79_GTL')

############
### Main ###
############

import argparse
parser = argparse.ArgumentParser(description='Generate AMBER parameters')
parser.add_argument('inF',help='input file')
parser.add_argument('prefix', help='output prefix')
args = parser.parse_args()

if os.path.exists(args.prefix+'.prmtop'):
  raise Exception('Ligand %s already parameterized'%(args.prefix))

ligand = oe.oechem.OEMol()

ifs = oe.oechem.oemolistream()
ofs = oe.oechem.oemolostream()

print 'Parameterizing '+args.inF
ifs.open(args.inF)
oe.oechem.OEReadMolecule(ifs,ligand)
oe.oechem.OECanonicalOrderAtoms(ligand)

# Save old conformer
ofs.open(args.prefix+'.original.mol2')
oe.oechem.OEWriteMolecule(ofs,ligand)
ofs.close()

# Generates conformers for partial charge calculations
print 'Generating new conformers'
ligand.DeleteConfs()
ligand.NewConf()
genConformers(ligand)

# Calculates partial charges with AM1BCC
if not os.path.exists(args.prefix+'.AM1BCC.mol2'):
  print 'Assigning AM1BCC charges'
  oe.oequacpac.OEAssignPartialCharges(ligand,oe.oequacpac.OECharges_AM1BCC)
  # Output mol2 files
  oe.oechem.OECanonicalOrderAtoms(ligand)
  # Get lowest energy conformer
  conf_energy = 1E10
  for current_conf in ligand.GetConfs():
    current_energy = current_conf.GetEnergy()
    if current_energy<conf_energy:
      conf = current_conf
      conf_energy = current_energy
  # Write the lowest energy conformer to a mol2 file
  ofs.open(args.prefix+'.AM1BCC.mol2')
  oe.oechem.OEWriteMolecule(ofs,conf)
  ofs.close()

#


# Check to make sure charging went okay
if not os.path.exists(args.prefix+'.AM1BCC.mol2'):
  raise Exception('AM1BCC charged mol2 file not found!')

AM1BCC_F = open(args.prefix+'.AM1BCC.mol2','r')
AM1BCC = AM1BCC_F.read()
AM1BCC_F.close()
if AM1BCC.find('NO_CHARGE')>-1:
  os.remove(args.prefix+'.AM1BCC.mol2')
  raise Exception('AM1BCC charging failed')

# Get net charge
AM1BCC = AM1BCC[AM1BCC.find('@<TRIPOS>ATOM')+14:]
AM1BCC = AM1BCC[:AM1BCC.find('@<TRIPOS>')-1].split('\n')
nc = 0.
for line in AM1BCC:
  nc = nc + float(line.split()[-1])
nc = int(nc)

# Prepare AMBER input files
tmpdir = os.path.join('/scratch','dm225',args.prefix)
os.system('mkdir -p '+tmpdir)
os.chdir(tmpdir)

command = '$AMBERHOME/bin/antechamber -i {0}/{1}.AM1BCC.mol2 -fi mol2 -o {1}.mol2 -fo mol2 -s 0 -rn {1} -nc {2}'.format(startdir,args.prefix,nc)
print 'Running ANTECHAMBER with '+command
os.system(command)

for FN in ['ANTECHAMBER_AC.AC','ANTECHAMBER_AC.AC0','ANTECHAMBER_BOND_TYPE.AC','ANTECHAMBER_BOND_TYPE.AC0','ATOMTYPE.INF']:
  if os.path.exists(FN):
    os.remove(FN)

command = '$AMBERHOME/bin/parmchk -i {0}.mol2 -f mol2 -o {0}.frcmod'.format(args.prefix)
print 'Running parmchk with '+command
os.system(command)

tleapF = open('{0}.tleap'.format(args.prefix),'w')
tleapF.write('''
source leaprc.gaff
loadamberparams {0}.frcmod

set default PBRadii bondi

ligand = loadmol2 {0}.mol2
saveoff ligand {0}.lib
loadoff {0}.lib

saveamberparm ligand {0}.prmtop {0}.inpcrd
savemol2 ligand {0}.mol2
savepdb ligand {0}.pdb

quit
'''.format(args.prefix))
tleapF.close()
print 'Running tleap'
os.system('$AMBERHOME/bin/tleap -f {0}.tleap > {0}.log'.format(args.prefix))

os.system('cp {0}/{2}.prmtop {1}/{2}.prmtop'.format(tmpdir,startdir,args.prefix))
os.system('cp {0}/{2}.inpcrd {1}/{2}.inpcrd'.format(tmpdir,startdir,args.prefix))
os.system('cp {0}/{2}.mol2 {1}/{2}.mol2'.format(tmpdir,startdir,args.prefix))
os.system('cp {0}/{2}.frcmod {1}/{2}.frcmod'.format(tmpdir,startdir,args.prefix))
os.system('cp {0}/{2}.lib {1}/{2}.lib'.format(tmpdir,startdir,args.prefix))
os.system('rm -rf {0}'.format(tmpdir))

print '---'

