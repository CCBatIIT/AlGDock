# Prepares a ligand for AlGDock
# Run from: [TARGET]/ligand/AlGDock_in/

try:
  import argparse
  parser = argparse.ArgumentParser()
  parser.add_argument('ligand_mol2', default=None,
    help='Input mol2 of the ligand with sybyl atom types')
  parser.add_argument('ligand_frcmod', default=None,
    help='Input frcmod of the ligand')
  parser.add_argument('receptor_pdb', default=None,
    help='Input PDB of the receptor with AMBER atom types')
  parser.add_argument('--complex_prefix', default=None,
    help='Prefix for the complex prmtop and inpcrd files')
  parser.add_argument('--debug', action='store_true')
  args = parser.parse_args()
except:
  import sys
  class args:
    ligand_mol2 = sys.argv[1]
    ligand_frcmod = sys.argv[2]
    receptor_pdb = sys.argv[3]

import os
for FN in [args.ligand_mol2, args.ligand_frcmod, args.receptor_pdb]:
  if not os.path.isfile(FN):
    raise Exception('Input file %s not found!'%FN)

import os, inspect
dirs = {'script':os.path.dirname(os.path.abspath(\
  inspect.getfile(inspect.currentframe())))}
execfile(os.path.join(dirs['script'],'_external_paths.py'))
command_paths = findPaths(['sander'])
dirs['amber'] = os.path.abspath(os.path.dirname(command_paths['sander'])[:-4])

ligand_prefix = os.path.basename(args.ligand_mol2)[:-5]
receptor_prefix = os.path.basename(args.receptor_pdb)[:-4]
if args.complex_prefix is None:
  complex_prefix = ligand_prefix + '-' + receptor_prefix
else:
  complex_prefix = args.complex_prefix

if not os.path.isdir(ligand_prefix):
  os.makedirs(ligand_prefix)
os.chdir(ligand_prefix)
if not os.path.isfile(os.path.join('..',ligand_prefix+'.mol2')):
  print '\n*** Writing mol2 file with amber atom types ***'
  command = dirs['amber']+'/bin/antechamber' + \
    ' -i {0} -fi mol2 -o {1}.mol2 -fo mol2 -rn {2}'.format(\
      args.ligand_mol2,os.path.join('..',ligand_prefix),ligand_prefix)
  os.system(command)
os.chdir('..')
if not args.debug:
  os.system('rm -rf '+ligand_prefix)

if not (os.path.isfile(os.path.join('..',complex_prefix+'.prmtop')) and \
        os.path.isfile(os.path.join('..',complex_prefix+'.inpcrd'))):
  print '\n*** Generating prmtop and inpcrd files ***'
  tleap_F = open(complex_prefix+'.tleap','w')
  tleap_F.write("""
source leaprc.ff14SB
set default PBRadii bondi

# Receptor
receptor = loadpdb {0}

# Ligand
source leaprc.gaff
loadamberparams {1}
ligand = loadmol2 {2}.mol2
saveoff ligand {2}.lib
loadoff {2}.lib

# Complex 
complex = combine {{receptor, ligand}}
saveamberparm complex {3}.prmtop {3}.inpcrd

quit
""".format(args.receptor_pdb, args.ligand_frcmod, ligand_prefix, complex_prefix))
  tleap_F.close()
  command = dirs['amber']+'/bin/tleap -f {0}.tleap'.format(complex_prefix)
  os.system(command)

if args.debug:
  if os.path.isfile('leap.log'):
    os.rename('leap.log', complex_prefix+'.leaplog')
else:
  for FN in [complex_prefix+'.tleap',ligand_prefix+'.mol2',ligand_prefix+'.lib', \
    'leap.log', \
    'ANTECHAMBER_AC.AC', 'ANTECHAMBER_AC.AC0', \
    'ANTECHAMBER_BOND_TYPE.AC', 'ANTECHAMBER_BOND_TYPE.AC0', 'ATOMTYPE.INF']:
    if os.path.isfile(FN):
      os.remove(FN)
