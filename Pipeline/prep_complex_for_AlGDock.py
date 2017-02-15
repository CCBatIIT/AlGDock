# Prepares a complex for AlGDock

try:
  import argparse
  parser = argparse.ArgumentParser()
  parser.add_argument('ligand_mol2', default=None,
    help='Input mol2 of the ligand with sybyl atom types')
  parser.add_argument('receptor_pdb', default=None,
    help='Input PDB of the receptor with AMBER atom types')
  parser.add_argument('complex_tarball', default=None,
    help='Prefix for the complex prmtop and inpcrd files')
  parser.add_argument('--debug', action='store_true')
  args = parser.parse_args()
except:
  import sys
  class args:
    ligand_mol2 = sys.argv[1]
    receptor_pdb = sys.argv[2]

import os
for FN in [args.ligand_mol2, args.receptor_pdb]:
  if not os.path.isfile(FN):
    raise Exception('Input file %s not found!'%FN)
args.ligand_mol2 = os.path.abspath(args.ligand_mol2)
args.receptor_pdb = os.path.abspath(args.receptor_pdb)
args.complex_tarball = os.path.abspath(args.complex_tarball)

import os, inspect
dirs = {'script':os.path.dirname(os.path.abspath(\
  inspect.getfile(inspect.currentframe())))}
execfile(os.path.join(dirs['script'],'_external_paths.py'))
command_paths = findPaths(['sander'])
dirs['amber'] = os.path.abspath(os.path.dirname(command_paths['sander'])[:-4])
dirs['temp'] = args.complex_tarball + '.tmp'

if not os.path.isdir(dirs['temp']):
  os.system('mkdir -p '+dirs['temp'])
os.chdir(dirs['temp'])

ligand_prefix = '.'.join(os.path.dirname(args.ligand_mol2).split('/')[-1].split('.')[:-1]) \
  + '.' + os.path.basename(args.ligand_mol2)[:-5]
# The receptor file name ends with '.pdb2pqr_amber.pqr',
# which is 18 characters long
receptor_prefix = os.path.basename(args.receptor_pdb)[:-18]
complex_prefix = ligand_prefix + '-' + receptor_prefix

if not os.path.isfile(ligand_prefix+'.mol2'):
  print '\n*** Writing mol2 file with amber atom types ***'
  command = dirs['amber']+'/bin/antechamber' + \
    ' -i {0} -fi mol2 -o {1}.mol2 -fo mol2 -rn LIG'.format(\
      args.ligand_mol2,ligand_prefix)
  os.system(command)
  if not os.path.isfile(ligand_prefix+'.mol2'):
    print command
    raise Exception('Could not write mol2 file')

if not os.path.isfile(ligand_prefix+'.frcmod'):
  print '\n*** Generating frcmod file ***'
  command = dirs['amber']+'/bin/parmchk' +\
    ' -i {0}.mol2 -f mol2 -o {0}.frcmod -a Y -w Y'.format(ligand_prefix)
  os.system(command)

if not (os.path.isfile(os.path.join(complex_prefix+'.prmtop')) and \
        os.path.isfile(os.path.join(complex_prefix+'.inpcrd')) and \
        os.path.isfile(os.path.join(complex_prefix+'.pdb'))):
  print '\n*** Generating prmtop and inpcrd and pdb files ***'
  tleap_F = open(complex_prefix+'.tleap','w')
  tleap_F.write("""
source leaprc.ff14SB
set default PBRadii mbondi2

# Receptor
receptor = loadpdb {0}

# Ligand
source leaprc.gaff2
loadamberparams {1}.frcmod
ligand = loadmol2 {1}.mol2
saveoff ligand {1}.lib
loadoff {1}.lib

# Complex 
complex = combine {{receptor, ligand}}
saveamberparm complex {2}.prmtop {2}.inpcrd
savepdb complex {2}.pdb

quit
""".format(args.receptor_pdb, ligand_prefix, complex_prefix))
  tleap_F.close()
  command = dirs['amber']+'/bin/tleap -f {0}.tleap'.format(complex_prefix)
  os.system(command)

if os.path.isfile(os.path.join(complex_prefix+'.pdb')):
  print '\n*** Setting fixed atoms in pdb file ***'
  command = 'python {0}/label_fixed_atoms.py {1}'
  command = command.format(dirs['script'], os.path.join(complex_prefix+'.pdb'))
  os.system(command)

# Compresses the complex files in a tarball
import tarfile
tarF = tarfile.open(args.complex_tarball,'w:gz')
tarF_contents = [complex_prefix+'.'+ext for ext in ['prmtop', 'inpcrd', 'pdb']]
for FN in tarF_contents:
  tarF.add(FN)
tarF.close()

os.chdir('..')

if not args.debug:
  os.system('rm -rf '+dirs['temp'])

