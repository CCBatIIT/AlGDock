# Prepares a ligand for AlGDock

try:
  import argparse
  parser = argparse.ArgumentParser()
  parser.add_argument('in_FN', default=None,
    help='Input mol2 with sybyl atom types')
  parser.add_argument('out_FN', default=None,
    help='Output tarball that includes prmtop and inpcrd files')
  parser.add_argument('--debug', action='store_true')
  args = parser.parse_args()
except:
  import sys
  class args:
    in_FN = sys.argv[1]
    out_FN = sys.argv[2]

import os
args.in_FN = os.path.abspath(args.in_FN)
args.out_FN = os.path.abspath(args.out_FN)

if not os.path.isfile(args.in_FN):
  raise Exception('Input file not found!')

import os, inspect
dirs = {'script':os.path.dirname(os.path.abspath(\
  inspect.getfile(inspect.currentframe())))}
execfile(os.path.join(dirs['script'],'_external_paths.py'))
command_paths = findPaths(['sander'])
dirs['amber'] = os.path.abspath(os.path.dirname(command_paths['sander'])[:-4])
dirs['temp'] = args.out_FN + '.tmp'

prefix = '.'.join(os.path.dirname(args.out_FN).split('/')[-1].split('.')[:-1]) \
  + '.' + os.path.basename(args.out_FN[:-7])
if not os.path.isdir(dirs['temp']):
  os.system('mkdir -p '+dirs['temp'])
os.chdir(dirs['temp'])

if not os.path.isfile(prefix+'.mol2'):
  print '\n*** Writing mol2 file with amber atom types ***'
  command = dirs['amber']+'/bin/antechamber' + \
    ' -i {0} -fi mol2 -o {1}.mol2 -fo mol2 -rn {1}'.format(\
      args.in_FN, prefix)
  os.system(command)

if not os.path.isfile(prefix+'.frcmod'):
  print '\n*** Generating frcmod file ***'
  command = dirs['amber']+'/bin/parmchk' +\
    ' -i {0}.mol2 -f mol2 -o {0}.frcmod -a Y -w Y'.format(prefix)
  os.system(command)

if not (os.path.isfile(os.path.join(prefix+'.prmtop')) and \
        os.path.isfile(os.path.join(prefix+'.inpcrd'))):
  print '\n*** Generating prmtop and inpcrd files ***'
  tleap_F = open(prefix+'.tleap','w')
  tleap_F.write("""
source leaprc.gaff2
loadamberparams {0}.frcmod

set default PBRadii mbondi2

ligand = loadmol2 {0}.mol2
saveoff ligand {0}.lib
loadoff {0}.lib

saveamberparm ligand {0}.prmtop {0}.inpcrd

quit
""".format(prefix))
  tleap_F.close()
  command = dirs['amber']+'/bin/tleap -f {0}.tleap'.format(prefix)
  os.system(command)
  if not (os.path.isfile(os.path.join(prefix+'.prmtop')) and \
          os.path.isfile(os.path.join(prefix+'.inpcrd'))):
    raise Exception('Unable to create prmtop and inpcrd for '+prefix)

if args.debug:
  if os.path.isfile('leap.log'):
    os.rename('leap.log',prefix+'.leaplog')
else:
  for FN in [prefix+'.tleap',prefix+'.mol2',prefix+'.lib','leap.log', \
    'ANTECHAMBER_AC.AC', 'ANTECHAMBER_AC.AC0', \
    'ANTECHAMBER_BOND_TYPE.AC', 'ANTECHAMBER_BOND_TYPE.AC0', 'ATOMTYPE.INF']:
    if os.path.isfile(FN):
      os.remove(FN)

db_FN = prefix.lower()+'.db'
if not os.path.isfile((db_FN)):
  print '\n*** Generating MMTK database ***'
  command = 'python '+dirs['script']+'/prmtop2database.py' + \
    ' {0}.prmtop {1}.db --inpcrd_FN {0}.inpcrd'.format(prefix,prefix.lower())
  os.system(command)
  if not os.path.isfile(db_FN):
    raise Exception('Unable to create db for '+prefix)

# Tests the ligand in MMTK
import MMTK

MMTK.Database.molecule_types.directory = os.getcwd()

from MMTK.ForceFields import Amber12SBForceField
ff = Amber12SBForceField(mod_files=[prefix+'.frcmod'])
molecule = MMTK.Molecule(db_FN)
universe = MMTK.Universe.InfiniteUniverse()
universe.addObject(molecule)
universe.setForceField(ff)
print 'MMTK Energy: %f'%universe.energy()

# Compresses the ligand files in a tarball
import tarfile
tarF = tarfile.open(args.out_FN,'w:gz')
tarF_contents = [prefix+'.prmtop', prefix+'.inpcrd', \
  prefix+'.frcmod', prefix.lower()+'.db']
for FN in tarF_contents:
  tarF.add(FN)
  os.remove(FN)
tarF.close()

os.chdir('..')
if not args.debug:
  os.system('rm -rf '+dirs['temp'])
