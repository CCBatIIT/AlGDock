# Expands a SMILES string to a 3D structure in mol2 format with am1bcc charges

try:
  import argparse
  parser = argparse.ArgumentParser()
  parser.add_argument('out_FN', default='AAA.mol2', help='Output file name')
  parser.add_argument('SMILES', default=None, help='SMILES string')
  args = parser.parse_args()
except:
  import sys
  class args:
    out_FN = sys.argv[1]
    SMILES = sys.argv[2]

print 'Generating %s based on %s'%(args.out_FN, args.SMILES)

import os, inspect
dirs = {}
dirs['script'] = os.path.dirname(os.path.abspath(\
  inspect.getfile(inspect.currentframe())))
execfile(os.path.join(dirs['script'],'_external_paths.py'))
command_paths = findPaths(['balloon','chimera'])
balloon_FN = os.path.join('balloon', os.path.basename(args.out_FN))

for dirN in [os.path.dirname(args.out_FN),'balloon']:
  if (dirN!='') and not os.path.isdir(dirN):
    os.system('mkdir -p '+dirN)

if not (os.path.isfile(balloon_FN) or \
        os.path.isfile(args.out_FN)):
  # Run Balloon to convert from a SMILES string to a 3D structure
  MMFF94_FN = os.path.join(os.path.dirname(command_paths['balloon']),'MMFF94.mff')
  command = command_paths['balloon'] + ' -f ' + MMFF94_FN + \
    ' --nconfs 1 --nGenerations 300 "' + args.SMILES + '" ' + balloon_FN
  os.system(command)

if os.path.isfile(os.path.basename(balloon_FN)[:-5]+'_bad.mol2'):
  print 'Balloon failed!'
  # Make the final output and empty file
  open(args.out_FN, 'a').close()

if os.path.isfile(balloon_FN):
  # Select the first model from the mol2 file
  F = open(balloon_FN, "r")
  mol2 = F.read()
  F.close()

  if mol2.count("@<TRIPOS>MOLECULE")>1:
    print 'Keeping first configuration in'+balloon_FN
    confs = mol2.strip().split("@<TRIPOS>MOLECULE")
    if confs[0]=='':
      confs.pop(0)
    F = open(balloon_FN,"w")
    F.write("@<TRIPOS>MOLECULE"+confs[0])
    F.close()

if not os.path.isfile(args.out_FN):
  # Run chimera to get AM1BCC charges
  prep_script = os.path.join(dirs['script'], '_prep_ligand.chimera.py')
  command = command_paths['chimera'] + " --nogui --script" + \
    " '%s --in_FN %s --out_FN %s'"%(prep_script, balloon_FN, args.out_FN)
  print command
  os.system(command)
if os.path.isfile(balloon_FN) and \
   os.path.isfile(args.out_FN):
  os.remove(balloon_FN)