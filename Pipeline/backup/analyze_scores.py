# To be run from the [TARGET]/ directory

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--dock6', default='dock6/', \
  help='Directory with dock6 results')
parser.add_argument('--AlGDock', default='AlGDock/dock/', \
  help='Directory with AlGDock results')
parser.add_argument('--ligand', default='ligand/', \
  help='Directory with ligand libraries')
parser.add_argument('--analysis', default='analysis/', \
  help='Directory with analysis results')
args = parser.parse_args()

# Check for the existence of directories
import os
if not os.path.isdir(args.analysis):
  os.system('mkdir -p '+args.analysis)

# Determine the library names
import glob
libraries = [os.path.basename(FN)[:-4] \
  for FN in glob.glob(os.path.join(args.ligand,'*.ism'))] + ['']

# Extract and sort dock6 scores.
# Store them in separate libraries and all together.
import numpy as np
if not np.array([os.path.isfile(os.path.join(\
                  args.analysis,key+'_dock6_scores.txt')) \
                    for key in libraries]).all():
  mol2_FNs = glob.glob(os.path.join(args.dock6,'*.mol2.gz'))
  mol2_FNs = [FN for FN in mol2_FNs if os.path.getsize(FN)>0]
  
  import gzip
  dock6_scores = []
  for mol2_FN in mol2_FNs:
    F = gzip.open(mol2_FN,'r')
    lines = [F.readline() for a in range(4)]
    F.close()
    lines.pop(0) # Blank
    lines.pop(1) # Cluster size
    if lines==['', '']:
      os.remove(mol2_FN)
      continue
    dock6_scores.append([mol2_FN[6:-8]]+[line.split()[-1] for line in lines])

  dock6_scores.sort(key=lambda record: float(record[2]))

  for key in libraries:
    F = open(os.path.join(args.analysis, key+'_dock6_scores.txt'),'w')
    F.write('\n'.join(['\t'.join([field for field in score]) for score in dock6_scores if score[0].find(key)!=-1]))
    F.close()
else:
  F = open(os.path.join(args.analysis, '_dock6_scores.txt'),'r')
  dock6_scores = [line.split('\t') for line in F.read().strip().split('\n')]
  F.close()

