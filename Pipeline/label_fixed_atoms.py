# labels the occupancy column with 1.0 for the receptor and 0.0 for the ligand

import os, sys, glob

import argparse
parser = argparse.ArgumentParser(description='Labels fixed atoms in a pdb file')
parser.add_argument('pdb', help='pdb file')
args = parser.parse_args()

F = open(args.pdb,'r')
lines = F.read().strip().split('\n')
F.close()

newlines = []
for line in lines:
  if len(line)<55:
    newlines.append(line)
  else:
    newlines.append('%s %3.2f  %3.2f'%(line[:55], \
        (not line[17:20]=='LIG'), 0.0))

F = open(args.pdb,'w')
F.write('\n'.join(newlines))
F.close()
