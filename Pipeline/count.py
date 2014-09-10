# To be run from the [TARGET] directory
import os
import glob
import numpy as np

def countLines(FN):
  if os.path.isfile(FN):
    F = open(FN,'r')
    nlines = len(F.read().strip().split('\n'))
    F.close
    return nlines
  else:
    raise Exception(FN+' is not a file!')

def num_nonzero(path):
  FNs = [FN for FN in glob.glob(path) if os.path.getsize(FN)>0]
  return len(FNs)

print "\nLigands"
nligands_all = 0
nligands_dock = 0
nligands_AlGDock = 0
ism_FNs = glob.glob('ligand/*.ism')
first_field_width = int(np.max([len(ism_FN[:-4]) for ism_FN in ism_FNs]))
print ' '*first_field_width+"\tDock\tAlGDock\tSMILES"
for ism_FN in ism_FNs:
  prefix = os.path.basename(ism_FN)[:-4]
  ism_nligands_all = countLines(ism_FN)
  ism_nligands_dock = num_nonzero('ligand/dock_in/%s*.mol2'%prefix)
  ism_nligands_AlGDock = num_nonzero('ligand/AlGDock_in/%s*.prmtop'%prefix)
  nligands_all += ism_nligands_all
  nligands_dock += ism_nligands_dock
  nligands_AlGDock += ism_nligands_AlGDock
  print '%*s\t%d\t%d\t%d'%(first_field_width,\
    prefix,ism_nligands_dock,ism_nligands_AlGDock,ism_nligands_all)
print '%*s\t%d\t%d\t%d'%(first_field_width,\
  'Total',nligands_dock,nligands_AlGDock,nligands_all)

print "\nReceptors"
print "\tHomology: ", num_nonzero('receptor/3-models/pdb_noH/*.pdb')
print "\tDock:     ", num_nonzero('receptor/dock_in/*.nrg')
print "\tAMBER:    ", num_nonzero('receptor/amber_in/*.prmtop')
print "\tAlGDock:  ", num_nonzero('receptor/AlGDock_in/*.PB.nc')

print "\nComplexes for AlGDock: ", \
  num_nonzero('complex/AlGDock_in/*.prmtop.gz')

print "\nCompleted Calculations"
print "\tDock:     ", num_nonzero('dock6/*.mol2.gz')
print "\tAlGDock:  ", num_nonzero('AlGDock/dock/*/f_RL.pkl.gz')

