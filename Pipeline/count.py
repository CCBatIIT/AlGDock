# To be run from the [TARGET] directory, or
# the parent directory with the multi_target option
import os, sys
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

def loopISM(ism_FNs, search_string, field_width):
  outstr = ''
  for ism_FN in ism_FNs:
    prefix = os.path.basename(ism_FN)[:-4]
    outstr += '%s\t%d\n'%(prefix.rjust(field_width), \
      num_nonzero(search_string%prefix))
  outstr += '%s\t%d\n'%('Total'.rjust(field_width), \
    num_nonzero(search_string%''))
  return outstr
  
if sys.argv[-1] == '--multi_target':
  dirs = sorted([d for d in glob.glob('*') if os.path.isdir(d)])
  print 'Parsing multiple targets'
else:
  dirs = ['.']
  
for dir in dirs:
  os.chdir(dir)

  print '-'*20
  print os.getcwd()
  
  nligands_all = 0
  nligands_dock = 0
  nligands_AlGDock = 0
  ism_FNs = glob.glob('ligand/*.ism')
  if ism_FNs==[]:
    print 'No ligands found for '+dir
    os.chdir('..')
    continue

  print "\nLigands"
  first_field_width = int(np.max([len(ism_FN[:-4]) for ism_FN in ism_FNs]))
  print ' '*first_field_width+"\tSMILES\tDock\tAlGDock"
  for ism_FN in ism_FNs:
    prefix = os.path.basename(ism_FN)[:-4]
    ism_nligands_all = countLines(ism_FN)
    ism_nligands_dock = num_nonzero('ligand/dock_in/%s*/*.mol2'%prefix)
    ism_nligands_AlGDock = num_nonzero('ligand/AlGDock_in/%s*/*.tar.gz'%prefix)
    nligands_all += ism_nligands_all
    nligands_dock += ism_nligands_dock
    nligands_AlGDock += ism_nligands_AlGDock
    print '%*s\t%d\t%d\t%d'%(first_field_width,\
      prefix,ism_nligands_all,ism_nligands_dock,ism_nligands_AlGDock)
  print '%*s\t%d\t%d\t%d'%(first_field_width,\
    'Total',nligands_all,nligands_dock,nligands_AlGDock)

  print "\nReceptors"
  print "\tHomology  ", num_nonzero('receptor/3-models/pdb_noH/*.pdb')
  print "\tDock      ", num_nonzero('receptor/dock_in/*.sph')
  print "\tAMBER     ", num_nonzero('receptor/amber_in/*.prmtop')
  print "\tAlGDock   ", num_nonzero('receptor/AlGDock_in/*.PB.nc')

  print "\nComplexes for AlGDock"
  print loopISM(ism_FNs, 'complex/AlGDock_in/%s*/*/*.tar.gz', first_field_width)

  print "\nCompleted Calculations"
  print "for dock"
  print loopISM(ism_FNs, 'dock6/%s*/*/*.nc', first_field_width)
  print "for AlGDock"
  print loopISM(ism_FNs, 'AlGDock/dock/%s*/*/*/f_RL.pkl.gz', first_field_width)

  os.chdir('..')
