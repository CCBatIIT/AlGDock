# To be run from the [TARGET] directory, or
# the parent directory with the multi_target option

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--clean_dock_not_complex', action='store_true', \
  default=False, \
  help='Removes docking results if there is not an AMBER complex')
parser.add_argument('--multi_target', action='store_true', default=False, \
  help='Performs counts for multiple targets')
args = parser.parse_args()

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

def numfiles(path, nonzero=True):
  FNs = glob.glob(path)
  if nonzero:
    FNs = [FN for FN in FNs if os.path.getsize(FN)>0]
  return len(FNs)

def FN_in_first_not_second(path1,path2):
  ext1 = path1.split('*')[-1]
  ext2 = path2.split('*')[-1]
  set1 = [os.path.basename(FN)[:-len(ext1)] \
    for FN in glob.glob(path1) if os.path.getsize(FN)>0]
  set2 = [os.path.basename(FN)[:-len(ext2)] \
    for FN in glob.glob(path2) if os.path.getsize(FN)>0]
  return [x for x in set1 if x not in set2]

def complex_nonoverlap(path1,path2):
  def complex_ids(path):
    ext = path.split('*')[-1]
    FNs = ['/'.join(FN[:-len(ext)].split('/')[-3:]) for FN in glob.glob(path)]
    FNs = [FN if os.path.basename(FN).find('-')==-1 else FN[:-2] for FN in FNs]
    return set(FNs)

  set1 = set()
  set2 = set()
  for p in path1:
    set1 = set1.union(complex_ids(p))
  for p in path2:
    set2 = set2.union(complex_ids(p))
  first_not_second = sorted([x for x in set1 if x not in set2])
  second_not_first = sorted([x for x in set2 if x not in set1])
  return (first_not_second, second_not_first)

def loopISM(ism_FNs, search_strings, field_width, nonzero=True):
  outstr = ''
  for ism_FN in ism_FNs:
    prefix = os.path.basename(ism_FN)[:-4]
    n = np.sum(np.array([numfiles(search_string%prefix, nonzero=nonzero) \
      for search_string in search_strings]))
    outstr += '%s\t%d\n'%(prefix.rjust(field_width), n)
  n = np.sum(np.array([numfiles(search_string%'', nonzero=nonzero) \
      for search_string in search_strings]))
  outstr += '%s\t%d\n'%('Total'.rjust(field_width), n)
  return outstr

start_dir = os.getcwd()
if args.multi_target:
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
    ism_nligands_dock = numfiles('ligand/dock_in/%s*/*.mol2'%prefix)
    ism_nligands_AlGDock = numfiles('ligand/AlGDock_in/%s*/*.tar.gz'%prefix)
    nligands_all += ism_nligands_all
    nligands_dock += ism_nligands_dock
    nligands_AlGDock += ism_nligands_AlGDock
    print '%*s\t%d\t%d\t%d'%(first_field_width,\
      prefix,ism_nligands_all,ism_nligands_dock,ism_nligands_AlGDock)
  print '%*s\t%d\t%d\t%d'%(first_field_width,\
    'Total',nligands_all,nligands_dock,nligands_AlGDock)

  print "\nReceptors"
  print "\tHomology  ", numfiles('receptor/3-models/pdb_noH/*.pdb')
  print "\tDock      ", numfiles('receptor/dock_in/*.sph')
  nonoverlap = FN_in_first_not_second('receptor/3-models/pdb_noH/*.pdb', \
                                   'receptor/dock_in/*.sph')
  if len(nonoverlap)>0:
    print '  in homology but not dock: ' + ', '.join(nonoverlap)
  print "\tAMBER     ", numfiles('receptor/amber_in/*.prmtop')
  nonoverlap = FN_in_first_not_second('receptor/3-models/pdb_noH/*.pdb', \
                                   'receptor/amber_in/*.prmtop')
  if len(nonoverlap)>0:
    print '  in homology but not AMBER: ' + ', '.join(nonoverlap)
  print "\tAlGDock   ", numfiles('receptor/AlGDock_in/*.PB.nc')
  nonoverlap = FN_in_first_not_second('receptor/3-models/pdb_noH/*.pdb', \
                                   'receptor/AlGDock_in/*.PB.nc')
  if len(nonoverlap)>0:
    print '  in homology but not AlGDock: ' + ', '.join(nonoverlap)

  print "\nComplexes for AlGDock"
  print loopISM(ism_FNs, ['complex/AlGDock_in/%s*/*/*.tar.gz'], \
    first_field_width)

  print "\nCompleted Calculations"
  print "for dock (*.mol2.gz)"
  print loopISM(ism_FNs, ['dock6/%s*/*/*.mol2.gz'], first_field_width, \
    nonzero=False)
  print "for dock (*.nc)"
  print loopISM(ism_FNs, ['dock6/%s*/*/*.nc'], first_field_width)
  print "for dock (all)"
  print loopISM(ism_FNs, ['dock6/%s*/*/*.mol2.gz', 'dock6/%s*/*/*.nc'], first_field_width, nonzero=False)

  print "for AlGDock"
  print loopISM(ism_FNs, ['AlGDock/dock/%s*/*/*/f_RL.pkl.gz'], \
    first_field_width)

  if args.clean_dock_not_complex:
    (complex_not_dock,dock_not_complex) = complex_nonoverlap(['complex/AlGDock_in/*/*/*.tar.gz'], ['dock6/*/*/*.mol2.gz', 'dock6/*/*/*.nc'])
    for p in dock_not_complex:
      os.system('rm dock6/%s*'%p)

  # print "\nComparing AlGDock and dock sets"
#  (AlGDock_not_nc,nc_not_AlGDock) = complex_nonoverlap(['AlGDock/dock/*/*/*/f_RL.pkl.gz'], ['dock6/*/*/*.nc'])
#  (AlGDock_not_dock,dock_not_AlGDock) = complex_nonoverlap(['AlGDock/dock/*/*/*/f_RL.pkl.gz'], ['dock6/*/*/*.mol2.gz', 'dock6/*/*/*.nc'])

#  import gzip, pickle
#  for path in AlGDock_not_nc:
#    F = gzip.open('AlGDock/dock/%s-0/f_RL.pkl.gz'%path)
#    dat = pickle.load(F)
#    F.close()
#    print dat[-1]['NAMD_OBC_MBAR'][-1]
#
#  for path in AlGDock_not_nc:
#    os.system('rm -rf AlGDock/dock/%s-0'%path)

  os.chdir(start_dir)
