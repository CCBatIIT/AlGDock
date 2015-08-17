# Reorganizes pipeline directory structure from before 8/14/2015

import os
import glob

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--prefix', default='*')
args = parser.parse_args()

if args.prefix=='*':
  # Find qsub_command
  import os, inspect
  dirs = {}
  dirs['script'] = os.path.dirname(os.path.abspath(\
    inspect.getfile(inspect.currentframe())))
  execfile(os.path.join(dirs['script'],'_external_paths.py'))
  command_paths = findPaths(['qsub_command'])

  import subprocess
  dirNs = [N for N in glob.glob('*') if os.path.isdir(N)]
  for dirN in dirNs:
    path = os.path.join(dirs['script'], 'reorganize_pipeline.py')
    command = 'python %s --prefix %s'%(path, dirN)
    print command
    subprocess.call(['python', command_paths['qsub_command'], dirN, command])
  import sys
  sys.exit()

cdir = os.path.abspath('.')

FNs = glob.glob('%s/ligand/*.ism'%args.prefix)
libraries = [FN.split('/')[-1][:-4] for FN in FNs]
print 'Libraries: ' + ' '.join(libraries)

# Rename dock6 ligand input
FNs = glob.glob('%s/ligand/dock_in/*.mol2'%args.prefix)

for FN in FNs:
  baseN = os.path.basename(FN)
  dirN = os.path.join(os.path.dirname(FN),baseN[:-7]+'__')
  if not os.path.isdir(dirN):
    os.system('mkdir -p '+dirN)
  os.rename(FN,os.path.join(dirN,baseN[-8:]))
  print 'mv %s %s'%(FN,os.path.join(dirN,baseN[-8:]))

# Rename dock6 output
FNs = glob.glob('%s/dock6/*.mol2.gz'%args.prefix)
for FN in FNs:
  baseN = os.path.basename(FN)
  if baseN.startswith('ancg-'):
    baseN = baseN[5:]
  baseN_split = baseN.split('-')
  ligand = '-'.join(baseN_split[:-1])
  receptor = '.'.join(baseN_split[-1].split('.')[:-2])
  dirN = os.path.join(os.path.dirname(FN), ligand[:-2]+'__', ligand.split('.')[-1])
  if not os.path.isdir(dirN):
    os.system('mkdir -p '+dirN)
  os.rename(FN, os.path.join(dirN, receptor+'.mol2.gz'))
  print 'mv %s %s'%(FN, os.path.join(dirN, receptor+'.mol2.gz'))

# Convert from mol2 to netcdf files
FNs = glob.glob('%s/dock6/*/*/*.mol2.gz'%args.prefix)
for FN in FNs:
  if not os.path.isfile(FN[:-8]+'.nc'):
    print 'Converting %s to nc'%FN
    os.system('python $ALGDOCKHOME/Pipeline/dock6_to_nc.py '+FN)
    # os.remove(FN)

# Rename db files to lowercase 
FNs = glob.glob('%s/ligand/AlGDock_in/*.db'%args.prefix)
for FN in FNs:
  baseN = os.path.basename(FN).lower()
  if baseN!=os.path.basename(FN):
    os.rename(FN, os.path.join(os.path.dirname(FN), baseN))
    print 'mv %s %s'%(FN, os.path.join(os.path.dirname(FN), baseN))

# Compress ligand files
FNs = [os.path.abspath(FN) for FN in glob.glob('%s/ligand/AlGDock_in/*.prmtop'%args.prefix)]
for FN in FNs:
  ligand = os.path.basename(FN)[:-7]
  dirN = os.path.join(os.path.dirname(FN), ligand[:-2]+'__')
  if not os.path.isdir(dirN):
    os.system('mkdir -p '+dirN) 
  os.chdir(os.path.dirname(FN))
  print 'Creating %s/%s.tar.gz'%(dirN, ligand.split('.')[-1])
  os.system('tar -cf %s/%s.tar.gz %s.* %s.*'%(dirN, ligand.split('.')[-1], ligand, ligand.lower()))
  os.system('rm %s.*'%ligand)
  os.system('rm %s.*'%ligand.lower())
os.chdir(cdir)

# Compress complex files
FNs = [os.path.abspath(FN) for FN in glob.glob('%s/complex/AlGDock_in/*.gz'%args.prefix)]
if len(FNs)>0:
  print 'Decompressing complex prmtop and inpcrd files'
  os.system('gunzip %s/complex/AlGDock_in/*.gz'%args.prefix)

FNs = [os.path.abspath(FN) for FN in glob.glob('%s/complex/AlGDock_in/*.prmtop'%args.prefix)]
for FN in FNs:
  baseN_split = os.path.basename(FN).split('-')
  ligand = '-'.join(baseN_split[:-1])
  receptor = '.'.join(baseN_split[-1].split('.')[:-1])
  dirN = os.path.abspath(os.path.join(os.path.dirname(FN), ligand[:-2]+'__', ligand.split('.')[-1]))
  if not os.path.isdir(dirN):
    os.system('mkdir -p '+dirN)
  os.chdir(os.path.dirname(FN))
  print 'Creating %s/%s.tar.gz'%(dirN,receptor)
  os.system('tar -cf %s/%s.tar.gz %s-%s.*'%(dirN,receptor,ligand,receptor))
  os.system('rm %s-%s.*'%(ligand,receptor))
os.chdir(cdir)

# Reorganize AlGDock cool directory
FNs = glob.glob('%s/AlGDock/cool/*/cool_progress.pkl.gz'%args.prefix)
for FN in FNs:
  dirN_o_split = os.path.dirname(FN).split('/')[-1].split('-')
  ligand = '-'.join(dirN_o_split[:-1])
  rep = dirN_o_split[-1]
  dirN = os.path.join('/'.join(os.path.dirname(FN).split('/')[:-1]), ligand[:-2]+'__', ligand.split('.')[-1]+'-'+rep)
  if not os.path.isdir(dirN):
    os.system('mkdir -p %s'%dirN)
  print 'Moving cooling directory %s to %s'%(os.path.dirname(FN),dirN)
  os.system('mv %s/* %s'%(os.path.dirname(FN),dirN))
  os.system('rm -rf %s'%os.path.dirname(FN))

# Reorganize AlGDock dock directory
FNs = glob.glob('%s/AlGDock/dock/*/dock_progress.pkl.gz'%args.prefix)
for FN in FNs:
  dirN_o_split = os.path.dirname(FN).split('/')[-1].split('-')
  ligand = '-'.join(dirN_o_split[:-2])
  receptor = dirN_o_split[-2]
  rep = dirN_o_split[-1]
  dirN = os.path.join('/'.join(os.path.dirname(FN).split('/')[:-1]), ligand[:-2]+'__', ligand.split('.')[-1], receptor+'-'+rep)
  if not os.path.isdir(dirN):
    os.system('mkdir -p %s'%dirN)
  print 'Moving cooling directory %s to %s'%(os.path.dirname(FN),dirN)
  os.system('mv %s/* %s'%(os.path.dirname(FN),dirN))
  os.system('rm -rf %s'%os.path.dirname(FN))
