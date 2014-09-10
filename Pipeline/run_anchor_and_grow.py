# Runs a set of anchor and grow calculation with UCSF DOCK 6
# To be run from the [TARGET]/dock6/ directory

# Constants
job_block = 5
# If it takes ~30 minutes per docking, then 5 docking should take 2.5 hours

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--ligand', \
  default='../ligand/dock_in', \
  help='The directory to look for ligand mol2 files or a single file (SYBYL atom types)')
parser.add_argument('--receptor', \
  default='../receptor/dock_in', \
  help='The directory to look for receptors or a single receptor site file (spheres). Grids should have the same prefix.')
parser.add_argument('--job_block', default=job_block, type=int, \
  help='Number of dockings per job')
parser.add_argument('--max_jobs', default=None, type=int)
parser.add_argument('--dry', action='store_true', default=False, \
  help='Does not actually submit the job to the queue')
args = parser.parse_args()

# Check for the existence of input files
import os, glob

if os.path.isfile(args.ligand):
  ligand_FNs = [args.ligand]
elif os.path.isdir(args.ligand):
  ligand_FNs = glob.glob(os.path.join(args.ligand,'*.mol2'))
else:
  raise Exception('Ligand input %s is not a file or directory!'%args.ligand)
ligand_FNs = [FN for FN in ligand_FNs if os.path.getsize(FN)>0]

if os.path.isfile(args.receptor):
  receptor_FNs = [args.receptor]
elif os.path.isdir(args.receptor):
  receptor_FNs = glob.glob(os.path.join(args.receptor,'*.sph'))
else:
  raise Exception('Receptor input %s is not a file or directory!'%args.receptor)
# Require nrg and bmp as well as sph files
receptor_FNs = [FN for FN in receptor_FNs \
  if os.path.isfile(FN[:-4]+'.nrg') and os.path.isfile(FN[:-4]+'.bmp')]

# Find anchor_and_grow.py
import os, inspect
dirs = {}
dirs['script'] = os.path.dirname(os.path.abspath(\
  inspect.getfile(inspect.currentframe())))
execfile(os.path.join(dirs['script'],'_external_paths.py'))
ancg_script = os.path.join(dirs['script'], 'anchor_and_grow.py')
command_paths = findPaths(['qsub_command'])

print 'Docking %d ligands to %d receptors'%(len(ligand_FNs),len(receptor_FNs))

import subprocess

command_list = []
out_FNs = []
out_remaps = []
code_list = []
job_count = 0
for receptor_FN in receptor_FNs:
  labels = {'receptor':os.path.basename(receptor_FN[:-4])}
  for ligand_FN in ligand_FNs:
    labels['ligand'] = os.path.basename(ligand_FN[:-5])
    labels['complex'] = labels['ligand']+'-'+labels['receptor']
    out_FN = 'ancg-'+labels['complex']+'.mol2.gz'
    if not os.path.exists(out_FN):
      command = 'python {0} {1} {2}'.format(
        os.path.join(dirs['script'],'anchor_and_grow.py'), ligand_FN, receptor_FN)
      command_list.append(command)
      out_FNs.append(out_FN)
      if os.path.basename(out_FN)!=out_FN:
        out_remaps.append(os.path.basename(out_FN))
        out_remaps.append(out_FN)
      code_list.append(labels['ligand'])
    ncommands = len(command_list)
    if ncommands==args.job_block or \
      ((ncommands>0) and (ligand_FN==ligand_FNs[-1])):
      command = '; '.join(command_list)
      name = labels['receptor'] + '-' + '.'.join(code_list)
      print command
      os.system(' '.join(['python',command_paths['qsub_command'],\
        name, "'"+command+"'", \
        '--input_files', ancg_script, \
        receptor_FN[:-4]+'.nrg', receptor_FN[:-4]+'.bmp', \
        '--output_files', ' '.join(out_FNs), \
        {True:'--output_remaps ' + ' '.join(out_remaps), \
         False:''}[len(out_remaps)>0], \
        {True:'--dry',False:''}[args.dry]]))
      command_list = []
      out_FNs = []
      out_remaps = []
      code_list = []
      job_count += 1
      print 'Submitted %d jobs'%(job_count)
      if (args.max_jobs is not None) and (job_count>=args.max_jobs):
        break
  if (args.max_jobs is not None) and (job_count>=args.max_jobs):
    break