# Prepares a directory of ligand mol2 (sybyl atom types) and a
#          a directory of receptor pqr (amber atom types) files
# for AMBER and AlGDock
# Run from the [TARGET]/ directory

# TODO: Generated a complex_fixed_atoms pdb file

# Constants
job_block = 10

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--ligand_in', default='ligand/dock_in/', \
  help='The directory to look for ligand mol2 files or a single file' + \
       ' (SYBYL atom types)')
parser.add_argument('--ligand_out', default='ligand/AlGDock_in/', \
  help='The directory to put ligand prmtop, inpcrd, frcmod, and db files')
parser.add_argument('--receptor', default='receptor/3-models/pqr/', \
  help='The directory to look for receptors or a single receptor pqr file.')
parser.add_argument('--complex_out', default='complex/AlGDock_in/', \
  help='The directory to put complex prmtop and inpcrd files')
parser.add_argument('--complex_list', default='dock6/scores.txt',
  help='A file that list labels of complexes to prepare' + \
       '(e.g. with the top docking scores)')
parser.add_argument('--prepare_first', default=None, type=int,
  help='The first PREPARE_FIRST molecules will be prepared for AlGDock' + \
        '(if an argument is passed)')
parser.add_argument('--job_block', default=job_block, type=int, \
  help='Number of ligands per job')
parser.add_argument('--max_jobs', default=None, type=int)
parser.add_argument('--dry', action='store_true', default=False, \
  help='Does not actually submit the job to the queue')
parser.add_argument('--debug', action='store_true', default=False)
args = parser.parse_args()

import os, inspect
dirs = {'target':os.getcwd()}
dirs['script'] = os.path.dirname(os.path.abspath(\
  inspect.getfile(inspect.currentframe())))
execfile(os.path.join(dirs['script'],'_external_paths.py'))
prep_script = os.path.join(dirs['script'], '_prep_ligand.chimera.py')
command_paths = findPaths(['qsub_command'])

import glob
if os.path.isfile(args.ligand_in):
  ligand_FNs = [args.ligand_in]
elif os.path.isdir(args.ligand_in):
  ligand_FNs = glob.glob(os.path.join(args.ligand_in,'*.mol2'))
else:
  raise Exception('Ligand input %s is not a file or directory!'%args.ligand_in)
ligand_FNs = [os.path.abspath(FN) for FN in ligand_FNs if os.path.getsize(FN)>0]

if os.path.isfile(args.receptor):
  receptor_FNs = [args.receptor]
elif os.path.isdir(args.receptor):
  receptor_FNs = glob.glob(os.path.join(args.receptor,'*.pdb2pqr_amber.pqr'))
else:
  raise Exception('Receptor input %s is not a file or directory!'%args.receptor)
receptor_FNs = [os.path.abspath(FN) for FN in receptor_FNs]

for dirN in [args.ligand_out, args.complex_out]:
  if not os.path.isdir(dirN):
    os.system('mkdir -p '+dirN)

print '%d ligand(s) and %d receptor(s) found'%(len(ligand_FNs),len(receptor_FNs))

if args.prepare_first is not None:
  if os.path.isfile(args.complex_list):
    F = open(args.complex_list,'r')
    lines = F.read().split('\n')
    top_complexes = [line.split()[0] for line in lines[:args.prepare_first]]
  else:
    raise Exception('Docking scores not found!')
else:
  top_complexes = None

# Submit a job for every args.job_block ligands
import subprocess
import numpy as np

command_list = []
job_count = 0
for receptor_FN in receptor_FNs:
  labels = {'receptor':os.path.basename(receptor_FN[:-18])}
  for ligand_FN in ligand_FNs:
    labels['ligand'] = os.path.basename(ligand_FN[:-5])
    labels['complex'] = labels['ligand']+'-'+labels['receptor']
    if (top_complexes is not None) and (labels['complex'] not in top_complexes):
      continue

    command = ''
    if not np.array([\
      os.path.exists(os.path.join(args.ligand_out,labels['ligand']+'.'+key))
        for key in ['prmtop','inpcrd','frcmod','db']]).all():
      command += 'cd {0}; python {1}/prep_ligand_for_AlGDock.py {2} {3}; cd {4}'
      command = command.format(args.ligand_out, dirs['script'], \
                               ligand_FN, \
                               {True:'--debug', False:''}[args.debug], \
                               dirs['target'])

    if not np.array([\
      os.path.exists(os.path.join(args.complex_out,labels['complex']+'.'+key))
        for key in ['prmtop.gz','inpcrd.gz']]).all():
      if command!='':
        command += '; '
      command += 'cd {0}; python {1}/prep_complex_for_AlGDock.py {2} {3} {4}' + \
        ' --complex_prefix {5} {6}; gzip {5}.prmtop; gzip {5}.inpcrd; cd {7}'
      command = command.format(
        args.complex_out, dirs['script'], ligand_FN,
        os.path.join(dirs['target'],args.ligand_out,labels['ligand']+'.frcmod'), \
        receptor_FN, labels['complex'],
        {True:'--debug', False:''}[args.debug], dirs['target'])

    if command!='':
      command_list.append(command)

    ncommands = len(command_list)
    if ncommands==args.job_block or \
      ((ncommands>0) and (ligand_FN==ligand_FNs[-1])):
      command = '; '.join(command_list)
      print 'Submitting: ' + command
      os.system(' '.join(['python',command_paths['qsub_command'],\
        labels['complex'], "'"+command+"'", '--ambertools', \
        {True:'--dry',False:''}[args.dry]]))
      command_list = []
      job_count += 1
      print 'Submitted %d jobs\n'%job_count
      if (args.max_jobs is not None) and (job_count>=args.max_jobs):
        break
  if (args.max_jobs is not None) and (job_count>=args.max_jobs):
    break
