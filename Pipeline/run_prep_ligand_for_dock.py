import os
import glob

# Constants
job_block = 10

key = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'
# Returns a 3-letter code corresponding to a base-10 number
def code(val):
    return key[(val/36/36)%36]+key[(val/36)%36]+key[val%36]
# Returns a base-10 number corresponding to a 3-letter code
def base10val(code):
    return key.find(code[0])*36*36 + key.find(code[1])*36 + key.find(code[2])

import inspect
dirs = {}
dirs['script'] = os.path.dirname(os.path.abspath(\
  inspect.getfile(inspect.currentframe())))
execfile(os.path.join(dirs['script'],'_external_paths.py'))
prep_script = os.path.join(dirs['script'], '_prep_ligand.chimera.py')
command_paths = findPaths(['qsub_command'])

import sys
if len(sys.argv)==1:
  sys.argv.append('DUD-E.active.ism')

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('library', default=None, \
  help='ISM file of the library')
parser.add_argument('--count', default=None, type=int, \
  help='Number of ligands to process')
parser.add_argument('--job_block', default=job_block, type=int, \
  help='Number of ligands per job')
parser.add_argument('--max_jobs', default=None, type=int)
parser.add_argument('--dry', action='store_true', default=False, \
  help='Does not actually submit the job to the queue')
args = parser.parse_args()

# Load the library
if not os.path.isfile(args.library):
  raise Exception('Library file missing!')
  
F = open(args.library,'r')
lines = F.read().strip().split('\n')
F.close()

# Choose ligands to prepare
nligands = len(lines)
if args.count is None:
  inds = range(nligands)
else:
  inds = range(0,nligands,int(float(nligands)/args.count))

print 'Preparing %d ligands'%nligands

# Submit a job for every args.job_block ligands
if not os.path.isdir('dock_in'):
  os.makedirs('dock_in')

import subprocess

command_list = []
out_FNs = []
out_remaps = []
code_list = []
job_count = 0
for ind in inds:
  code_i = code(ind)
  out_FN = 'dock_in/{0}.{1}.mol2'.format(args.library[:-4],code_i)
  if not os.path.exists(out_FN):
    command = 'python {0} {1} "{2}"'.format(\
      os.path.join(dirs['script'],'prep_ligand_for_dock.py'), \
      out_FN, lines[ind].split()[0])
    command_list.append(command)
    out_FNs.append(out_FN)
    out_remaps.append(os.path.basename(out_FN))
    out_remaps.append(out_FN)
    code_list.append(code_i)
  ncommands = len(command_list)
  if ncommands==args.job_block or ((ncommands>0) and (ind==inds[-1])):
    command = '; '.join(command_list)
    name = args.library[:-4] + '.' + '.'.join(code_list)
    print 'Submitting: ' + command
    os.system(' '.join(['python',command_paths['qsub_command'],\
      name, "'"+command+"'",\
      '--input_files', prep_script,\
      '--output_files', ' '.join(out_FNs), \
      '--output_remaps', ' '.join(out_remaps), \
      {True:'--dry',False:''}[args.dry]]))
    command_list = []
    out_FNs = []
    out_remaps = []
    code_list = []
    job_count += 1
    print 'Submitted %d jobs\n'%job_count
    if (args.max_jobs is not None) and (job_count>args.max_jobs):
      break
