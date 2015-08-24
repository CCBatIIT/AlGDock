job_block = 50

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--job_block', default=job_block, type=int, \
  help='Number of dockings per job')
parser.add_argument('--dry', action='store_true', default=False, \
  help='Does not actually submit the job to the queue')  
args = parser.parse_args()

# Find dock6_to_nc.py
import os, inspect
dirs = {}
dirs['script'] = os.path.dirname(os.path.abspath(\
  inspect.getfile(inspect.currentframe())))
execfile(os.path.join(dirs['script'],'_external_paths.py'))
dock6_to_nc_script = os.path.join(dirs['script'], 'dock6_to_nc.py')
command_paths = findPaths(['qsub_command'])

import glob

# Convert from mol2 to netcdf files
command_list = []
FNs = [FN for FN in glob.glob('dock6/*/*/*.mol2.gz') \
  if os.path.getsize(FN)>0]
FNs_c = []
outFNs_c = []
for FN in FNs:
  if not os.path.isfile(FN[:-8]+'.nc'):
    command_list.append('python {0} {1}'.format(dock6_to_nc_script, FN))
    FNs_c.append(FN)
    outFNs_c.append(FN[:-8]+'.nc')
  ncommands = len(command_list)
  if ncommands==args.job_block or ((ncommands>0) and (FN==FNs[-1])):
    command = '; '.join(command_list)
    print command
    os.system(' '.join(['python',command_paths['qsub_command'],\
        'dock6_to_nc', "'"+command+"'", \
        '--input_files', dock6_to_nc_script, ' '.join(FNs_c) + \
        '--output_files', ' '.join(outFNs_c), \
        {True:'--dry',False:''}[args.dry]]))
    command_list = []
    FNs_c = []
    outFNs_c = []
