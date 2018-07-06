#!/usr/bin/python

# Submits a command to the queue
import argparse
parser = argparse.ArgumentParser(description='Run a command on the queue')
parser.add_argument('name', help='Job name')
parser.add_argument('command', help='Command to execute')
parser.add_argument('--mem', type=int, default=2, help='Amount of memory to allocate (GB)')
parser.add_argument('--comment', default='', help='Adds a comment to the end of the script')
parser.add_argument('--dry', action='store_true', default=False, \
  help='Does not actually submit the job to the queue')
parser.add_argument('--no_release', action='store_true', default=False, \
  help='Does not release held jobs')

# For OSG
parser.add_argument('--disk', default='350MB', help='Amount of disk space to allocate')
parser.add_argument('--input_files', default=[], nargs='+', help='Input files for the job')
parser.add_argument('--output_files', default=[], nargs='+', help='Output files for the job')
parser.add_argument('--output_remaps', default=[], nargs='+', help='New output file names')
parser.add_argument('--min_job_time', default=20, type=int, \
  help='Will resubmit job if there is an error and the job takes less than this time (in minutes)')
parser.add_argument('--max_job_time', default=12, type=int, \
  help='Will hold job that takes longer than this time (in hours)')
# For CCB and DSCR
parser.add_argument('--nodes', type=int, default=1, help='Number of nodes to run the job on')
parser.add_argument('--ppn', type=int, default=1, help='Number of processors per node to run the job on')
parser.add_argument('--ambertools', action='store_true', default=False, \
  help='Load the ambertools/16 module')
parser.add_argument('--email', default='', help='Adds email to job')
parser.add_argument('--email_options', default='abe', help='Options for email notifications. When job begins (b), job ends (e), and/or aborted (a)')
args = parser.parse_args()

# Find unique name for the submission script
import os
curdir = os.getcwd()
n = 0
submit_FN = os.path.join(curdir,'jobs','%s-%d.submit'%(args.name,n))
while os.path.exists(submit_FN):
  n = n + 1
  submit_FN = os.path.join(curdir,'jobs','%s-%d.submit'%(args.name,n))
sh_FN = os.path.join(curdir,'jobs','%s-%d.sh'%(args.name,n))
out_FN = os.path.join(curdir,'jobs','%s-%d.out'%(args.name,n))
err_FN = os.path.join(curdir,'jobs','%s-%d.err'%(args.name,n))

# Looks for maximum job time, in minutes
max_time = 168*60
arg_list = args.command.split(' ')
try:
  max_time = int(arg_list[arg_list.index('--max_time')+1])
except ValueError:
  pass

# Sets up the submission and execution scripts
submit_script = ''
execute_script = ''
if os.path.exists('/share/apps/algdock'): # CCB Cluster
  cluster = 'CCB'
  
  # Split the command onto multiple lines
  command_list = args.command.split(';')
  command = '\n'.join([c.strip() for c in command_list])

  # By default, use miniconda
  if command.find('python')>-1:
    modules = 'module load miniconda/2\n'
    modules += 'module load gcc/6.2\n'
  else:
    modules = ''

  if command.find('chimera')>-1:
    modules += 'module load chimera/1.11\n'

  if command.find('modeller')>-1:
    modules += 'module load anaconda/2\n'

  if command.find('cores')>-1:
    cores = command[command.find('cores')+5:]
    cores = cores[:cores.find('\n')].strip()
    cores = cores.split(' ')[0]
    if cores!='':
      args.ppn = int(cores)

#  if args.ambertools:
#    modules += 'module load ambertools/16\n'
  
  email_specified = ''
  if args.email == '':
    email_specified = '#'

# ONLY USE IF THERE ARE ISSUES WITH SPECIFIC NODES!
# Selects certain specific nodes:
# nodeslist=['compute-1-1','compute-1-2','compute-1-6','compute-1-7','compute-1-8']
#  import numpy as np
#  args.nodes = nodeslist[np.random.randint(0,len(nodeslist))]

  max_time = '%02d:%02d:00'%(int(max_time/60), max_time%60)

  # Write script
  submit_script = '''#!/bin/bash
#
#PBS -S /bin/bash
#PBS -N {0}
#PBS -l mem={1}GB,nodes={2}:ppn={3},walltime={4}
#PBS -d {5}
#PBS -o {6}
#PBS -e {7}
#PBS -q default
{11}#PBS -M {12}
{11}#PBS -m {13}

hostname

{8}
{9}

# {10}
'''.format(args.name, args.mem, args.nodes, args.ppn, max_time, \
           curdir, out_FN, err_FN, \
           modules, command, args.comment, \
           email_specified, args.email, args.email_options)
elif os.path.exists('/pylon2') or os.path.exists('/oasis/projects/nsf'): # Bridges Cluster
  if os.path.exists('/pylon2'):
    cluster = 'Bridges'
  else:
    cluster = 'Comet'
  
  # Split the command onto multiple lines
  command_list = args.command.split(';')
  command = '\n'.join([c.strip() for c in command_list])

  if command.find('cores')>-1:
    cores = command[command.find('cores')+5:]
    cores = cores[:cores.find('\n')].strip()
    cores = cores.split(' ')[0]
    if cores!='':
      args.ppn = int(cores)

  email_specified = ''
  if args.email == '':
    email_specified = '#'

  max_time = min(max_time, 48*60)
  max_time = '%02d:%02d:00'%(int(max_time/60), max_time%60)

  # Write script
  submit_script = '''#!/bin/bash
#
#SBATCH --job-name={0}
#SBATCH --mem={1}
#SBATCH --nodes={2}
#SBATCH --ntasks-per-node={3}
#SBATCH -t {4}
#SBATCH --partition={5}
#SBATCH --workdir={6}
#SBATCH -o {7}
#SBATCH -e {8}

hostname

{9}

# {10}
'''.format(args.name, args.mem*1000, args.nodes, args.ppn, max_time, \
  {'Bridges':'RM-shared', 'Comet':'shared'}[cluster], \
  curdir, out_FN, err_FN, command, args.comment)
elif os.path.exists('/g/g19/minh1'): # LLNL cab
  cluster = 'cab'

  # Split the command onto multiple lines
  command_list = [c.strip() for c in args.command.strip().split(';')]
  command = '\n'.join([('srun -n 1 ' + c) for c in command_list if c!=''])

  if command.find('cores')>-1:
    cores = command[command.find('cores')+5:]
    cores = cores[:cores.find('\n')].strip()
    cores = cores.split(' ')[0]
    if cores!='':
      args.ppn = int(cores)

  max_time = min(max_time, 16*60)
  max_time = '%02d:%02d:00'%(int(max_time/60), max_time%60)

  submit_script = '''#!/bin/bash
#SBATCH -J {0}
#SBATCH -N {1}
#SBATCH -n {2}
#SBATCH -t {3}
#SBATCH -p pbatch
#SBATCH -D {4}
#SBATCH -o {5}
#SBATCH -e {6}

hostname

{7}

wait
'''.format(args.name, args.nodes, args.ppn, max_time, \
           curdir, out_FN, err_FN, command)
# weekday queue
#SBATCH -t 16:00:00
# for debugging:
#SBATCH -t 30:00
#SBATCH -p pdebug
elif os.path.exists('/stash') or os.path.exists('/home/dminh/perl5'):
  cluster = 'OSG' # OSG connect or OTSGrid
  
  # Split the command onto multiple lines
  command_list = args.command.split(';')
  command = '\n'.join([c.strip() for c in command_list])
  
  # Determine the input files
  #   All specified input files
  input_files = set([os.path.abspath(FN) for FN in args.input_files])
  #   The _external_paths.py script (6/13/2018 - this may not be necessary)
#  input_files = input_files.union(\
#    ['/home/daveminh/public/AlGDock-0.0.1/Pipeline/_external_paths.py'])
  #   Any files mentioned in the command
  command_list = command.split(' ')
  input_files = input_files.union(\
    set([os.path.abspath(FN) for FN in command_list if os.path.isfile(FN)]))
  #   Use the basename for any file, the current directory for any directory
  command = ' '.join([\
    {True:os.path.basename(item), False:
      {True:'.', False:item}[os.path.isdir(item)]}[os.path.isfile(item)] \
        for item in command_list])
  #   Any output file that already exists
  input_files = input_files.union(\
    set([os.path.abspath(FN) for FN in args.output_files if os.path.isfile(FN)]))
  transfer_input_files = ', '.join(input_files)

  # Format the output files for the script
  touches = ""
  if len(args.output_remaps)>0:
    transfer_output_files = 'transfer_output_remaps = "' + \
      '; '.join(['%s = %s'%(FNo,FNn) for (FNo,FNn) in zip(\
        args.output_remaps[::2],args.output_remaps[1::2])]) + '"'
    for FN in args.output_remaps[::2]:
      touches += "if [ ! -e {0} ]\n  then\n    touch {0}\nfi\n".format(FN)

  if command.find('$ALGDOCK')!=-1:
    hold_string = 'on_exit_hold = (ExitCode == 100)'
    requirements_string = '&& (HAS_CVMFS_oasis_opensciencegrid_org =?= TRUE)'
  else:
    hold_string = '''
# stay in queue if there was an error and 
# the job ran for less than min_job_time minutes
on_exit_hold = (ExitCode != 0) && ((CurrentTime - JobStartDate) < ({0}*60))
'''.format(args.min_job_time)
    requirements_string = ''

  # Write the submission script
  submit_script = """Universe       = vanilla

Executable     = {0}
Error   = jobs/{1}.$(Cluster)-$(Process).err
Output  = jobs/{1}.$(Cluster)-$(Process).out
Log     = jobs/{1}.$(Cluster).log
Requirements = (FileSystemDomain != "") && (OpSys == "LINUX" ) && (Arch == "X86_64") {2}
request_disk = {3}
request_memory = {4}GB

# File transfer
should_transfer_files = YES
transfer_input_files = {5}
{6}
when_to_transfer_output = ON_EXIT_OR_EVICT
want_graceful_removal = (ExitCode == 100)

{7}
# protect against hung jobs (taking more than max_job_time hours)
periodic_hold = (JobStatus==2) && ((CurrentTime - EnteredCurrentStatus) > {8}*60*60)
# make sure the job is being retried and rematched
{9}
+ProjectName="AlGDock"
Queue 1
""".format(sh_FN, args.name, requirements_string, args.disk, args.mem, \
    transfer_input_files, transfer_output_files, \
    hold_string, args.max_job_time, \
    {True:'',
    False:'periodic_release = ((CurrentTime - EnteredCurrentStatus) > 10 * 60) && (NumJobStarts < 40)'}[\
      args.no_release])

  if command.find('$ALGDOCK')!=-1:
    command = """

module load libgfortran

# Download data
wget --no-verbose --no-check-certificate http://stash.osgconnect.net/+daveminh/algdock.tar.gz
# wget --no-check-certificate 'https://docs.google.com/uc?export=download&id=0ByidOA_rkLLSbXl3WnZ3MmtWbnc' -O google_drive_downloader.tar.gz
# tar -xvf google_drive_downloader.tar.gz
# ./google_drive_downloader/google_drive_downloader
tar xzf algdock.tar.gz

# Modify paths
echo "
search_paths = {'MMTK':['$WORK_DIR/AlGDock/MMTK']}
" | cat AlGDock/AlGDock/_external_paths.py - > AlGDock/AlGDock/paths.py
mv AlGDock/AlGDock/paths.py AlGDock/AlGDock/_external_paths.py
export ALGDOCK=$WORK_DIR/AlGDock/BindingPMF

""" + command + """

rm -rf AlGDock namd* sander* ambpdb* molsurf* APBS*
rm -f algdock.tar.gz*
rm -f *.inpcrd* *.prmtop* *.frcmod *.pdb *.db
rm -f *.out *.namd *.dcd
rm -f .lock

"""
  execute_script = """#!/bin/bash

WORK_DIR=`pwd`
echo Working in $WORK_DIR
echo Directory before command:
ls -ltr

"""+command+touches+"""

echo Directory after command:
ls -ltr

"""

  if command.find('$ALGDOCK')!=-1 and command.find('timed')!=-1:
      # -s means file is not zero size
      dock_result = command[command.find('f_RL'):]
      dock_result = dock_result[:dock_result.find(' ')]
      execute_script += """
if [ ! -s {0} ]
  then
    exit 100
fi""".format(dock_result)
else:
  cluster = None
  commands = args.command.split(';')
  submit_script = '\n'.join([c.strip() for c in commands])

# Write and submit scripts
if submit_script!='':
  if not os.path.isdir('jobs'):
    os.makedirs('jobs')
  submit_F = open(submit_FN,'w')
  submit_F.write(submit_script)
  submit_F.close()

if execute_script!='':
  if not os.path.isdir('jobs'):
    os.makedirs('jobs')
  sh_F = open(sh_FN,'w')
  sh_F.write(execute_script)
  sh_F.close()

if (not args.dry) and cluster in ['OSG','CCB','Bridges','Comet']:
  print 'Submitting job script: ' + submit_FN

print('Job name: ' + args.name)

if not args.dry:
  if cluster=='CCB':
    os.system('qsub %s'%submit_FN)
  elif cluster in ['Bridges','Comet','cab']:
    os.system('sbatch %s'%submit_FN)
  elif cluster=='OSG':
    os.system('condor_submit %s'%submit_FN)
  else:
    os.system(args.command)
