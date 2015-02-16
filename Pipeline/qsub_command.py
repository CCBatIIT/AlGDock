#!/usr/bin/python

# Submits a command to the queue
import argparse
parser = argparse.ArgumentParser(description='Run a command on the queue')
parser.add_argument('name', help='Job name')
parser.add_argument('command', help='Command to execute')
parser.add_argument('--mem', type=int, default=1, help='Amount of memory to allocate (GB)')
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
  help='Load the ambertools/14 module')
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

# Sets up the submission and execution scripts
submit_script = ''
execute_script = ''
if os.path.exists('/home/dminh/scripts/qsub_command.py'): # CCB Cluster
  cluster = 'CCB'
  
  # Split the command onto multiple lines
  command_list = args.command.split(';')
  command = '\n'.join([c.strip() for c in command_list])

  # By default, use Enthought Canopy python
  if command.find('python')>-1:
    modules = 'module load canopy/1.5.0\n'
  else:
    modules = ''

  if command.find('chimera')>-1:
    modules += 'module load chimera/1.8.1\n'

  if command.find('modeller')>-1:
    modules += 'module load modeller/9.14\n'

  if args.ambertools:
    modules += 'module load ambertools/14\n'

  # Write script
  submit_script = '''#!/bin/bash
#
#PBS -S /bin/bash
#PBS -N {0}
#PBS -l mem={1}GB,nodes={2}:ppn={3},walltime=168:00:00
#PBS -d {4}
#PBS -o {5}
#PBS -e {6}

{7}
{8}

# {9}
'''.format(args.name, args.mem, args.nodes, args.ppn, \
           curdir, out_FN, err_FN, \
           modules, command, args.comment)
elif os.path.exists('/home/dbchem/dm225/scripts/qsub_command.py'): # DSCR Cluster
  cluster = 'DSCR'
  if not (args.nodes==1):
    nodes = "\n#$ -pe threaded %s\n"%args.nodes
  else:
    nodes = ""

  submit_script = '''#!/bin/bash
#
#$ -S /bin/bash -cwd
#$ -N {0}
#$ -l mem_free={1}G{2}
#$ -o {3} -j y

{4}

# {5}
'''.format(args.name, args.mem, nodes,
           out_FN, args.command, args.comment)
elif os.path.exists('/stash'):   # Open Science Grid
  cluster = 'OSG'
  
  # Split the command onto multiple lines
  command_list = args.command.split(';')
  command = '\n'.join([c.strip() for c in command_list])
  
  # Determine the input files
  #   All specified input files
  input_files = set([os.path.abspath(FN) for FN in args.input_files])
  #   The _external_paths.py script
  input_files = input_files.union(\
    ['/home/daveminh/public/AlGDock-0.0.1/Pipeline/_external_paths.py'])
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
  else:
    hold_string = '''
# stay in queue if there was an error and 
# the job ran for less than min_job_time minutes
on_exit_hold = (ExitCode != 0) && ((CurrentTime - JobStartDate) < ({0}*60))
'''.format(args.min_job_time)

  # Write the submission script
  submit_script = """Universe       = vanilla

Executable     = {0}
Error   = jobs/{1}.$(Cluster)-$(Process).err
Output  = jobs/{1}.$(Cluster)-$(Process).out
Log     = jobs/{1}.$(Cluster).log
Requirements = (FileSystemDomain != "") && (OpSys == "LINUX" ) && (Arch == "X86_64") 
request_disk = {2}
request_memory = {3}GB

# File transfer
should_transfer_files = YES
transfer_input_files = {4}
{5}
when_to_transfer_output = ON_EXIT_OR_EVICT
want_graceful_removal = (ExitCode == 100)

{6}
# protect against hung jobs (taking more than max_job_time hours)
periodic_hold = (JobStatus==2) && ((CurrentTime - EnteredCurrentStatus) > {7}*60*60)
# make sure the job is being retried and rematched
{8}
+ProjectName="AlGDock"
Queue 1
""".format(sh_FN, args.name, args.disk, args.mem, \
    transfer_input_files, transfer_output_files, \
    hold_string, args.max_job_time, \
    {True:'',
    False:'periodic_release = ((CurrentTime - EnteredCurrentStatus) > 60)'}[\
      args.no_release])

  if command.find('$ALGDOCK')!=-1:
    command = """

# Download data
wget --no-verbose --no-check-certificate http://stash.osgconnect.net/+daveminh/algdock.tar.gz
tar xzf algdock.tar.gz

# Modify paths
echo "
search_paths = {
  'gaff.dat':[None],
  'namd':[None],
  'sander':[None],
  'apbs':[None],
  'ambpdb':[None],
  'molsurf':[None],
  'MMTK':['$WORK_DIR/AlGDock/MMTK'],
  'vmd':[None]}
" | cat AlGDock/AlGDock/_external_paths.py - > AlGDock/AlGDock/paths.py
mv AlGDock/AlGDock/paths.py AlGDock/AlGDock/_external_paths.py
export ALGDOCK=$WORK_DIR/AlGDock/BindingPMF

""" + command + """

rm -rf AlGDock namd2* sander* ambpdb* molsurf* APBS*
rm algdock.tar.gz*
rm *.inpcrd.gz *.prmtop.gz
rm *.out *.namd *.dcd
rm .lock

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
      execute_script += """
if [ ! -s f_RL.pkl.gz ]
  then
    exit 100
fi"""
else:
  cluster = None
  submit_script = args.command

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

if (not args.dry) and cluster in ['OSG','CCB','DSCR']:
  print 'Submitting job script: ' + submit_FN

print('Job name: ' + args.name)
print('Script contents: ' + submit_script)

if not args.dry:
  if cluster=='OSG':
    os.system('condor_submit %s'%submit_FN)
  elif cluster=='CCB' or cluster=='DSCR':
    os.system('qsub %s'%submit_FN)
  else:
    os.system(args.command)
