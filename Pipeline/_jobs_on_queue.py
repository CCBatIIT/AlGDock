# A function that returns names of the jobs on the queue
def jobs_on_queue():
  onq = []
  import os
  import subprocess

  if os.path.exists('/home/daveminh/stash'):
    qstat = subprocess.Popen(['condor_q','-long','daveminh'], \
      stdout=subprocess.PIPE).stdout.read().split('\n')
    onq = [os.path.basename(l.split('"')[-2]) \
      for l in qstat if l.startswith('Cmd')]
    onq = ['-'.join(line.split('-')[:-1]) for line in onq] # Get rid of job number
  elif os.path.exists('/home/dminh/scripts/qsub_command.py'): # CCB cluster
    qstat = subprocess.Popen(['qstat','-f'], \
      stdout=subprocess.PIPE).stdout.read().split('\n')
    onq = [line.strip().split()[-1] for line in qstat if line.find('Job_Name')>-1]
  return onq

#  elif exists('/home/dbchem/dm225/scripts/qsub_command.py'):
#    cluster = 'DSCR'
#    qstat = subprocess.Popen(['qstat','-xml'],stdout=subprocess.PIPE).stdout.read().split('\n')
#    onq = [line.strip()[9:-14] for line in qstat if line.strip().startswith('<JB_name>')]
#    onq = ['-'.join(line.split('-')[:-1]) for line in onq] # Get rid of job number
