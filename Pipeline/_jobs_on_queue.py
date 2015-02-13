# A function that returns names of the jobs on the queue
def jobs_on_queue():
  onq = []
  import os
  import subprocess

  if os.path.exists('/home/daveminh/stash'):
    condor_q = subprocess.Popen(['condor_q','-long','daveminh'], \
      stdout=subprocess.PIPE).stdout.read().split('\n\n')
    condor_q = [dict([(line[:line.find('=')-1],line[line.find('=')+2:]) for line in c.split('\n')]) for c in condor_q if c!='']
    onq = [q['Cmd'][1:-1] for q in condor_q if q['JobStatus']!='3'] # 3 is removed
    onq = [os.path.basename(line) for line in onq] # Get rid of directory name
    onq = ['-'.join(line.split('-')[:-1]) for line in onq] # Get rid of job number
  elif os.path.exists('/home/dminh/scripts/qsub_command.py'): # CCB cluster
    qstat = subprocess.Popen(['qstat','-f'], \
      stdout=subprocess.PIPE).stdout.read().split('\n')
    onq = [line.strip().split()[-1] for line in qstat if line.find('Job_Name')>-1]
  return onq
