# Calculates alchemical grids
# Run from: receptor/AlGDock_in

import os, inspect
dirs = {}
dirs['script'] = os.path.dirname(os.path.abspath(\
  inspect.getfile(inspect.currentframe())))
execfile(os.path.join(dirs['script'],'_external_paths.py'))
command_paths = findPaths(['qsub_command'])

# Parse arguments
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--prmtop', default='../amber_in',\
  help='Directory with AMBER prmtop and inpcrd files (or single pair)')
parser.add_argument('--spacing', nargs=3, type=float, \
  default = [0.25, 0.25, 0.25],
  help='Grid spacing (overrides header)')
parser.add_argument('--counts', nargs=3, type=int, help='Number of point in each direction (overrides header)')
parser.add_argument('--max_jobs', default=None, type=int)
parser.add_argument('--dry', action='store_true', default=False, \
  help='Does not actually submit the job to the queue')
parser.add_argument('--pylab', action='store_true')
args = parser.parse_args()

if os.path.isfile(args.prmtop):
  prmtopFNs = [args.prmtop]
elif os.path.isdir(args.prmtop):
  import glob
  prmtopFNs = glob.glob(os.path.join(args.prmtop,'*.prmtop'))
else:
  raise Exception('Invalid prmtop file/directory!')
prmtopFNs = [FN for FN in prmtopFNs if os.path.isfile(FN[:-7]+'.inpcrd')]

print '%d prmtop/inpcrd pairs found'%len(prmtopFNs)

job_count = 0
for prmtop_FN in prmtopFNs:
  inpcrd_FN = prmtop_FN[:-7]+'.inpcrd'
  prefix = os.path.basename(prmtop_FN)[:-7]
  out_FNs = {'PB':prefix+'.PB.nc'}
  for key in ['ele','LJa','LJr']:
    # name includes prefix, key, and spacing (in pm)
    out_FNs[key] = '%s.%s.%d.nc'%(prefix,key,int(args.spacing[0]*100))
  import numpy as np
  if np.array([os.path.isfile(out_FNs[key]) for key in out_FNs.keys()]).all():
    continue # Grids already calculated
  
  print 'Calculating interaction grids for %s, with spacing of %f A'%(\
    prefix, args.spacing[0])
  jobname = '%s.%d'%(prefix, args.spacing[0]*100)
  command = 'python {0}/alchemicalGrids.py --prmtop_FN {1} --inpcrd_FN {2}' + \
    ' --pqr_FN {3}.pqr --PB_FN {3}.PB.nc --ele_FN {3}.ele.{4}.nc' + \
    ' --LJa_FN {3}.LJa.{4}.nc --LJr_FN {3}.LJr.{4}.nc' + \
    ' --spacing {5[0]} {5[1]} {5[2]}' + \
    {True:'', False:' --counts {6[0]} {6[1]} {6[2]}'}[args.counts is None]
  command = command.format(dirs['script'], prmtop_FN, inpcrd_FN, prefix, \
    int(args.spacing[0]*100), args.spacing, args.counts)
  print command

  print 'Submitting: ' + command
  import subprocess
  subprocess.call(['python',command_paths['qsub_command'],\
    jobname, command, '--ambertools'] + {True:['--dry'],False:[]}[args.dry])

  job_count += 1
  if (args.max_jobs is not None) and (job_count>=args.max_jobs):
    break
