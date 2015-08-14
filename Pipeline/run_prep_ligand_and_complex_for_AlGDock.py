# Prepares a directory of ligand mol2 (sybyl atom types) and a
#          a directory of receptor pqr (amber atom types) files
# for AMBER and AlGDock
# Run from the [TARGET]/ directory

# Constants
job_block = 10

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--ligand_in', default='ligand/dock_in/', \
  help='The root directory to look for ligand mol2 files or a single file' + \
       ' (SYBYL atom types)')
parser.add_argument('--ligand_out', default='ligand/AlGDock_in/', \
  help='The root directory to put ligand prmtop, inpcrd, frcmod, and db files')
parser.add_argument('--receptor', default='receptor/3-models/pqr/', \
  help='The directory to look for receptors or a single receptor pqr file.')
parser.add_argument('--complex_out', default='complex/AlGDock_in/', \
  help='The root directory to put complex prmtop and inpcrd files')
parser.add_argument('--library_requirement', default=None, \
  help='The library name must contain the string LIBRARY_REQUIREMENT')
#  parser.add_argument('--complex_list', default='dock6/scores.txt',
#    help='A file that list labels of complexes to prepare' + \
#         '(e.g. with the top docking scores)')
#  parser.add_argument('--prepare_first', default=None, type=int,
#    help='The first PREPARE_FIRST molecules will be prepared for AlGDock' + \
#          '(if an argument is passed)')
parser.add_argument('--prepare_number', default=None, type=int,
  help='Up to PREPARE_NUMBER molecules will be prepared for AlGDock' + \
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

for path in ['ligand_in','ligand_out','receptor','complex_out']:
  setattr(args,path,os.path.abspath(getattr(args,path)))

import glob
if os.path.isfile(args.ligand_in):
  in_ligand_FNs = [args.ligand_in]
elif os.path.isdir(args.ligand_in):
  in_ligand_FNs = glob.glob(os.path.join(args.ligand_in,'*/*.mol2'))
  in_ligand_FNs = sorted([os.path.abspath(FN) for FN in in_ligand_FNs])
else:
  raise Exception('Ligand input %s is not a file or directory!'%args.ligand_in)
in_ligand_FNs = [os.path.abspath(FN) for FN in in_ligand_FNs if os.path.getsize(FN)>0]

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

print '%d ligand(s) and %d receptor(s) found'%(len(in_ligand_FNs),len(receptor_FNs))

if (args.library_requirement is not None):
  in_ligand_FNs = [FN for FN in in_ligand_FNs \
    if FN[len(args.ligand_in)+1:].find(args.library_requirement)>-1]
  print '%d ligand(s) meet the library requirement'%(len(in_ligand_FNs))

if (args.prepare_number is not None) and (args.prepare_number < len(in_ligand_FNs)):
  import numpy as np
  inds_o = set([int(np.floor(ind)) \
    for ind in np.linspace(0,len(in_ligand_FNs)-1,args.prepare_number)])
  in_ligand_FNs = [in_ligand_FNs[ind] for ind in inds_o]
  print '%d ligand(s) will be prepared'%(len(in_ligand_FNs))

# TODO: Prepare top ligands from dock6
#  if (args.prepare_first is not None) or (args.prepare_number is not None):
#    if os.path.isfile(args.complex_list):
#      F = open(args.complex_list,'r')
#      lines = F.read().split('\n')
#      if args.prepare_first is not None:
#        lines = lines[:args.prepare_first]
#      elif (args.prepare_number is not None) and (args.prepare_number < len(lines)):
#        import numpy as np
#        inds_o = [int(np.floor(ind)) for ind in np.linspace(0,len(lines)-1,args.prepare_number)]
#        inds = []
#        for ind in inds_o:
#          if not ind in inds:
#            inds.append(ind)
#        lines = [lines[ind] for ind in inds]
#      top_complexes = [line.split()[0] for line in lines]
#    else:
#      raise Exception('Docking scores not found!')
#  else:
#    top_complexes = None

# Submit a job for every args.job_block ligands
command_list = []
job_count = 0
for receptor_FN in receptor_FNs:
  labels = {'receptor':os.path.basename(receptor_FN[:-18])}
  for in_ligand_FN in in_ligand_FNs:
    labels['relpath_ligand'] = in_ligand_FN[len(args.ligand_in)+1:]
    labels['key'] = os.path.basename(in_ligand_FN[:-5])
    labels['library'] = '.'.join(os.path.dirname(labels['relpath_ligand']).split('.')[:-1])

    labels['complex'] = labels['library']+'.'+labels['key']+'-'+labels['receptor']

#    if (top_complexes is not None) and (labels['complex'] not in top_complexes):
#      continue

    command = ''
    out_ligand_FN = os.path.join(args.ligand_out,
      labels['library']+'.'+labels['key'][:-2]+'__', labels['key']+'.tar.gz')
    if not os.path.isfile(out_ligand_FN):
      command += 'python {0}/prep_ligand_for_AlGDock.py {1} {2}'
      command = command.format(dirs['script'], in_ligand_FN, out_ligand_FN, \
                               {True:'--debug', False:''}[args.debug])

    complex_FN = os.path.join(args.complex_out, \
      labels['library']+'.'+labels['key'][:-2]+'__', labels['key'], \
      labels['receptor']+'.tar.gz')
    if not os.path.exists(complex_FN):
      if command!='':
        command += '; '
      command += 'python {0}/prep_complex_for_AlGDock.py {1} {2} {3} {4}'
      command = command.format(
        dirs['script'], in_ligand_FN, receptor_FN, complex_FN, \
        {True:'--debug', False:''}[args.debug])

    if command!='':
      command_list.append(command)

    ncommands = len(command_list)
    if ncommands==args.job_block or \
      ((ncommands>0) and (in_ligand_FN==in_ligand_FNs[-1])):
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
