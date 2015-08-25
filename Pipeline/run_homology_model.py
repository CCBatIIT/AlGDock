# Builds homology models

import os, inspect
dirs = {}
dirs['script'] = os.path.dirname(os.path.abspath(\
  inspect.getfile(inspect.currentframe())))
execfile(os.path.join(dirs['script'],'_external_paths.py'))
command_paths = findPaths(['qsub_command'])

# Parse arguments
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--sequence_ali', default='../1-search/seq.ali')
parser.add_argument('--template', \
  default='../2-binding_site/receptor_trans',\
  help='Directory with template PDB files (or single file)')
parser.add_argument('--reference', action='store_true', \
  help='Only builds a homology model based on the reference model')
parser.add_argument('--max_jobs', default=None, type=int)
parser.add_argument('--dry', action='store_true', default=False, \
  help='Does not actually submit the job to the queue')
parser.add_argument('--pylab', action='store_true')
args = parser.parse_args()

if args.reference:
  if os.path.isfile('../search_options.py'):
    execfile('../search_options.py')
    ref_pdb_chain_id = ref_pdb_id + ref_chain_id
    print 'Reference pdb is %s, chain %s'%(ref_pdb_id, ref_chain_id)
  else:
    raise Exception('Reference model not selected!')

if os.path.isfile(args.template):
  templateFNs = [args.template]
elif os.path.isdir(args.template):
  import glob
  templateFNs = glob.glob(os.path.join(args.template,'*.pdb'))
else:
  raise Exception('Invalid template!')

job_count = 0
for templateFN in templateFNs:
  basename = os.path.basename(templateFN)
  pdb_chain_id = basename[:5]
  if args.reference and pdb_chain_id!=ref_pdb_chain_id:
    continue

  out_FN = os.path.join('pdb_noH',pdb_chain_id+'.0.pdb')
  if os.path.isfile(out_FN):
    continue

  if not os.path.isdir(pdb_chain_id):
    os.makedirs(pdb_chain_id)
  os.chdir(pdb_chain_id)

  print 'Working in %s, building homology model(s) based on %s'%(\
    os.getcwd(),pdb_chain_id)
  jobname = pdb_chain_id
  command = '{0} {1}/homology_model.modeller.py ' + \
            '--sequence_ali ../{2} --template_pdb ../{3}'
  command = command.format('mod9.15', dirs['script'], \
                           args.sequence_ali, templateFN)
  print command

  print 'Submitting: ' + command
  import subprocess
  subprocess.call(['python',command_paths['qsub_command'],\
    jobname, command] + {True:['--dry'],False:[]}[args.dry])

  os.chdir('..')

  job_count += 1
  if (args.max_jobs is not None) and (job_count>=args.max_jobs):
    break
