# Builds one homology model

# Parse arguments
#  import argparse
#  parser = argparse.ArgumentParser()
#  parser.add_argument('sequence_ali', default=None, help='Location of seq.ali')
#  parser.add_argument('template_pdb', default=None, help='Template')
#  parser.add_argument('--pylab', action='store_true')
#  args = parser.parse_args()

import optparse
parser = optparse.OptionParser()
parser.add_option('--sequence_ali', default='../1-search/seq.ali', \
  help='Location of seq.ali')
parser.add_option('--template_pdb', default=None, help='Template')
(args, options) = parser.parse_args()

import os
if not os.path.isfile(args.sequence_ali):
  raise Exception('Sequence missing or is not a file!')
if not os.path.isfile(args.template_pdb):
  raise Exception('Template PDB file missing or is not a file!')

print 'Generating homology model for sequence in %s based on %s'%(\
  args.sequence_ali, args.template_pdb)

basename = os.path.basename(args.template_pdb)
pdb_id = basename[:4]
chain_id = basename[4]

import time
start_time = time.time()

import modeller
import modeller.automodel

env = modeller.environ()

alignment_FN = '%s%s.ali'%(pdb_id,chain_id)
if not os.path.exists(alignment_FN):
  aln = modeller.alignment(env)
  mdl = modeller.model(env, file=args.template_pdb)
  # Model from PDB file
  aln.append_model(mdl, align_codes=pdb_id+chain_id, atom_files=args.template_pdb)
  # Reference sequence
  aln.append(file=args.sequence_ali, align_codes='all')
  aln.align2d()
  cdir = os.getcwd()
  aln.write(file=alignment_FN, alignment_format='PIR')
  print 'Completed 2D sequence alignment in %f s'%(time.time()-start_time)

# Read the aligned sequence of the crystal structure
F = open(alignment_FN,'r')
seqs = F.read().split('>P1;')
F.close()

knowns = seqs[1][:seqs[1].find('\n')]
seqid = seqs[2][:seqs[2].find('\n')]
sequence_pdb = ''.join(seqs[1].split('\n')[2:])

print 'The aligned sequence of the PDB file is ',sequence_pdb
print sequence_pdb.replace('-',' ').strip().find(' ')

# Do homology modelling
start_time = time.time()

loopExists = sequence_pdb.replace('-',' ').strip().find(' ')!=-1
if loopExists:
  print 'There are one or more loops in ', pdb_id, 'chain ', chain_id
  a = modeller.automodel.loopmodel(env,
    alnfile=alignment_FN,             # alignment filename
    knowns=knowns,                    # codes of the templates
    sequence=seqid,                   # code of the target
    loop_assess_methods=(modeller.automodel.assess.DOPE)) # assessment
  a.loop.starting_model = 0           # First loop model
  a.loop.ending_model   = 9           # Last loop model
  a.loop.md_level = modeller.automodel.refine.fast # Loop model refinement level
else:
  a = modeller.automodel.automodel(env,
    alnfile=alignment_FN,             # alignment filename
    knowns=knowns,                    # codes of the templates
    sequence=seqid)                   # code of the target

a.starting_model= 1                 # index of the first model
a.ending_model  = 1                 # index of the last model
# (determines how many models to calculate)
a.md_level = None

a.make()                            # do homology modeling
print 'Completed homology modeling in %f s'%(time.time()-start_time)

# Organize files
import glob, shutil

for dest_dir in ['pdb_noH','sequence_alignment']:
  if not os.path.isdir('../'+dest_dir):
    os.makedirs('../'+dest_dir)

pdb_FNs = []
if loopExists:
  # Rank models by DOPE score
  print 'Ranking models by DOPE score'
  ok_models = [x for x in a.loop.outputs if x['failure'] is None]
  key = 'DOPE score'
  ok_models.sort(key=lambda a: a[key])
  pdb_FNs = [m['name'] for m in ok_models]
else:
  pdb_FNs = glob.glob('*.B*.pdb')
for rep in range(len(pdb_FNs)):
  shutil.move(pdb_FNs[rep],'../pdb_noH/%s%s.%d.pdb'%(pdb_id,chain_id,rep))

FNs = glob.glob('*.ali')
for FN in FNs:
  shutil.move(FN, os.path.join('..','sequence_alignment',os.path.basename(FN)))

clearFNs = glob.glob(seqid+'*')
for FN in clearFNs:
  os.remove(FN)

os.chdir('..')
os.system('rm -rf %s%s'%(pdb_id,chain_id))