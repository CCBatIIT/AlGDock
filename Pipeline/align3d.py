# Downloads pdb files
# Uses ProDy to align a set of crystal structures to a reference structures.
# Also performs a principal components analysis.

import os, inspect, shutil, pickle

sequence = ''
ref_pdb_id = ''
ref_chain_id = 'A'
ref_res_id_range = (1,-1)
exclude = []

script_dir = os.path.dirname(os.path.abspath(\
  inspect.getfile(inspect.currentframe())))
exec(open(os.path.join(script_dir,'_load_profile.py')).read())

# Parse arguments
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--sequence', default=sequence)
parser.add_argument('--ref_pdb_id', default=ref_pdb_id.lower())
parser.add_argument('--ref_chain_id', default=ref_chain_id.upper())
parser.add_argument('--ref_res_id_range', default=ref_res_id_range)
parser.add_argument('--min_seq_identity', type=int, default=min_seq_identity)
parser.add_argument('--min_equivalent_positions', type=int, \
  default=min_equivalent_positions)
parser.add_argument('--pylab', action='store_true')
args = parser.parse_args()
del sequence, ref_pdb_id, ref_chain_id, min_seq_identity, min_equivalent_positions

if (args.sequence==''):
  # Read from seq.ali
  if os.path.isfile('seq.ali'):
    F = open('seq.ali','r')
    seq = F.read().split('\n')
    seq.pop(0)
    seq.pop(0)
    seq = ''.join(seq)
    args.sequence = seq[:-1]
  if (args.sequence==''):
    raise Exception('No sequence specified')

# Print arguments
print "align3d.py"
print "Sequence: ", args.sequence
print "Reference PDB Identifier:", args.ref_pdb_id
print "Reference Chain Identifier:", args.ref_chain_id
print "Reference Resid Range:", args.ref_res_id_range
print "Minimum sequence identity:", args.min_seq_identity
print "Minimum equivalent positions:", args.min_equivalent_positions

# Determine the PDB files and chains to download
chain_hits = []
pdb_hits = []
for prf in profile.items():
  if (prf[1][2]>=args.min_seq_identity) and \
     (prf[1][1]>=args.min_equivalent_positions):
    chain_hits.append(prf[0])
    if prf[0][0] not in pdb_hits:
      pdb_hits.append(prf[0][0])

print "\n*** Downloading original PDB files ***"
if not os.path.isdir('pdb_original'):
  os.makedirs('pdb_original')
for pdb_id in pdb_hits:
  if not os.path.isfile('pdb_original/%s.pdb'%pdb_id):
    print 'Downloading '+pdb_id
    os.system('curl --compressed http://www.rcsb.org/pdb/files/%s.pdb.gz > pdb_original/%s.pdb.gz'%(pdb_id,pdb_id))
    os.system('gunzip pdb_original/%s.pdb.gz'%pdb_id)

print "\n*** Separating chains ***"
refF = open('reference.pdb','w')
if not os.path.isdir('chains_original'):
  os.mkdir('chains_original')
for pdb_id in pdb_hits:
  inF = open('pdb_original/%s.pdb'%pdb_id,'r')
  lines = inF.read().split('\n')
  inF.close()
  outFs = {}
  res_id_history = {}
  adjusted_res_id_history = {}
  for line in lines:
    if line.startswith('ATOM  ') or line.startswith('HETATM'):
      res_name = line[17:20]
      chain_id = line[21]
      res_id = int(line[22:26])
      if ((pdb_id,chain_id) in chain_hits):
        # Open the pdb file in chains_original
        if not (chain_id in outFs):
          outFN = os.path.join('chains_original', pdb_id+chain_id+'.pdb')
          outFs[chain_id] = open(outFN,'w')
          res_id_history[chain_id] = []
          adjusted_res_id_history[chain_id] = []
          first_res_id = res_id # Subtract this to make the first residue 1
        # Renumber the residues
        if not res_id in res_id_history[chain_id]:
          res_id_history[chain_id].append(res_id)
          # If it is a new residue, find a res_id that
          #   aligns with the target sequence, and
          #   does not overlap with previous res_id
          if line.startswith('ATOM  '):
            adjusted_res_id = res_id - (first_res_id - 1) + \
              (args.ref_res_id_range[0] - 1) + \
              (profile[(pdb_id,chain_id)][0][0] - 1) - \
              (profile[(pdb_id,chain_id)][0][2] - 1)
          else:
            adjusted_res_id = res_id
          if (adjusted_res_id<1) or \
             (adjusted_res_id in adjusted_res_id_history[chain_id]):
            while (adjusted_res_id<1) or \
                  (adjusted_res_id in adjusted_res_id_history[chain_id]):
              adjusted_res_id += 1
            if not res_name.startswith('HOH'):
              print 'In %s, chain %s, residue %s resid %d reassigned to %d'%(\
              pdb_id, chain_id, res_name, res_id, adjusted_res_id)
          adjusted_res_id_history[chain_id].append(adjusted_res_id)
        else:
          # If it is an old residue, use the last adjusted_res_id value
          adjusted_res_id = adjusted_res_id_history[chain_id][-1]
        line = line[:22] + '%4d'%adjusted_res_id + line[26:]
        # Write atoms
        if line.startswith('ATOM'):
          if (adjusted_res_id>=(args.ref_res_id_range[0] + \
                                (profile[(pdb_id,chain_id)][0][2]-1))) and \
             (adjusted_res_id<=args.ref_res_id_range[1]):
            outFs[chain_id].write(line+'\n')
            if (pdb_id==args.ref_pdb_id) and (chain_id==args.ref_chain_id):
              refF.write(line+'\n')
        else:
          outFs[chain_id].write(line+'\n')

for chain_id in outFs.keys():
  outFs[chain_id].close()
if refF is not None:
  refF.close()

print "There are %d PDB entries and %d chains"%(len(pdb_hits),len(chain_hits))
print "with a minimum sequence identity of %d"%args.min_seq_identity
print "and at least %d equivalent positions "%args.min_equivalent_positions

print "\n*** Aligning chains to reference structure ***"
import prody

# Load reference structure
if os.path.isfile('reference.pdb'):
  reference_structure = prody.parsePDB('reference.pdb',
    subset='ca', chain=args.ref_chain_id)
  if reference_structure is None:
    raise Exception('Error loading reference structure')
else:
  raise Exception('Reference structure not found')
reference_hierview = reference_structure.getHierView()
reference_chain = reference_hierview[args.ref_chain_id]

# Generate an ensemble of structures from the pdb files
ens_FN = 'prody.ens.npz'
mappings_FN = 'prody.mappings.pkl'
if os.path.isfile(ens_FN) and os.path.isfile(mappings_FN):
  ensemble = prody.loadEnsemble(ens_FN)
  mappings = pickle.load(open(mappings_FN,'r'))
else:
  ensemble = prody.PDBEnsemble()
  # Set ensemble atoms
  ensemble.setAtoms(reference_chain)
  # Set reference coordinates
  ensemble.setCoords(reference_chain.getCoords())

  mappings = []
  new_pdb_FNs = []
  for (pdb_id,chain_id) in chain_hits:
    pdb_FN = 'chains_original/%s%s.pdb'%(pdb_id,chain_id)
    if pdb_id in exclude:
      continue
    structure = prody.parsePDB(pdb_FN, subset='ca', chain=chain_id)
    if structure is None:
      prody.plog('Failed to parse ' + pdb_FN)
      continue
    mapping = prody.mapOntoChain(structure, reference_chain, seqid=args.min_seq_identity)
    if len(mapping) == 0:
      prody.plog('Failed to map', pdb_id)
      continue
    mappings.append(mapping)
    atommap = mapping[0][0]
    ensemble.addCoordset(atommap, weights=atommap.getFlags('mapped'))
    new_pdb_FNs.append(pdb_FN)

  print '%d chains map onto the reference sequence\n'%(len(mappings))
  ensemble.iterpose()
  prody.saveEnsemble(ensemble, filename=ens_FN[:-8])
  pickle.dump(mappings, open(mappings_FN,'w'))

  output_dir = os.path.join('chains_aligned')
  if not os.path.exists(output_dir):
    os.mkdir(output_dir)
  for pdb_FN, conf in zip(new_pdb_FNs,ensemble):
    trans = conf.getTransformation()
    chain_coords = prody.parsePDB(pdb_FN)
    trans.apply(chain_coords)
    outFN = os.path.join(output_dir,os.path.basename(pdb_FN))
    prody.writePDB(outFN, chain_coords)

print 'There are %d structures in the ensemble\n'%len(ensemble)

print "\n*** Principal Components Analysis ***"
pca_FN = os.path.join('prody.pca.npz')
if os.path.exists(pca_FN):
  pca = prody.loadModel(pca_FN)
else:
  pca = prody.PCA()
  pca.buildCovariance(ensemble) # Build covariance matrix
  pca.calcModes() # Calculate modes
  prody.saveModel(pca, filename=pca_FN[:-8])

if not os.path.isdir('figures'):
  os.makedirs('figures')

import matplotlib.pyplot as plt

if not os.path.isfile('rmsd.png'):
  rmsd = prody.calcRMSD(ensemble)
  plt.clf()
  plt.plot(rmsd);
  plt.xlabel('Conformation index');
  plt.ylabel('RMSD (A)');
  plt.title('RMSD %f (%f)'%(rmsd.mean(), rmsd.std()))
  plt.savefig('figures/rmsd.png')

if not os.path.isfile('blastPCA.png'):
  pc_ind0 = 0
  pc_ind1 = 1
  xtal_projection = prody.calcProjection(ensemble, pca[:20], rmsd=False)

  plt.clf()
  plt.plot(xtal_projection[:,pc_ind0],xtal_projection[:,pc_ind1],'ks')
#  titles = ['%s%s'%(pdb_id,chain_id) for (pdb_id,chain_id) in chain_hits]
  titles = []
  for map in mappings:
    titles += [map[0][0]._title[13:17]+'%d'%r \
               for r in range(len(map[0][0].getCoordsets()))]

  for label, x, y in zip(titles, \
                         xtal_projection[:,pc_ind0], \
                         xtal_projection[:,pc_ind1]):
    plt.annotate(label,
                 xy = (x, y), xytext = (-20, 20),
                 textcoords = 'offset points', ha = 'right', va = 'bottom',
                 bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
                 arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
  plt.savefig('figures/blastPCA.png')
