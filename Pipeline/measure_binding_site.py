# Loads and clusters potential ligands from aligned chains
# Measures a putative binding site from the largest cluster
# Translates the ligands, receptors, and complex configurations

import os, sys

# Ligands to include, specifying (pdb_id, chain_id, rename, res_id)
selection = []
# Ligands to exclude, specifying (pdb_id, chain_id, rename, res_id)
exclude = []
# Ligands to exclude, only specifying resname
exclude_resname = ['HOH']
# The minimum number of atoms in a ligand
minimum_natoms = 8
# Ligands will be clustered using hierarchical clustering with
# clustering_method = (method, number_of_clusters)
# or not clustered at all.
# If clustering is None, then a plot will be generated
# showing a number of possible clusterings.
# If is specified, then colors of the clusters will be shown.
clustering_method = None
# If site_R is None, the site radius will come from rounding up the
# distance from the site center to the closest ligand center of mass
site_R = None

if os.path.isfile('../xtal_ligand_selection.py'):
  execfile('../xtal_ligand_selection.py')

if len(selection)>0:
  exclude = []
  exclude_resname = []

# Parse arguments
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--source_directory', default='../1-search/chains_aligned/')
parser.add_argument('--dest_directory', default='ligand_aligned')
parser.add_argument('--pylab', action='store_true')
args = parser.parse_args()

# Loads records for
# receptors and
# ligands that are selected or not excluded
ligands = {} # Ligand records
receptors = {}

import glob
chainFNs = glob.glob(os.path.join(args.source_directory,'*'))
for chainFN in chainFNs:
  basename = os.path.basename(chainFN)
  pdb_id = basename[:4]
  chain_id = basename[4]
  F = open(chainFN,'r')
  lines = F.read().split('\n')
  F.close()
  for line in lines:
    if line.startswith('ATOM  '):
      key = (pdb_id,chain_id)
      if not key in receptors.keys():
        receptors[key] = {'pdb':[]}
      receptors[key]['pdb'].append(line)
    if line.startswith('HETATM'):
      resname = line[17:20].strip()
      res_id = int(line[23:26])
      key = (pdb_id,chain_id,resname,res_id)
      if (len(selection)>0) and (key not in selection):
        continue
      if resname in exclude_resname:
        continue
      if key in exclude:
        continue
      if not key in ligands.keys():
        ligands[key] = {'pdb':[]}
      ligands[key]['pdb'].append(line)

# Extract receptor coordinates
import numpy as np

for key in receptors.keys():
  receptors[key]['crd'] = np.array([\
    (float(line[30:38]),float(line[38:46]),float(line[46:54])) \
      for line in receptors[key]['pdb']])

# Remove ligands which have fewer than minimum_natoms atoms
for key in ligands.keys():
  if len(ligands[key]['pdb'])<minimum_natoms:
    del ligands[key]

# Calculate the center of mass for each ligand
import inspect
dirs = {}
dirs['script'] = os.path.dirname(os.path.abspath(\
  inspect.getfile(inspect.currentframe())))
execfile(os.path.join(dirs['script'],'_masses.py'))
execfile(os.path.join(dirs['script'],'_external_paths.py'))
command_paths = findPaths(['dock6'])
dirs['dock6'] = os.path.abspath(os.path.dirname(command_paths['dock6']))

for key in ligands.keys():
  ligands[key]['crd'] = np.array([\
    (float(line[30:38]),float(line[38:46]),float(line[46:54])) \
      for line in ligands[key]['pdb']])
  ligMasses = [masses[line[76:78].strip().capitalize()] \
    for line in ligands[key]['pdb']]
  ligands[key]['com'] = \
    [sum(ligands[key]['crd'][:,c]*ligMasses)/sum(ligMasses) for c in range(3)]
  ligands[key]['dmax_COM2atom'] = \
    np.max(np.sqrt(np.sum((ligands[key]['crd']-ligands[key]['com'])**2,1)))
com = np.array([ligands[key]['com'] for key in ligands.keys()])

# Cluster ligand centers of mass
import scipy.cluster.hierarchy as hierarchy

z = hierarchy.linkage(com)
assignments = {}

clr = ['#2200CC' ,'#D9007E' ,'#FF6600' ,'#FFCC00' ,'#ACE600' ,'#0099CC' ,
       '#8900CC' ,'#FF0000' ,'#FF9900' ,'#FFFF00' ,'#00CC01' ,'#0055CC']

if not os.path.isdir('figures'):
  os.makedirs('figures')

import matplotlib.pyplot as plt
plt.clf()

if clustering_method is None:
  # Visualize possible clusters
  fig, axes23 = plt.subplots(2, 3)
  for method, axes in zip(['single', 'complete'], axes23):
    z = hierarchy.linkage(com, method=method)

    axes[0].plot(range(1, len(z)+1), z[::-1, 2])
    knee = np.diff(z[::-1, 2], 2)
    axes[0].plot(range(2, len(z)), knee)

    if len(knee)>0:
      # Good cluster sizes
      num_clust1 = knee.argmax() + 2
      knee[knee.argmax()] = 0
      num_clust2 = knee.argmax() + 2

      axes[0].text(num_clust1, z[::-1, 2][num_clust1-1], \
        'possible\n<- knee point')

      part1 = hierarchy.fcluster(z, num_clust1, 'maxclust')
      part2 = hierarchy.fcluster(z, num_clust2, 'maxclust')

      assignments[(method,num_clust1)] = part1
      assignments[(method,num_clust2)] = part2

      for part, ax in zip([part1, part2], axes[1:]):
        for cluster in set(part):
          if cluster<12: # TO DO: Use more colors
            ax.scatter(com[part==cluster,0], com[part==cluster,1],
                         color=clr[cluster])

      m = '\n(method: {})'.format(method)
      plt.setp(axes[0], title='Screeplot{}'.format(m), xlabel='partition',
               ylabel='{}\ncluster distance'.format(m))
      plt.setp(axes[1], title='{} Clusters'.format(num_clust1))
      plt.setp(axes[2], title='{} Clusters'.format(num_clust2))
  part = np.ones(com.shape[0],dtype=int)
else:
  z = hierarchy.linkage(com, method=clustering_method[0])
  part = hierarchy.fcluster(z, clustering_method[1], 'maxclust')
  fig = plt.figure()
  ax = fig.add_axes()
  for cluster in set(part):
    if cluster<12: # TO DO: Use more colors
      plt.scatter(com[part==cluster,0], com[part==cluster,1],
                     color=clr[cluster])
fig.savefig('figures/clusters.png')

# Keep ligands in the largest cluster
biggest_cluster = np.argmax(np.bincount(part))
ligand_keys = ligands.keys()
selected_ligand_keys = [ligand_keys[n] for n in range(len(ligands)) \
  if part[n]==biggest_cluster]
selected_com = np.array([ligands[key]['com'] for key in selected_ligand_keys])
del com
selected_receptor_keys = []
for key in selected_ligand_keys:
  if not (key[0],key[1]) in selected_receptor_keys:
    selected_receptor_keys.append((key[0],key[1]))

# Prepare output
def tee(val, F):
  sys.stdout.write(val)
  F.write(val)

logF = open('measured_binding_site.py','w')

tee('# Selected %d ligands and %d receptors\n'%(\
  len(selected_ligand_keys),len(selected_receptor_keys)), logF)

# Write separate pdb files for all the ligands
if not os.path.isdir('ligand_aligned'):
  os.mkdir('ligand_aligned')
else:
  os.system('rm ligand_aligned/*')
for key in selected_ligand_keys:
  FN = '{0[0]}{0[1]}_{0[2]}_{0[3]}.pdb'.format(key)
  F = open(os.path.join('ligand_aligned',FN),'w')
  F.write('\n'.join(ligands[key]['pdb']))
  F.close()

# Measure the range of the binding site
com_min = np.array([min(selected_com[:,c]) for c in range(3)])
com_max = np.array([max(selected_com[:,c]) for c in range(3)])
site_center_aligned = (com_min+com_max)/2.
dmax_COM2site_center = np.max(np.sqrt(\
  np.sum((selected_com-site_center_aligned)**2,1)))
if site_R is None:
  site_R = np.ceil(dmax_COM2site_center)

# Measure the size of the grid
# The grid will be larger than the site by at least
# the maximum distance from any ligand atom to its center of mass
dmax_COM2atom = np.max(np.array([ligands[key]['dmax_COM2atom'] \
  for key in selected_ligand_keys]))
half_edge_length = np.ceil(dmax_COM2site_center + dmax_COM2atom)
origin_aligned = site_center_aligned - half_edge_length

tee('\n# For translated systems,\n', logF)
tee('# Minimum center of mass: \ncom_min = ' + \
  repr(list(com_min-origin_aligned))+'\n', logF)
tee('# Maximum center of mass: \ncom_max = ' + \
  repr(list(com_max-origin_aligned))+'\n', logF)
tee('# Site center: ' + \
  repr(list(site_center_aligned-origin_aligned)) + '\n', logF)
tee('# Maximum ligand COM distance from site center: ' + \
  repr(dmax_COM2site_center) + '\n', logF)
tee('# Site radius\n', logF)
tee('site_R = ' + repr(site_R) + '\n',logF)

tee('# Maximum distance from ligand COM to any atom: ' + \
  repr(dmax_COM2atom) + '\n', logF)
tee('# in ligand '+'{0[0]}{0[1]}_{0[2]}_{0[3]}'.format(selected_ligand_keys[\
  np.argmax(np.array([ligands[key]['dmax_COM2atom'] \
    for key in selected_ligand_keys]))]) + '\n', logF)
tee('half_edge_length = ' + repr(half_edge_length) + '\n', logF)

# Translate coordinates
for key in selected_ligand_keys:
  ligands[key]['crd_trans'] = ligands[key]['crd']-origin_aligned
for key in selected_receptor_keys:
  receptors[key]['crd_trans'] = receptors[key]['crd']-origin_aligned

# Create directories for the translated files
for dest_directory in ['ligand_trans','receptor_trans','complex_trans']:
  if not os.path.isdir(dest_directory):
    os.mkdir(dest_directory)
  else:
    os.system('rm %s/*'%dest_directory)

# Write separate pdb files for all the ligands and complexes
for key in selected_ligand_keys:
  FN = '{0[0]}{0[1]}_{0[2]}_{0[3]}.pdb'.format(key)
  for (dest_directory, source_pdb, source_crds) in [\
      ('ligand_trans', ligands[key]['pdb'], ligands[key]['crd_trans']), \
      ('complex_trans', \
        receptors[(key[0],key[1])]['pdb'] + ligands[key]['pdb'], \
        np.vstack((receptors[(key[0],key[1])]['crd_trans'],
                   ligands[key]['crd_trans'])))]:
    pdb_trans = ['{0}{1[0]:8.3f}{1[1]:8.3f}{1[2]:8.3f}{2}'.format(\
      source_pdb[n][:30],\
      source_crds[n],\
      source_pdb[n][54:]) for n in range(len(source_pdb))]
    F = open(os.path.join(dest_directory,FN),'w')
    F.write('\n'.join(pdb_trans))
    F.close()

for key in selected_receptor_keys:
  FN = '{0[0]}{0[1]}.pdb'.format(key)
  dest_directory = 'receptor_trans'
  source_pdb = receptors[key]['pdb']
  source_crds = receptors[key]['crd_trans']
  pdb_trans = ['{0}{1[0]:8.3f}{1[1]:8.3f}{1[2]:8.3f}{2}'.format(\
    source_pdb[n][:30],\
    source_crds[n],\
    source_pdb[n][54:]) for n in range(len(source_pdb))]
  F = open(os.path.join(dest_directory,FN),'w')
  F.write('\n'.join(pdb_trans))
  F.close()

logF.close()

# Save the ligand centers of mass to a file
comF = open('ligand_trans/com.py','w')
comF.write('\n'.join(['shape sphere radius 0.25 center {0[0]},{0[1]},{0[2]}'.\
  format(com-origin_aligned) for com in selected_com]))
comF.close()

# Create a box showing the range of ligand centers of mass
print
showboxF = open('showbox.in','w')
showboxF.write('''N
U
{0} {0} {0}
{1[0]} {1[1]} {1[2]}
ligand_trans/com_box.pdb
'''.format(half_edge_length, (com_max-com_min)))
showboxF.close()
command = '{0}/showbox < showbox.in'.format(dirs['dock6'])
os.system(command)
os.remove('showbox.in')
