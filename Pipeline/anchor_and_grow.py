# Runs an anchor and grow calculation with UCSF DOCK 6
# To be run from the [TARGET]/dock/ directory

try:
  import argparse
  parser = argparse.ArgumentParser()
  parser.add_argument('ligand', \
    default='../ligand/dock_in/DUD-E.active.AAA.mol2', \
    help='Ligand mol2 file (SYBYL atom types)')
  parser.add_argument('receptor', \
    default='../receptor/dock_in/3cqwA.0.sph', \
    help='Receptor site file (spheres). Grids should have the same prefix.')
  args = parser.parse_args()
except:
  import sys
  class args:
    ligand = sys.argv[1]
    receptor = sys.argv[2]

# Check for the existence of input files
import os
for FN in [args.ligand, args.receptor, \
    args.receptor[:-4]+'.bmp', args.receptor[:-4]+'.nrg']:
  if not os.path.isfile(FN):
    raise Exception(FN + ' is missing!')

# Find dock6
import os, inspect
dirs = {}
dirs['script'] = os.path.dirname(os.path.abspath(\
  inspect.getfile(inspect.currentframe())))
execfile(os.path.join(dirs['script'],'_external_paths.py'))
command_paths = findPaths(['dock6'])
dirs['dock6'] = os.path.abspath(os.path.join(\
  os.path.dirname(command_paths['dock6']),'..'))

# Create directories if necessary
labels = {'ligand':os.path.basename(args.ligand)[:-5],
          'receptor':os.path.basename(args.receptor)[:-4]}
labels['complex'] = labels['ligand']+'-'+labels['receptor']
ancg_in = {'ligand':os.path.abspath(args.ligand),
           'receptor_prefix':os.path.abspath(args.receptor[:-4]),
           'grid_prefix':os.path.abspath(args.receptor[:-4]),
           'dir_dock6':dirs['dock6'],
           'outfile_prefix':labels['complex']}

F = open(labels['complex']+'.in','w')
F.write('''
ligand_atom_file                                             ''' + ancg_in['ligand'] + '''
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               no
use_database_filter                                          no
orient_ligand                                                yes
automated_matching                                           yes
receptor_site_file                                           ''' + ancg_in['receptor_prefix'] + '''.sph
max_orientations                                             5000
critical_points                                              no
chemical_matching                                            no
use_ligand_spheres                                           no
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
flexible_ligand                                              yes
user_specified_anchor                                        no
limit_max_anchors                                            no
min_anchor_size                                              40
pruning_use_clustering                                       yes
pruning_max_orients                                          1000
pruning_clustering_cutoff                                    1000
pruning_conformer_score_cutoff                               25.0
use_clash_overlap                                            no
write_growth_tree                                            no
bump_filter                                                  yes
bump_grid_prefix                                             ''' + ancg_in['receptor_prefix'] + '''
max_bumps_anchor                                             12
max_bumps_growth                                             12
score_molecules                                              yes
contact_score_primary                                        no
contact_score_secondary                                      no
grid_score_primary                                           yes
grid_score_secondary                                         no
grid_score_rep_rad_scale                                     1
grid_score_vdw_scale                                         1
grid_score_es_scale                                          1
grid_score_grid_prefix                                       ''' + ancg_in['receptor_prefix'] + '''
multigrid_score_secondary                                    no
dock3.5_score_secondary                                      no
continuous_score_secondary                                   no
descriptor_score_secondary                                   no
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
SASA_descriptor_score_secondary                              no
amber_score_secondary                                        no
minimize_ligand                                              yes
minimize_anchor                                              yes
minimize_flexible_growth                                     yes
use_advanced_simplex_parameters                              no
simplex_max_cycles                                           1
simplex_score_converge                                       0.1
simplex_cycle_converge                                       1.0
simplex_trans_step                                           1.0
simplex_rot_step                                             0.1
simplex_tors_step                                            10.0
simplex_anchor_max_iterations                                500
simplex_grow_max_iterations                                  500
simplex_grow_tors_premin_iterations                          0
simplex_random_seed                                          0
simplex_restraint_min                                        no
atom_model                                                   all
vdw_defn_file                                                ''' + ancg_in['dir_dock6'] + '''/parameters/vdw_AMBER_parm99.defn
flex_defn_file                                               ''' + ancg_in['dir_dock6'] + '''/parameters/flex.defn
flex_drive_file                                              ''' + ancg_in['dir_dock6'] + '''/parameters/flex_drive.tbl
ligand_outfile_prefix                                        ''' + ancg_in['outfile_prefix'] + '''
write_orientations                                           no
num_scored_conformers                                        1000
write_conformations                                          yes
cluster_conformations                                        yes
cluster_rmsd_threshold                                       2.0
rank_ligands                                                 no
''')
F.close()

command = dirs['dock6'] + '/bin/dock6 -i '+labels['complex']+'.in; ' + \
  'gzip '+labels['complex']+'_scored.mol2; ' + \
  'mv '+labels['complex']+'_scored.mol2.gz ancg-'+labels['complex']+'.mol2.gz; ' + \
  'rm '+labels['complex']+'_conformers.mol2; ' + \
  'rm '+labels['complex']+'.in'
print command
os.system(command)
