# The example is 1of6

# run_types = [None]
run_types = ['all']
# run_types = ['cool']
# run_types = ['cool','dock','postprocess','free_energies']
# run_types = ['timed']

import AlGDock.BindingPMF_plots
import os, shutil, glob

#  from AlGDock.BindingPMF_arguments import allowed_phases
#  phases = allowed_phases

phases = ['NAMD_Gas', 'NAMD_OBC']

self = None
for run_type in run_types:
  del self
  self = AlGDock.BindingPMF_plots.BPMF_plots(\
    dir_dock='dock', dir_cool='cool',\
    ligand_tarball='prmtopcrd/ligand.tar.gz', \
    ligand_database='ligand.db', \
    forcefield='prmtopcrd/gaff2.dat', \
    ligand_prmtop='ligand.prmtop', \
    ligand_inpcrd='ligand.trans.inpcrd', \
    ligand_mol2='ligand.mol2', \
    ligand_rb='ligand.rb', \
    receptor_tarball='prmtopcrd/receptor.tar.gz', \
    receptor_prmtop='receptor.prmtop', \
    receptor_inpcrd='receptor.trans.inpcrd', \
    receptor_fixed_atoms='receptor.pdb', \
    complex_tarball='prmtopcrd/complex.tar.gz', \
    complex_prmtop='complex.prmtop', \
    complex_inpcrd='complex.trans.inpcrd', \
    complex_fixed_atoms='complex.pdb', \
    score = 'prmtopcrd/anchor_and_grow_scored.mol2', \
    pose=-1, \
    rmsd=True, \
    dir_grid='grids', \
    protocol='Adaptive', cool_therm_speed=5.0, dock_therm_speed=0.5, \
    T_HIGH=600.0, T_SIMMIN=300.0, T_TARGET=300.0, \
    sampler='HMC', \
    MCMC_moves=1, \
    sampling_importance_resampling = True, \
    solvation = 'Fractional', \
    seeds_per_state=10, steps_per_seed=200, darts_per_seed=0, \
    sweeps_per_cycle=25, attempts_per_sweep=100, \
    steps_per_sweep=100, darts_per_sweep=0, \
    cool_repX_cycles=3, dock_repX_cycles=4, \
    site='Cylinder', site_center=[1.74395, 1.74395, 1.44395], \
    site_direction=[0, 0, 1.0], \
    site_max_Z=2.04395, \
    site_max_R=0.6, \
    phases=phases, \
    cores=-1, \
    random_seed=-1, \
    max_time=240, \
    keep_intermediate=True)
  self._run(run_type)

# To test timed runs:
#    run_type='timed',
#    max_time=0.25)

# To use spherical binding site
#    site='Sphere', site_center=[1.74395, 1.74395, 1.74395], \
#    site_max_R=0.6, \
#    site_density=10., \

# To use cylindrical binding site
#    site='Cylinder', site_center=[1.74395, 1.74395, 1.44395], \
#    site_direction=[0, 0, 1.0], \
#    site_max_Z=2.04395, \
#    site_max_R=0.6, \

# Reasonable testing protocols:
# protocol='Geometric', cool_therm_speed=0.05, dock_therm_speed=0.013,\
# protocol='Adaptive', cool_therm_speed=0.5, dock_therm_speed=0.5

# Reasonable sampler options:
# sampler='HMC', \
