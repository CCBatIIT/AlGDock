import AlGDock.BindingPMF_plots
import os, shutil, glob

for run_type in ['cool','dock','postprocess','free_energies']:
  self = AlGDock.BindingPMF_plots.BPMF_plots(\
    dir_dock='dock', dir_cool='cool',\
    ligand_tarball='prmtopcrd/ligand.tar.gz', \
    ligand_database='prmtopcrd/ligand.db', \
    forcefield='prmtopcrd/gaff.dat', \
    ligand_prmtop='prmtopcrd/ligand.prmtop', \
    ligand_inpcrd='prmtopcrd/ligand.trans.inpcrd', \
    receptor_tarball='prmtopcrd/receptor.tar.gz', \
    receptor_prmtop='prmtopcrd/receptor.prmtop', \
    receptor_inpcrd='prmtopcrd/receptor.trans.inpcrd', \
    receptor_fixed_atoms='prmtopcrd/receptor.pdb', \
    complex_tarball='prmtopcrd/complex.tar.gz', \
    complex_prmtop='prmtopcrd/complex.prmtop', \
    complex_inpcrd='prmtopcrd/complex.trans.inpcrd', \
    complex_fixed_atoms='prmtopcrd/complex.pdb', \
    score = 'prmtopcrd/anchor_and_grow_scored.mol2', \
    dir_grid='grids', \
    protocol='Adaptive', cool_therm_speed=1.5, dock_therm_speed=1.5,\
    sampler='NUTS', \
    MCMC_moves=1, \
    seeds_per_state=10, steps_per_seed=200,
    sweeps_per_cycle=10, attempts_per_sweep=100, steps_per_sweep=50,
    cool_repX_cycles=2, dock_repX_cycles=3, \
    site='Sphere', site_center=[1.91650, 1.91650, 1.91650], site_max_R=0.01, \
    site_density=10., \
    phases=['NAMD_Gas','NAMD_GBSA'], \
    cores=-1, \
    rmsd=True)
  self._run(run_type)
  del self

# site='Sphere', site_center=[1.91650, 1.91650, 1.91650], site_max_R=0.01, \
# score = 'prmtopcrd/anchor_and_grow_scored.mol2', \