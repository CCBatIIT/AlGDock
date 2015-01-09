import AlGDock.BindingPMF
import os, shutil, glob

for run_type in ['cool','dock','postprocess','free_energies']:
  self = AlGDock.BindingPMF.BPMF(\
    dir_dock='dock', dir_cool='cool',\
    ligand_database='prmtopcrd/ligand.db', \
    forcefield='prmtopcrd/gaff.dat', \
    ligand_prmtop='prmtopcrd/ligand.prmtop', \
    ligand_inpcrd='prmtopcrd/ligand.trans.inpcrd', \
    receptor_prmtop='prmtopcrd/receptor.prmtop', \
    receptor_inpcrd='prmtopcrd/receptor.trans.inpcrd', \
    receptor_fixed_atoms='prmtopcrd/receptor.pdb', \
    complex_prmtop='prmtopcrd/complex.prmtop', \
    complex_inpcrd='prmtopcrd/complex.trans.inpcrd', \
    complex_fixed_atoms='prmtopcrd/complex.pdb', \
    dir_grid='grids', \
    protocol='Adaptive', cool_therm_speed=1.25, dock_therm_speed=1.25,\
    no_protocol_refinement=True,\
    sampler='NUTS', \
    cool_seeds_per_state=20, dock_seeds_per_state=10,\
    steps_per_seed=250,
    sweeps_per_cycle=5, steps_per_sweep=250,
    cool_repX_cycles=3, dock_repX_cycles=4, \
    site='Sphere', site_center=[1.91650, 1.91650, 1.91650], \
    site_max_R=0.01, site_density=10., \
    phases=['NAMD_Gas'], \
    MCMC_moves=1, \
    cores=-1, \
    score='prmtopcrd/anchor_and_grow_scored.mol2', \
    rmsd=True, \
    run_type=run_type)
  del self
