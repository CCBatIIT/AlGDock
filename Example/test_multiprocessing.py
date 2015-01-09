import AlGDock.BindingPMF
import os, shutil, glob
import time

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
  rmsd=True, \
  dir_grid='grids', \
  protocol='Adaptive', cool_therm_speed=1.0, dock_therm_speed=1.5,\
  no_protocol_refinement=True,\
  sampler='NUTS', seeds_per_state=10, steps_per_seed=500, sweeps_per_cycle=50, \
  steps_per_sweep=100, cool_repX_cycles=3, dock_repX_cycles=2, \
  site='Sphere', site_center=[1.91650, 1.91650, 1.91650], \
  site_max_R=0.01, site_density=10., \
  receptor_NAMD_GBSA=[-9053.10071], receptor_NAMD_Gas=[0.00000], \
  MCMC_moves=1, \
  cores=32, \
  phases=['NAMD_Gas','NAMD_GBSA'],
  run_type='all')

for (dock_cycle,cores) in zip(range(3,9),[32,16,8,4,2,1]):
  start_time = time.time()
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
    rmsd=True, \
    dir_grid='grids', \
    protocol='Adaptive', cool_therm_speed=1.0, dock_therm_speed=1.5,\
    no_protocol_refinement=True,\    
    sampler='NUTS', seeds_per_state=10, steps_per_seed=500, sweeps_per_cycle=50, \
    steps_per_sweep=100, cool_repX_cycles=3, dock_repX_cycles=dock_cycle, \
    site='Sphere', site_center=[1.91650, 1.91650, 1.91650], \
    site_max_R=0.01, site_density=10., \
    receptor_NAMD_GBSA=[-9053.10071], receptor_NAMD_Gas=[0.00000], \
    MCMC_moves=1, \
    cores=cores, \
    phases=['NAMD_Gas','NAMD_GBSA'],
    run_type='all')
  del self
  print 'Time for one docking cycle with %d cores: %f s'%(\
    cores, time.time()-start_time)

