# The example is 1of6

import AlGDock.BindingPMF_plots
import os, shutil, glob

self = AlGDock.BindingPMF_plots.BPMF_plots(\
  dir_dock='dock', dir_cool='cool',\
  ligand_tarball='prmtopcrd/ligand.tar.gz', \
  ligand_database='ligand.db', \
  forcefield='prmtopcrd/gaff.dat', \
  ligand_prmtop='ligand.prmtop', \
  ligand_inpcrd='ligand.trans.inpcrd', \
  receptor_tarball='prmtopcrd/receptor.tar.gz', \
  receptor_prmtop='receptor.prmtop', \
  receptor_inpcrd='receptor.trans.inpcrd', \
  receptor_fixed_atoms='receptor.pdb', \
  complex_tarball='prmtopcrd/complex.tar.gz', \
  complex_prmtop='complex.prmtop', \
  complex_inpcrd='complex.trans.inpcrd', \
  complex_fixed_atoms='complex.pdb', \
  score = 'prmtopcrd/anchor_and_grow_scored.mol2', \
  dir_grid='grids', \
  protocol='Adaptive', cool_therm_speed=1.5, dock_therm_speed=1.5,\
  sampler='NUTS', \
  MCMC_moves=1, \
  seeds_per_state=10, steps_per_seed=200, darts_per_seed=10, \
  sweeps_per_cycle=25, attempts_per_sweep=100, \
  steps_per_sweep=50, darts_per_sweep=5, \
  cool_repX_cycles=3, dock_repX_cycles=4, \
  site='Sphere', site_center=[1.74395, 1.74395, 1.74395], site_max_R=0.6, \
  site_density=10., \
  phases=['NAMD_Gas','NAMD_GBSA'], \
  cores=-1, \
  rmsd=True)

AlGDock_path = os.path.join(os.path.dirname(AlGDock.__file__),'BindingPMF.py')
import time
energyString = 'AlGDock modified %s\n'%time.ctime(os.path.getmtime(AlGDock_path))

for a in [0.75, 0.8]:
  self._set_universe_evaluator(self._lambda(a=a))
  terms = self.universe.energyTerms()
  for key in terms.keys():
    energyString += '%s: %f\t'%(key,terms[key])
  energyString += '\n'

F = open('grid_energies.txt','w')
F.write(energyString)
F.close()
