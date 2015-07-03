import AlGDock.BindingPMF_plots
import os, shutil, glob

self = AlGDock.BindingPMF_plots.BPMF_plots(\
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

# These are cartesian coordinates of docked poses
(confs, Es) = self._get_confs_to_rescore(minimize=False)

# Convert to bond-angle-torsion (BAT)
# Figure out how to move from one docked pose to another. This is your MCMC move set.

import AlGDock.BAT
c = AlGDock.BAT.converter()
c.BAT() # Calculates BAT coordinates

# Basically need to loope over docked poses to get "darts."
# _get_confs_to_rescore will get you docked poses in AlGDock

# What pose are you in?
# Your current configuration is: self.universe.configuration().array
# Get rmsd values between your current pose and pose of interest. Here a quick one-liner that will do it:
# np.array([np.sqrt(((conf_c - confs[c][self.molecule.heavy_atoms,:])**2).sum()/self.molecule.nhatoms) for c in range(len(confs))])

# Select a new pose to dart to with a probability that satisfies detailed balance.
# See _MC_translate_rotate in AlGDock.BindingPMF.BPMF for an example for how to do an MCMC move.
