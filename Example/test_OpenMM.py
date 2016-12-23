# The example is 1of6

run_types = [None]
# run_types = ['all']
# run_types = ['cool']
# run_types = ['cool','dock','postprocess','free_energies']

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
    pose=-1, \
    rmsd=True, \
    dir_grid='grids', \
    protocol='Geometric', cool_therm_speed=0.05, dock_therm_speed=0.013,\
    sampler='HMC', \
    MCMC_moves=1, \
    seeds_per_state=10, steps_per_seed=200, darts_per_seed=0, \
    sweeps_per_cycle=25, attempts_per_sweep=100, \
    steps_per_sweep=50, darts_per_sweep=0, \
    cool_repX_cycles=3, dock_repX_cycles=4, \
    site='Sphere', site_center=[1.74395, 1.74395, 1.74395], \
    site_max_R=0.6, \
    site_density=10., \
    phases=phases, \
    cores=-1, \
    random_seed=-1)
  self._run(run_type)

# Good thermodynamic speeds for testing:
# protocol='Adaptive', cool_therm_speed=0.5, dock_therm_speed=0.5

#  self._set_universe_evaluator(self._lambda(1.0,process='cool'))
#  (en,grad)=self.universe.energyAndGradients()

#  import time
#  t = time.time()
#  for r in range(1000):
#    (en,grad)=self.universe.energyAndGradients()
#  print en
#  print grad
#  print 'Time for MMTK:', time.time() - t
#
#  from AlGDock.ForceFields.OpenMM.OpenMM import OpenMMForceField
#  self._forceFields['OpenMM'] = OpenMMForceField(self._FNs['prmtop']['L'], self.molecule.prmtop_atom_order, self.molecule.inv_prmtop_atom_order, implicitSolvent='OpenMM_Gas')
#  self.universe.setForceField(self._forceFields['OpenMM'])
#  (en,grad)=self.universe.energyAndGradients()
#
#  import time
#  t = time.time()
#  for r in range(1000):
#    (en,grad)=self.universe.energyAndGradients()
#  print en
#  print grad
#  print 'Time for OpenMM Gas:', time.time() - t
#
from AlGDock.ForceFields.OpenMM.OpenMM import OpenMMForceField
self._forceFields['OpenMM'] = OpenMMForceField(self._FNs['prmtop']['L'], self.molecule.prmtop_atom_order, self.molecule.inv_prmtop_atom_order, implicitSolvent='OpenMM_OBC2')
self.universe.setForceField(self._forceFields['OpenMM'])
(en,grad)=self.universe.energyAndGradients()

print en
print grad.array
#
#  import time
#  t = time.time()
#  for r in range(1000):
#    (en,grad)=self.universe.energyAndGradients()
#  print en
#  print grad
#  print 'Time for OpenMM GBSA:', time.time() - t

# This works but is much slower.
# Time for MMTK: 0.0425839424133
# Time for OpenMM Gas: 1.64846014977
# Time for OpenMM GBSA: 2.23015904427
# It might be better to recode GB
# https://github.com/pandegroup/openmm/blob/master/platforms/cpu/src/CpuGBSAOBCForce.cpp
