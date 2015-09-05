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
  seeds_per_state=10, steps_per_seed=200,
  sweeps_per_cycle=10, attempts_per_sweep=100, steps_per_sweep=50,
  cool_repX_cycles=2, dock_repX_cycles=3, \
  site='Sphere', site_center=[1.91650, 1.91650, 1.91650], site_max_R=0.01, \
  site_density=10., \
  phases=['NAMD_Gas','NAMD_GBSA','GBSA'], \
  cores=-1, \
  rmsd=True)

prmtopFN = self._FNs['prmtop']['R']
inpcrdFN = self._FNs['inpcrd']['R']

import simtk.openmm.app
prmtop = simtk.openmm.app.AmberPrmtopFile(prmtopFN)
inpcrd = simtk.openmm.app.AmberInpcrdFile(inpcrdFN)
system_Gas = prmtop.createSystem(nonbondedMethod=simtk.openmm.app.NoCutoff, \
  constraints=HBonds)
system_GBSA = prmtop.createSystem(nonbondedMethod=simtk.openmm.app.NoCutoff, \
  constraints=HBonds, implicitSolvent=simtk.openmm.app.OBC2)
integrator_Gas = simtk.openmm.LangevinIntegrator(300*simtk.unit.kelvin, \
  1/simtk.unit.picosecond, 0.002*simtk.unit.picoseconds)
# OpenMM GBSA model can be None, HCT (igb=1), OBC1 (igb=2), OBC2 (igb=5), GBn (igb=7), GBn2 (igb=8).
integrator_GBSA = simtk.openmm.LangevinIntegrator(300*simtk.unit.kelvin, \
  1/simtk.unit.picosecond, 0.002*simtk.unit.picoseconds)
simulation_Gas = simtk.openmm.app.Simulation(prmtop.topology, system_Gas, integrator_Gas)
simulation_Gas.context.setPositions(inpcrd.positions)
simulation_GBSA = simtk.openmm.app.Simulation(prmtop.topology, system_GBSA, integrator_GBSA)
simulation_GBSA.context.setPositions(inpcrd.positions)
print simulation_Gas.context.getState(getEnergy=True).getPotentialEnergy()
print simulation_GBSA.context.getState(getEnergy=True).getPotentialEnergy()
