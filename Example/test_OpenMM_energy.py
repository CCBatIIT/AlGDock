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

import simtk.openmm.app

def setup_OMM_simulation(moiety):
  prmtopFN = self._FNs['prmtop'][moiety]
  inpcrdFN = self._FNs['inpcrd'][moiety]
  
  prmtop = simtk.openmm.app.AmberPrmtopFile(prmtopFN)
  inpcrd = simtk.openmm.app.AmberInpcrdFile(inpcrdFN)
  # A typical constraint is simtk.openmm.app.HBonds
  # OpenMM GBSA model can be None, HCT (igb=1), OBC1 (igb=2), OBC2 (igb=5), GBn (igb=7), GBn2 (igb=8).
  # A typical implicitSolvent is simtk.openmm.app.OBC2
  OMM_system = prmtop.createSystem(nonbondedMethod=simtk.openmm.app.NoCutoff, \
    constraints=None, implicitSolvent=None)
  dummy_integrator = simtk.openmm.LangevinIntegrator(300*simtk.unit.kelvin, \
    1/simtk.unit.picosecond, 0.002*simtk.unit.picoseconds)
  return simtk.openmm.app.Simulation(prmtop.topology, OMM_system, dummy_integrator)

# Test ligand gas energies
OMM_L_simulation = setup_OMM_simulation('L')
Es_L = []
confs = self.confs['dock']['samples'][-1][-1]
for conf in confs:
  OMM_L_simulation.context.setPositions(conf[self.molecule.prmtop_atom_order,:])
  Es_L.append(OMM_L_simulation.context.getState(getEnergy=True).getPotentialEnergy()/simtk.unit.kilojoule*simtk.unit.mole)

import numpy as np
Es_L = np.array(Es_L)
print 'Energy differences between OpenMM and NAMD Gas energies'
print Es_L-self.dock_Es[-1][-1]['LNAMD_Gas'][:,-1]

# Test complex energies
OMM_R_simulation = setup_OMM_simulation('R')
OMM_R_simulation.context.setPositions(self.confs['receptor'])
Es_R = OMM_R_simulation.context.getState(getEnergy=True).getPotentialEnergy()/simtk.unit.kilojoule*simtk.unit.mole

OMM_RL_simulation = setup_OMM_simulation('RL')
Es_RL = []
for conf in confs:
  full_conf = np.vstack((self.confs['receptor'],conf[self.molecule.prmtop_atom_order,:]))
  OMM_RL_simulation.context.setPositions(full_conf)
  Es_RL.append(OMM_RL_simulation.context.getState(getEnergy=True).getPotentialEnergy()/simtk.unit.kilojoule*simtk.unit.mole)
Es_RL = np.array(Es_RL)

Psi = Es_RL - Es_R - Es_L
print 'Difference between OpenMM and NAMD Gas interaction energies'
print Psi/self.RT_TARGET-self.stats_RL['Psi_NAMD_Gas'][-1]
