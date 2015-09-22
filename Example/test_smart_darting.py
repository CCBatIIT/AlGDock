import AlGDock.BindingPMF_plots
import os, shutil, glob
import MMTK
from MMTK import Units, Configuration
import numpy as np

R = 8.3144621*Units.J/Units.mol/Units.K

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
  phases=['NAMD_Gas','NAMD_GBSA'], \
  cores=-1, \
  rmsd=True)

# These are cartesian coordinates of docked poses
(confs, Es) = self._get_confs_to_rescore(minimize=False)

# Use a force field without grids
self._set_universe_evaluator(self._lambda(0.0, 'dock', site=False, MM=True))

#  conf_energies = []
#  for conf in confs:
#    self.universe.setConfiguration(MMTK.Configuration(self.universe,conf))
#    conf_energies.append(self.universe.energy())
#  # Probabilty of jumping to a conformation k is proportional to exp(-E/RT).
#  weights = np.exp(-np.array(conf_energies)/(R*1000.0))
#  weights = weights/sum(weights)
#
#  import AlGDock.RigidBodies
#  BAT_util = AlGDock.RigidBodies.identifier(self.universe, self.molecule)
#  confs_BAT = [np.array(BAT_util.BAT(Cartesian=confs[n],extended=True)) \
#    for n in range(len(confs))]
#  # Perturb external coordinates and primary torsions
#  BAT_to_perturb = range(6)+list(len(confs_BAT[0])-BAT_util.ntorsions+np.array(sorted(list(set(BAT_util._firstTorsionInd)))))
#  # darts[j][k] will jump from conformation j to conformation k
#  darts = [[confs_BAT[j][BAT_to_perturb]-confs_BAT[k][BAT_to_perturb] \
#    for j in range(len(confs))] for k in range(len(confs))]
#
#  def closest_pose(conf, confs):
#    return np.argmin(np.array([np.sqrt(((\
#      confs[c][self.molecule.heavy_atoms,:] - conf)**2).sum()\
#      /self.molecule.nhatoms) for c in range(len(confs))]))
#
#  def p_attempt(weights,dart_from,dart_to):
#    return weights[dart_to]/(1.-weights[dart_from])
#
#  from NUTS import NUTSIntegrator  # @UnresolvedImport
#  NUTS_sampler = NUTSIntegrator(self.universe)
#
#  acc = 0.
#  trials = 25
#  energies = []
#  closest_poses = []
#  for rep in range(10):
#    NUTS_sampler(steps=50,T=300)
#    xo_Cartesian = self.universe.copyConfiguration()
#    xo_BAT = np.array(BAT_util.BAT(extended=True))
#    eo = self.universe.energy()
#    closest_pose_o = closest_pose(xo_Cartesian.array[self.molecule.heavy_atoms,:], confs)
#
#    for t in range(trials):
#      # Choose a pose to dart towards
#      dart_towards = closest_pose_o
#      while dart_towards==closest_pose_o:
#        dart_towards = np.random.choice(len(weights), p=weights)
#      # Generate a trial move
#      xn_BAT = np.copy(xo_BAT)
#      xn_BAT[BAT_to_perturb] = xo_BAT[BAT_to_perturb] + darts[closest_pose_o][dart_towards]
#      xn_Cartesian = Configuration(self.universe, BAT_util.Cartesian(xn_BAT)) # Sets the universe
#      en = self.universe.energy()
#      closest_pose_n = closest_pose(\
#        self.universe.configuration().array[self.molecule.heavy_atoms,:], confs)
#
#      if (en<eo) or (np.random.random()<np.exp(-(en-eo)/self.RT)*\
#          p_attempt(weights,closest_pose_n,closest_pose_o)/\
#          p_attempt(weights,closest_pose_o,dart_towards)):
#        
#        if en>(1000):
#          print 'eo: %f, en: %f'%(eo,en)
#          print 'closest_pose_o %d, closest_pose_n %d'%(\
#            closest_pose_o,closest_pose_n)
#          print 'p_attempt, forward %f, backwards %f'%(\
#            p_attempt(weights,closest_pose_o,dart_towards), \
#            p_attempt(weights,closest_pose_n,closest_pose_o))
#          raise Exception('High energy pose!')
#        acc += 1
#        
#        xo_Cartesian = xn_Cartesian
#        xo_BAT = xn_BAT
#        eo = en
#        closest_pose_o = closest_pose_n
#        
#      else:
#        self.universe.setConfiguration(xo_Cartesian)
#
#    energies.append(eo)
#    closest_poses.append(closest_pose_o)
#
#  print 'Acceptance rate: ', acc/trials/10.
#  #  import matplotlib.pyplot as plt
#  #  plt.plot(energies,'.-')
#  #  plt.show()

BPMF = self

from AlGDock.Integrators.SmartDarting.SmartDarting import SmartDartingIntegrator
self = SmartDartingIntegrator(BPMF.universe, BPMF.molecule, False)
print self.set_confs(confs)
print 'Coordinates to perturb:', ' '.join(['%d'%d for d in self._BAT_to_perturb])
print 'Acceptance rate, no extended coordinates: ', self(T=300, ntrials=100)[2]

del self
self = SmartDartingIntegrator(BPMF.universe, BPMF.molecule, True)
print self.set_confs(confs)
print 'Coordinates to perturb:', ' '.join(['%d'%d for d in self._BAT_to_perturb])
print 'Acceptance rate:    extended coordinates:', self(T=300, ntrials=100)[2]
