# This module implements a Smart Darting "integrator"

from MMTK import Configuration, Dynamics, Environment, Features, Trajectory, Units
import MMTK_dynamics
import numpy as np
import AlGDock.RigidBodies

import random

R = 8.3144621*Units.J/Units.mol/Units.K

#
# Smart Darting integrator
#
class SmartDartingIntegrator(Dynamics.Integrator):
  def __init__(self, universe, molecule, extended, confs=None, **options):
    """
    confs - configurations to dart to
    extended - whether or not to use external coordinates
    """
    Dynamics.Integrator.__init__(self, universe, options)
    # Supported features: none for the moment, to keep it simple
    self.features = []

    self.molecule = molecule
    self.extended = extended
    
    # Converter between Cartesiand BAT coordinates
    self._BAT_util = AlGDock.RigidBodies.identifier(self.universe, self.molecule)
    # BAT coordinates to perturb: external coordinates and primary torsions
    dof = self.universe.configuration().array.shape[0]*3 if extended \
      else self.universe.configuration().array.shape[0]*3 - 6
    self._BAT_to_perturb = range(6) if extended else []
    self._BAT_to_perturb += list(dof - self._BAT_util.ntorsions + \
      np.array(sorted(list(set(self._BAT_util._firstTorsionInd)))))

    if confs is None:
      self.confs = None
    else:
      self.set_confs(confs)

  def set_confs(self, confs, rmsd_threshold=0.1, period_threshold=0.25, \
      append=False):

    nconfs_attempted = len(confs)
    if append and (self.confs is not None):
      nconfs_o = len(self.confs)
      confs = confs + self.confs
    else:
      nconfs_o = 0
    
    # Calculate energies of all configurations
    conf_energies = []
    for conf in confs:
      self.universe.setConfiguration(Configuration(self.universe,conf))
      conf_energies.append(self.universe.energy())
    self.universe.setConfiguration(Configuration(self.universe,confs[0]))

    # Sort by increasing energy
    conf_energies, confs = (list(l) \
      for l in zip(*sorted(zip(conf_energies, confs), key=lambda p:p[0])))
    
    if self.extended:
      # Keep only unique configurations, using rmsd as a threshold
      inds_to_keep = [0]
      for j in range(len(confs)):
        min_rmsd = np.min([((confs[j][self.molecule.heavy_atoms,:] - \
          confs[k][self.molecule.heavy_atoms,:])**2).sum()/self.molecule.nhatoms \
          for k in inds_to_keep])
        if min_rmsd>rmsd_threshold:
          inds_to_keep.append(j)
      confs = [confs[i] for i in inds_to_keep]
      conf_energies = [conf_energies[i] for i in inds_to_keep]

    confs_BAT = [np.array(self._BAT_util.BAT(Cartesian=confs[n], \
      extended=self.extended)) for n in range(len(confs))]
    confs_BAT_tp = [confs_BAT[c][self._BAT_to_perturb] \
      for c in range(len(confs_BAT))]

    if not self.extended:
      # Keep only unique configurations, using period as a threshold
      inds_to_keep = [0]
      for j in range(len(confs)):
        diffs = np.array([(confs_BAT_tp[j] - confs_BAT_tp[k]) \
          for k in inds_to_keep])/(2*np.pi)
        diffs = np.min([diffs,1-diffs],0)
        if (np.max(np.abs(diffs),1)>period_threshold).all():
          inds_to_keep.append(j)
      confs = [confs[i] for i in inds_to_keep]
      confs_BAT = [confs_BAT[i] for i in inds_to_keep]
      confs_BAT_tp = [confs_BAT_tp[i] for i in inds_to_keep]
      conf_energies = [conf_energies[i] for i in inds_to_keep]

    self.confs = confs
    self.confs_ha = [self.confs[c][self.molecule.heavy_atoms,:] \
      for c in range(len(self.confs))]

    # Probabilty of jumping to a conformation k is proportional to exp(-E/R*1000).
    weights = np.exp(-np.array(conf_energies)/(R*1000.))
    self.weights = weights/sum(weights)

    self.confs_BAT = confs_BAT
    self.confs_BAT_tp = confs_BAT_tp
    # self.darts[j][k] will jump from conformation j to conformation k
    self.darts = [[self.confs_BAT[j][self._BAT_to_perturb] - \
      self.confs_BAT[k][self._BAT_to_perturb] \
      for j in range(len(confs))] for k in range(len(confs))]
    
    return '  Set smart darting configurations. ' + \
      'Started with %d, attempted %d, ended with %d configurations.'%(\
      nconfs_o,nconfs_attempted,len(self.confs))

  def _closest_pose_Cartesian(self, conf_ha):
    # Closest pose has smallest sum of square distances between heavy atom coordinates
    return np.argmin(np.array([((self.confs_ha[c] - conf_ha)**2).sum() \
      for c in range(len(self.confs))]))

  def _closest_pose_BAT(self, conf_BAT_tp):
    # Closest pose has smallest sum of square distances between torsion angles
    # For only torsion angles, differences in units of periods (between 0 and 1)
    diffs = np.abs(np.array([self.confs_BAT_tp[c] - conf_BAT_tp \
      for c in range(len(self.confs_BAT_tp))]))/(2*np.pi)
    # Wraps around the period
    diffs = np.min([diffs,1-diffs],0)
    return np.argmin(np.sum(diffs**2,1))

  def _p_attempt(self, dart_from, dart_to):
    return self.weights[dart_to]/(1.-self.weights[dart_from])

  def __call__(self, **options):
    if (self.confs is None) or len(self.confs)<3:
      return ([self.universe.copyConfiguration()], [self.universe.energy()], 0.0, 0.0)
    
    # Process the keyword arguments
    self.setCallOptions(options)
    # Check if the universe has features not supported by the integrator
    Features.checkFeatures(self, self.universe)
  
    RT = R*self.getOption('T')
    ntrials = self.getOption('ntrials')

    acc = 0.
    energies = []
    closest_poses = []

    xo_Cartesian = self.universe.copyConfiguration()
    xo_BAT = np.array(self._BAT_util.BAT(extended=self.extended))
    eo = self.universe.energy()
    if self.extended:
      closest_pose_o = self._closest_pose_Cartesian(\
        xo_Cartesian.array[self.molecule.heavy_atoms,:])
    else:
      closest_pose_o = self._closest_pose_BAT(xo_BAT[self._BAT_to_perturb])
      
    for t in range(ntrials):
      # Choose a pose to dart towards
      dart_towards = closest_pose_o
      while dart_towards==closest_pose_o:
        dart_towards = np.random.choice(len(self.weights), p=self.weights)
      # Generate a trial move
      xn_BAT = np.copy(xo_BAT)
      xn_BAT[self._BAT_to_perturb] = xo_BAT[self._BAT_to_perturb] + self.darts[closest_pose_o][dart_towards]
      xn_Cartesian = Configuration(self.universe,self._BAT_util.Cartesian(xn_BAT)) # Sets the universe
      en = self.universe.energy()
      if self.extended:
        closest_pose_n = self._closest_pose_Cartesian(\
          xn_Cartesian.array[self.molecule.heavy_atoms,:])
      else:
        closest_pose_n = self._closest_pose_BAT(xn_BAT[self._BAT_to_perturb])
        
      # Accept or reject the trial move
      if (en<eo) or (np.random.random()<np.exp(-(en-eo)/RT)*\
          self._p_attempt(closest_pose_n,closest_pose_o)/\
          self._p_attempt(closest_pose_o,dart_towards)):
        
        if en>(1000):
          print 'eo: %f, en: %f'%(eo,en)
          print 'closest_pose_o %d, closest_pose_n %d'%(\
            closest_pose_o,closest_pose_n)
          print '_p_attempt, forward %f, backwards %f'%(\
            self._p_attempt(closest_pose_o,dart_towards), \
            self._p_attempt(closest_pose_n,closest_pose_o))
          raise Exception('High energy pose!')
        
        xo_Cartesian = xn_Cartesian
        xo_BAT = xn_BAT
        eo = en
        closest_pose_o = closest_pose_n
        acc += 1
      else:
        self.universe.setConfiguration(xo_Cartesian)

    return ([np.copy(xo_Cartesian)], [en], float(acc)/float(ntrials), 0.0)