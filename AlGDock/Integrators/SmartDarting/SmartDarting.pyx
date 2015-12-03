# This module implements a Smart Darting "integrator"
# It is a little bit faster than SmartDarting.py,
# but does not seem to work on the cluster.

import cython

import numpy as np
cimport numpy as np

cimport MMTK_trajectory_generator
from MMTK import Units
from MMTK import Features

import MMTK_trajectory
import MMTK_forcefield

cdef extern from "stdlib.h":

    ctypedef long size_t
    cdef void *malloc(size_t size)
    cdef void free(void *ptr)

from MMTK.ParticleProperties import Configuration, ParticleVector

include "MMTK/python.pxi"
include "MMTK/numeric.pxi"
include "MMTK/core.pxi"
include "MMTK/universe.pxi"
include "MMTK/trajectory.pxi"
include "MMTK/forcefield.pxi"

from cpython cimport bool
cimport BAT

R = 8.3144621*Units.J/Units.mol/Units.K

#
# Smart Darting integrator
#
cdef class SmartDartingIntegrator(MMTK_trajectory_generator.EnergyBasedTrajectoryGenerator):

  cdef BAT.converter _BAT_util
  cdef readonly object _BAT_to_perturb
  cdef object molecule, confs, confs_ha, \
    confs_BAT, confs_BAT_tp, darts
  cdef bool extended
  cdef np.ndarray weights
  cdef energy_data energy

  default_options = {'first_step': 0, 'steps': 100, 'delta_t': 1.*Units.fs,
                     'background': False, 'threads': None,
                     'actions': []}

  available_data = ['configuration', 'energy', 'time']

  restart_data = ['configuration', 'velocities', 'energy']

  def __init__(self, universe, molecule, extended, confs=None, **options):
    """
    confs - configurations to dart to
    extended - whether or not to use external coordinates
    """
    MMTK_trajectory_generator.EnergyBasedTrajectoryGenerator.__init__(
        self, universe, options, "SmartDarting integrator")
    # Supported features: none for the moment, to keep it simple
    self.features = []

    self.molecule = molecule
    self.extended = extended
    
    # Converter between Cartesian and BAT coordinates
    self._BAT_util = BAT.converter(self.universe, self.molecule)
    # BAT coordinates to perturb: external coordinates and primary torsions
    dof = self.universe.configuration().array.shape[0]*3 if extended \
      else self.universe.configuration().array.shape[0]*3 - 6
    self._BAT_to_perturb = range(6) if extended else []
    self._BAT_to_perturb += self._BAT_util.getFirstTorsionInds(extended)

    if confs is None:
      self.confs = None
    else:
      self.set_confs(confs)

  def set_confs(self, confs, rmsd_threshold=0.05, period_threshold=0.33, \
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

    # Sort by increasing energy
    conf_energies, confs = (list(l) \
      for l in zip(*sorted(zip(conf_energies, confs), key=lambda p:p[0])))

    # Only keep configurations with energy with 15 kJ/mol of the lowest energy
    confs = [confs[i] for i in range(len(confs)) \
      if (conf_energies[i]-conf_energies[0])<15.]
    conf_energies = [conf_energies[i] for i in range(len(confs)) \
      if (conf_energies[i]-conf_energies[0])<15.]

    if self.extended:
      # Keep only unique configurations, using rmsd as a threshold
      inds_to_keep = [0]
      for j in range(len(confs)):
        min_rmsd = np.min([np.sqrt((confs[j][self.molecule.heavy_atoms,:] - \
          confs[k][self.molecule.heavy_atoms,:])**2).sum()/self.molecule.nhatoms \
          for k in inds_to_keep])
        if min_rmsd>rmsd_threshold:
          inds_to_keep.append(j)
      confs = [confs[i] for i in inds_to_keep]
      conf_energies = [conf_energies[i] for i in inds_to_keep]

    confs_BAT = [self._BAT_util.BAT(confs[n], extended=self.extended) \
      for n in range(len(confs))]
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

    if len(self.confs)>1:
      # Probabilty of jumping to a conformation k
      # is proportional to exp(-E/(R*600.)).
      logweight = np.array(conf_energies)/(R*600.)
      weights = np.exp(-logweight+min(logweight))
      self.weights = weights/sum(weights)

      self.confs_BAT = confs_BAT
      self.confs_BAT_tp = confs_BAT_tp
      # self.darts[j][k] will jump from conformation j to conformation k
      self.darts = [[self.confs_BAT[j][self._BAT_to_perturb] - \
        self.confs_BAT[k][self._BAT_to_perturb] \
        for j in range(len(confs))] for k in range(len(confs))]

    # Set the universe to the lowest-energy configuration
    self.universe.setConfiguration(Configuration(self.universe,np.copy(confs[0])))

    return '  started with %d, attempted %d, ended with %d smart darting targets'%(\
      nconfs_o,nconfs_attempted,len(self.confs))

  # TODO: Try Cython loops in the following functions:

  # - No bound checks on index operations
  # - No support for negative indices
  # - Division uses C semantics
  @cython.boundscheck(False)
  @cython.wraparound(False)
  @cython.cdivision(True)
  cdef _closest_pose_Cartesian(self, conf_ha):
    # Closest pose has smallest sum of square distances between heavy atom coordinates
    return np.argmin(np.array([((self.confs_ha[c] - conf_ha)**2).sum() \
      for c in range(len(self.confs))]))

  # - No bound checks on index operations
  # - No support for negative indices
  # - Division uses C semantics
  @cython.boundscheck(False)
  @cython.wraparound(False)
  @cython.cdivision(True)
  cdef _closest_pose_BAT(self, conf_BAT_tp):
    # Closest pose has smallest sum of square distances between torsion angles
    # For only torsion angles, differences in units of periods (between 0 and 1)
    cdef np.ndarray[np.float_t, ndim=2] diffs
    diffs = np.abs(np.array([self.confs_BAT_tp[c] - conf_BAT_tp \
      for c in range(len(self.confs_BAT_tp))]))/(2*np.pi)
    # Wraps around the period
    diffs = np.min([diffs,1-diffs],0)
    return np.argmin(np.sum(diffs**2,1))

  def __call__(self, **options):
    if (self.confs is None) or len(self.confs)<2:
      return ([self.universe.configuration().array], [self.universe.energy()], 0, 0, 0.0)
    
    # Process the keyword arguments
    self.setCallOptions(options)
    # Check if the universe has features not supported by the integrator
    Features.checkFeatures(self, self.universe)
  
    cdef float RT = R*self.getOption('T')
    cdef int ntrials = self.getOption('ntrials')

    # Seed the random number generator
    if 'random_seed' in self.call_options.keys():
      np.random.seed(self.getOption('random_seed'))

    cdef float acc = 0.
    energies = []
    closest_poses = []

    cdef np.ndarray[np.float_t, ndim=2] xo_Cartesian, xn_Cartesian
    cdef np.ndarray[np.float_t] xo_BAT, xn_BAT
    cdef float eo, en
    cdef int closest_pose_o, closest_pose_n

    xo_Cartesian = np.copy(self.universe.configuration().array)
    xo_BAT = self._BAT_util.BAT(xo_Cartesian, extended=self.extended)
    eo = self.universe.energy()
    if self.extended:
      closest_pose_o = self._closest_pose_Cartesian(\
        xo_Cartesian[self.molecule.heavy_atoms,:])
    else:
      closest_pose_o = self._closest_pose_BAT(xo_BAT[self._BAT_to_perturb])

    cdef int t, dart_towards
    for t in range(ntrials):
      # Choose a pose to dart towards
      dart_towards = closest_pose_o
      while dart_towards==closest_pose_o:
        dart_towards = np.random.choice(len(self.weights), p=self.weights)
      # Generate a trial move
      xn_BAT = np.copy(xo_BAT)
      xn_BAT[self._BAT_to_perturb] = xo_BAT[self._BAT_to_perturb] + self.darts[closest_pose_o][dart_towards]
      
      # Check that the trial move is closest to dart_towards
      if self.extended:
        xn_Cartesian = self._BAT_util.Cartesian(xn_BAT)
        closest_pose_n = self._closest_pose_Cartesian(\
          xn_Cartesian[self.molecule.heavy_atoms,:])
        if (closest_pose_n!=dart_towards):
          continue
      else:
        closest_pose_n = self._closest_pose_BAT(xn_BAT[self._BAT_to_perturb])
        if (closest_pose_n!=dart_towards):
          continue
        xn_Cartesian = self._BAT_util.Cartesian(xn_BAT)

      # Determine energy of new state
      self.universe.setConfiguration(Configuration(self.universe, xn_Cartesian))
      en = self.universe.energy()

      # Accept or reject the trial move
      if (closest_pose_n==dart_towards) and (abs(en-eo)<1000) and \
          ((en<eo) or (np.random.random()<np.exp(-(en-eo)/RT))):
        xo_Cartesian = xn_Cartesian
        xo_BAT = xn_BAT
        eo = 1.*en
        closest_pose_o = closest_pose_n
        acc += 1

    self.universe.setConfiguration(Configuration(self.universe, xo_Cartesian))
    return ([np.copy(xo_Cartesian)], [eo], acc, ntrials, 0.0)