import os
import cPickle as pickle
import gzip
import copy

from AlGDock.IO import load_pkl_gz
from AlGDock.IO import write_pkl_gz
from AlGDock.logger import NullDevice

import sys
import time
import numpy as np

from collections import OrderedDict

from AlGDock import dictionary_tools
from AlGDock import path_tools

import MMTK
import MMTK.Units
from MMTK.ParticleProperties import Configuration
from MMTK.ForceFields import ForceField

import Scientific
try:
  from Scientific._vector import Vector
except:
  from Scientific.Geometry.VectorModule import Vector

import pymbar.timeseries

import multiprocessing
from multiprocessing import Process

try:
  import requests  # for downloading additional files
except:
  print '  no requests module for downloading additional files'


class LigandPreparation():
  def __init__(self, ):
    pass

  def _ramp_T(self, T_START, T_LOW=20., normalize=False):
    self.log.recordStart('T_ramp')

    # First minimize the energy
    from MMTK.Minimization import SteepestDescentMinimizer  # @UnresolvedImport
    minimizer = SteepestDescentMinimizer(self.top.universe)

    original_stderr = sys.stderr
    sys.stderr = NullDevice()  # Suppresses warnings for minimization

    x_o = np.copy(self.top.universe.configuration().array)
    e_o = self.top.universe.energy()
    for rep in range(5000):
      minimizer(steps=10)
      x_n = np.copy(self.top.universe.configuration().array)
      e_n = self.top.universe.energy()
      diff = abs(e_o - e_n)
      if np.isnan(e_n) or diff < 0.05 or diff > 1000.:
        self.top.universe.setConfiguration(
          Configuration(self.top.universe, x_o))
        break
      else:
        x_o = x_n
        e_o = e_n

    sys.stderr = original_stderr
    self.log.tee("  minimized to %.3g kcal/mol over %d steps" % (e_o, 10 *
                                                                 (rep + 1)))

    # Then ramp the energy to the starting temperature
    from AlGDock.Integrators.HamiltonianMonteCarlo.HamiltonianMonteCarlo \
      import HamiltonianMonteCarloIntegrator
    sampler = HamiltonianMonteCarloIntegrator(self.top.universe)

    e_o = self.top.universe.energy()
    T_LOW = 20.
    T_SERIES = T_LOW * (T_START / T_LOW)**(np.arange(30) / 29.)
    for T in T_SERIES:
      delta_t = 2.0 * MMTK.Units.fs
      steps_per_trial = 10
      attempts_left = 10
      while attempts_left > 0:
        random_seed = int(T*10000) + attempts_left + \
          int(self.top.universe.configuration().array[0][0]*10000)
        if self.args.random_seed == 0:
          random_seed += int(time.time() * 1000)
        random_seed = random_seed % 32767
        (xs, energies, acc, ntrials, delta_t) = \
          sampler(steps = 2500, steps_per_trial = 10, T=T,\
                  delta_t=delta_t, random_seed=random_seed)
        attempts_left -= 1
        acc_rate = float(acc) / ntrials
        if acc_rate < 0.4:
          delta_t -= 0.25 * MMTK.Units.fs
        else:
          attempts_left = 0
        if delta_t < 0.1 * MMTK.Units.fs:
          delta_t = 0.1 * MMTK.Units.fs
          steps_per_trial = max(int(steps_per_trial / 2), 1)
      fmt = "  T = %d, delta_t = %.3f fs, steps_per_trial = %d, acc_rate = %.3f"
      if acc_rate < 0.01:
        print self.top.universe.energyTerms()
      self.log.tee(fmt % (T, delta_t * 1000, steps_per_trial, acc_rate))
    if normalize:
      self.top.universe.normalizePosition()
    e_f = self.top.universe.energy()

    self.log.tee("  ramped temperature from %d to %d K in %s, "%(\
      T_LOW, T_START, HMStime(self.log.timeSince('T_ramp'))) + \
      "changing energy to %.3g kcal/mol"%(e_f))
