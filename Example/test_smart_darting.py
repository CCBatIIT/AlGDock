# The example is 1of6

execfile('test_python.py')

import MMTK
from MMTK import Units, Configuration
import numpy as np

R = 8.3144621*Units.J/Units.mol/Units.K

# Use a force field without grids
self._set_universe_evaluator({'MM':True, 'T':300.*MMTK.Units.K, \
                            'delta_t':1.5*MMTK.Units.fs,
                            'a':0.0, 'crossed':False})
BPMF = self

# python version
# from AlGDock.Integrators.SmartDarting.SmartDarting import SmartDartingIntegrator

#  cython version
#  from SmartDarting import SmartDartingIntegrator
#  confs = BPMF._get_confs_to_rescore(site=False, minimize=True)[0]
#
#  self = SmartDartingIntegrator(BPMF.universe, BPMF.molecule, extended=False)
#  print self.set_confs(confs)
#  print 'Coordinates to perturb:', ' '.join(['%d'%d for d in self._BAT_to_perturb])
#  print 'Acceptance rate, no extended coordinates: ', self(T=300, ntrials=100)[2]
#
#  self = SmartDartingIntegrator(BPMF.universe, BPMF.molecule, extended=True)
#  print self.set_confs(confs)
#  print 'Coordinates to perturb:', ' '.join(['%d'%d for d in self._BAT_to_perturb])
#  print 'Acceptance rate,    extended coordinates: ', self(T=300, ntrials=100)[2]


RT = R*300.
ntrials = self.getOption('ntrials')

acc = 0.
energies = []
closest_poses = []

xo_Cartesian = np.copy(self.universe.configuration().array)
xo_BAT = self._BAT_util.BAT(xo_Cartesian, extended=self.extended)
eo = self.universe.energy()
print eo
if self.extended:
  closest_pose_o = self._closest_pose_Cartesian(\
    xo_Cartesian[self.molecule.heavy_atoms,:])
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
  print en

  # Accept or reject the trial move
  if (abs(en-eo)<1000) and \
      ((en<eo) or (np.random.random()<np.exp(-(en-eo)/RT))):
    xo_Cartesian = xn_Cartesian
    xo_BAT = xn_BAT
    eo = 1.*en
    closest_pose_o = closest_pose_n
    acc += 1
    print 'Accepted'

self.universe.setConfiguration(Configuration(self.universe, xo_Cartesian))
