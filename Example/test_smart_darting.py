# The example is 1of6

execfile('test_python.py')
BPMF = self

import MMTK
from MMTK import Units, Configuration
import numpy as np

R = 8.3144621*Units.J/Units.mol/Units.K

#  # Use a force field without grids
#  self._set_universe_evaluator({'MM':True, 'T':300.*MMTK.Units.K, \
#                              'delta_t':1.5*MMTK.Units.fs,
#                              'a':0.0, 'crossed':False})

#  from AlGDock.Integrators.SmartDarting.SmartDarting import SmartDartingIntegrator
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

self = BPMF
self._set_universe_evaluator(self._lambda(0.25,'dock'))

from AlGDock.Integrators.SmartDarting.SmartDarting import SmartDartingIntegrator
confs = BPMF._get_confs_to_rescore(site=True, minimize=True)[0]

self = SmartDartingIntegrator(BPMF.universe, BPMF.molecule, extended=True)
print self.set_confs(confs)
print 'Coordinates to perturb:', ' '.join(['%d'%d for d in self._BAT_to_perturb])
print 'Acceptance rate,    extended coordinates: ', self(T=300, ntrials=100)[2]

confs = []
BPMF.sampler['dock'](steps=50, T=300)

# TODO: This 'dart' can rotate the molecule too much!
# TODO: Make sure there are adequate copies in SmartDarting.py and ExternalMC.py

xo_Cartesian = np.copy(self.universe.configuration().array)
xo_BAT = self._BAT_util.BAT(xo_Cartesian, extended=self.extended)
eo = self.universe.energy()
if self.extended:
  closest_pose_o = self._closest_pose_Cartesian(\
    xo_Cartesian[self.molecule.heavy_atoms,:])
else:
  closest_pose_o = self._closest_pose_BAT(xo_BAT[self._BAT_to_perturb])

# Choose a pose to dart towards
dart_towards = closest_pose_o
while dart_towards==closest_pose_o:
  dart_towards = np.random.choice(len(self.weights), p=self.weights)
# Generate a trial move
xn_BAT = np.copy(xo_BAT)
xn_BAT[self._BAT_to_perturb] = xo_BAT[self._BAT_to_perturb] + self.darts[closest_pose_o][dart_towards]
xn_Cartesian = self._BAT_util.Cartesian(xn_BAT) # Also sets the universe
en = self.universe.energy()

#  import AlGDock.IO
#  IO_dcd = AlGDock.IO.dcd(self.molecule)
#  IO_dcd.write('smart_darting.dcd', \
#    self.confs + [xo_Cartesian, xn_Cartesian])
#  self._BAT_util.showMolecule(dcdFN='smart_darting.dcd')
#  os.remove('smart_darting.dcd')

for conf in confs:
  self.universe.setConfiguration(Configuration(self.universe,conf))
  print self.universe.energyTerms()
