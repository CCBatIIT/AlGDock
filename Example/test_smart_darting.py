# The example is 1of6

execfile('test_python.py')
BPMF = self

import MMTK
from MMTK import Units, Configuration
import numpy as np

R = 8.3144621*Units.J/Units.mol/Units.K

# Use a force field without grids
self._set_universe_evaluator({'MM':True, 'T':300.*MMTK.Units.K, \
                            'delta_t':1.5*MMTK.Units.fs,
                            'a':0.0, 'crossed':False})

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
self._set_universe_evaluator(self._lambda(1.0,'dock'))

from AlGDock.Integrators.SmartDarting.SmartDarting import SmartDartingIntegrator
confs = BPMF._get_confs_to_rescore(site=True, minimize=True)[0]

self = SmartDartingIntegrator(BPMF.universe, BPMF.molecule, extended=True)
print self.set_confs(confs)
print 'Coordinates to perturb:', ' '.join(['%d'%d for d in self._BAT_to_perturb])
print 'Acceptance rate,    extended coordinates: ', self(T=300, ntrials=100)[2]

import AlGDock.IO
IO_dcd = AlGDock.IO.dcd(molecule)
IO_dcd.write('smart_darting.dcd', self.confs)
self._BAT_util.showMolecule(dcdFN='smart_darting.dcd')
os.remove('smart_darting.dcd')