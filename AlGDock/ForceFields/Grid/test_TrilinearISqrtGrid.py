import AlGDock

from MMTK import *
import TrilinearISqrtGrid
from MMTK.ForceFields.ForceFieldTest import gradientTest

universe = InfiniteUniverse()

universe.atom1 = Atom('C', position=Vector(1.02, 0.512, 1.513))
universe.atom1.test_charge = 1.
universe.atom2 = Atom('C', position=Vector(1.553, 1.724, 1.464))
universe.atom2.test_charge = -0.2

ForceField = TrilinearISqrtGrid.TrilinearISqrtGridForceField(
  '../../../Example/grids/LJa.nc',
  1.0,'test_charge')

universe.setForceField(ForceField)

print 'Energy Terms:'
print universe.energyTerms()
e, g = universe.energyAndGradients()
print 'Gradient on Atom 1'
print g[universe.atom1]
print 'Gradient on Atom 2'
print g[universe.atom2]

print 'Gradient Test'
gradientTest(universe)

import time
start_time = time.time()
for t in range(100000):
  e, g = universe.energyAndGradients()
print 'Time to do 100000 energy and gradient evaluations'
print time.time()-start_time
