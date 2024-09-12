import AlGDock

from MMTK import *
from Cylinder import CylinderForceField
from MMTK.ForceFields.ForceFieldTest import gradientTest, forceConstantTest

import numpy as N

try:
  from Scientific._vector import Vector
except:
  from Scientific.Geometry.VectorModule import Vector

universe = InfiniteUniverse()
universe.atom1 = Atom('C', position=Vector((2, 4.1, 2.1)))

ff1 = CylinderForceField(\
  origin=N.array((2., 2., 1.)), \
  direction=N.array((0,0,1.0)), \
  max_Z=2., max_R=2.)
universe.setForceField(ff1)

e, g = universe.energyAndGradients()
print universe.energyTerms()
print e

gradientTest(universe, delta=0.001)
#forceConstantTest(universe)
