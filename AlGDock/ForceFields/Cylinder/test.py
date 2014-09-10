from MMTK import *
from Cylinder import CylinderForceField
from MMTK.ForceFields.ForceFieldTest import gradientTest, forceConstantTest

import numpy as N
from Scientific_vector import Vector

universe = InfiniteUniverse()
universe.atom1 = Atom('C', position=Vector((-0.1, -1., 2.)))

ff1 = CylinderForceField(origin=N.array((0., 2., 2.)), direction=N.array((1.0,0,0)), max_X=2., max_R=2.)
universe.setForceField(ff1)

e, g = universe.energyAndGradients()
print universe.energyTerms()
print e

gradientTest(universe, delta=0.001)
#forceConstantTest(universe)
