from MMTK import *
from Sphere import SphereForceField
from MMTK.ForceFields.ForceFieldTest import gradientTest, forceConstantTest

import numpy as N
from Scientific._vector import Vector

universe = InfiniteUniverse()
universe.atom1 = Atom('C', position=Vector((-3, 1., 0.)))

ff1 = SphereForceField(center=N.array((0., 0., 1.)), max_R=2.)
universe.setForceField(ff1)

e, g = universe.energyAndGradients()
print universe.energyTerms()
print e

gradientTest(universe, delta=0.001)
#forceConstantTest(universe)
