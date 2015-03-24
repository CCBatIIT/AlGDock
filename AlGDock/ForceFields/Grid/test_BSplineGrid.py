import numpy as np
import pylab as pl
import AlGDock
import numpy as np
from MMTK import *
import BSplineGrid
import matplotlib as plt
from MMTK.ForceFields.ForceFieldTest import gradientTest
from MMTK.ForceFields.ForceFieldTest import forceConstantTest

universe = InfiniteUniverse()

universe.atom1 = Atom('C', position=Vector(1.02, 0.512, 1.513))
universe.atom1.test_charge = 1.

universe.atom2 = Atom('C', position=Vector(1.553, 1.724, 1.464))
universe.atom2.test_charge = -0.2
'''
universe.atom3 = Atom('C', position=Vector(1.0, 1.724, 1.464))
universe.atom3.test_charge = 1
'''
ForceField = BSplineGrid.BSplineGridForceField(
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
print 'Second derivative Test'
forceConstantTest(universe)

#import time
#start_time = time.time()

#for t in range(100000):
#   e, g = universe.energyAndGradients()
#print 'Time to do 100000 energy and gradient evaluations'
#print time.time()-start_time


#print '------------------Graph-------------------'
#step = 10000
#x=np.linspace(1.1,1.25,step)
#y=np.linspace(1.1,1.25,step)
#i=0;
#for i in range(step):
#    universe.atom1.setPosition(Vector(x[i],1.5,1.5))
#    y[i]=universe.energy()
#pl.plot(x, y,'b-',label =u"B spline ")
#pl.legend()
#pl.show()
