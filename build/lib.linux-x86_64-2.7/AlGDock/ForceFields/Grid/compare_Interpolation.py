import AlGDock

from MMTK import *
import Interpolation
from MMTK.ForceFields.ForceFieldTest import gradientTest

universe = InfiniteUniverse()

universe.atom1 = Atom('C', position=Vector(1.1, 0.5, 1.5))
universe.atom1.test_charge = 1.
universe.atom2 = Atom('C', position=Vector(1.553, 1.724, 1.464))
universe.atom2.test_charge = -0.2

param_sets = [\
  {'interpolation_type':'Trilinear', 'inv_power':4, 'energy_thresh':-1.0}]

steps = 100000
import numpy as np

from collections import OrderedDict
Es = OrderedDict()
x=np.linspace(1.35,1.6,steps)

for params in param_sets:
  print
  print params
  print
  
  ForceField = Interpolation.InterpolationForceField(\
    '../../../Example/grids/LJa.nc',
    name='Name',
    interpolation_type=params['interpolation_type'],
    inv_power=params['inv_power'],
    energy_thresh=params['energy_thresh'],
    scaling_property='test_charge')

  print "With C"
  universe.setForceField(ForceField)

  universe.atom1.setPosition(Vector(x[0],0.5,1.5))

#  print 'Energy Terms:'
#  print universe.energyTerms()
#  e, g = universe.energyAndGradients()

  print 'Gradient Test'
  gradientTest(universe)

  import time
  start_time = time.time()
  
  for n in range(steps):
    universe.atom1.setPosition(Vector(x[n],0.5,1.5))
    e, g = universe.energyAndGradients()
  print 'Time to do %d energy and gradient evaluations: %f s'%(\
    steps, time.time()-start_time)

  print "Without C"
  ForceField.use_C = False

  universe.setForceField(ForceField)

  universe.atom1.setPosition(Vector(x[0],0.5,1.5))

#  print 'Energy Terms:'
#  print universe.energyTerms()
#  e, g = universe.energyAndGradients()

  print 'Gradient Test'
  gradientTest(universe)

  import time
  start_time = time.time()
  
  for n in range(steps):
    universe.atom1.setPosition(Vector(x[n],0.5,1.5))
    e, g = universe.energyAndGradients()
  print 'Time to do %d energy and gradient evaluations: %f s'%(\
    steps, time.time()-start_time)


#  import matplotlib.pyplot as plt
#  for key in Es.keys():
#    plt.plot(x,Es[key])
#  plt.legend(Es.keys())
#  plt.savefig('test_Interpolation_results.jpg')

# These are results from the cython module
#  Energy Terms:
#  {'Interpolation': -0.5035438268671999}
#  Gradient on Atom 1
#  [2.4492360000000035, 4.4652519999999996, 2.356524000000002]
#  Gradient on Atom 2
#  [-1.1748288704000007, 0.5334514112000035, -0.6846025087999981]
#  Gradient Test
#  Energy:  -0.503543826867
#  Atom carbon
#  [2.4492360000000035, 4.4652519999999996, 2.356524000000002]
#  [3.106016000001266, 6.027314000000561, 2.9514220000015357]
#  Atom carbon
#  [-1.1748288704000007, 0.5334514112000035, -0.6846025087999981]
#  [-1.174828870399991, 0.53345141119987, -0.6846025088003138]
#  Time to do 50000 energy and gradient evaluations: 0.999300 s
