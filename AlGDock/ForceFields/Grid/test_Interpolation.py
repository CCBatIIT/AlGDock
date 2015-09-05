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
  {'interpolation_type':'Trilinear', 'inv_power':None, 'energy_thresh':-1.0},
  {'interpolation_type':'Trilinear', 'inv_power':-2, 'energy_thresh':-1.0},
  {'interpolation_type':'Trilinear', 'inv_power':-3, 'energy_thresh':-1.0},
  {'interpolation_type':'Trilinear', 'inv_power':None, 'energy_thresh':10.0},
  {'interpolation_type':'BSpline', 'inv_power':None, 'energy_thresh':-1.0},
  {'interpolation_type':'BSpline', 'inv_power':-3, 'energy_thresh':-1.0},
  {'interpolation_type':'CatmullRom', 'inv_power':None, 'energy_thresh':-1.0},
  {'interpolation_type':'CatmullRom', 'inv_power':-3, 'energy_thresh':-1.0}]

steps = 50000
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
    interpolation_type=params['interpolation_type'],
    inv_power=params['inv_power'],
    energy_thresh=params['energy_thresh'],
    scaling_property='test_charge')

  universe.setForceField(ForceField)

  universe.atom1.setPosition(Vector(x[0],0.5,1.5))

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
  
  key = params['interpolation_type']
  if params['inv_power'] is not None:
    key += ', x**%d'%params['inv_power']
  if params['energy_thresh']>0:
    key += ', e<%f'%params['energy_thresh']

  Es[key] = np.zeros((steps,1))
  for n in range(steps):
    universe.atom1.setPosition(Vector(x[n],0.5,1.5))
    e, g = universe.energyAndGradients()
    Es[key][n] = e
  print 'Time to do %d energy and gradient evaluations: %f s'%(\
    steps, time.time()-start_time)

import matplotlib.pyplot as plt
for key in Es.keys():
  plt.plot(x,Es[key])
plt.legend(Es.keys())
plt.savefig('test_Interpolation_results.jpg')