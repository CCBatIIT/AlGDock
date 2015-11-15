# This module implements am External Monte Carlo move "integrator"

from MMTK import Configuration, Dynamics, Environment, Features, Trajectory, Units
import MMTK_dynamics
import numpy as np

import random

R = 8.3144621*Units.J/Units.mol/Units.K

def random_rotate():
  """
  Return a random rotation matrix
  """
  u = np.random.uniform(size=3)

  # Random quaternion
  q = np.array([np.sqrt(1-u[0])*np.sin(2*np.pi*u[1]),
               np.sqrt(1-u[0])*np.cos(2*np.pi*u[1]),
               np.sqrt(u[0])*np.sin(2*np.pi*u[2]),
               np.sqrt(u[0])*np.cos(2*np.pi*u[2])])
  
  # Convert the quaternion into a rotation matrix 
  rotMat = np.array([[q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3],
                     2*q[1]*q[2] - 2*q[0]*q[3],
                     2*q[1]*q[3] + 2*q[0]*q[2]],
                    [2*q[1]*q[2] + 2*q[0]*q[3],
                     q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3],
                     2*q[2]*q[3] - 2*q[0]*q[1]],
                    [2*q[1]*q[3] - 2*q[0]*q[2],
                     2*q[2]*q[3] + 2*q[0]*q[1],
                     q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3]]])
  return rotMat

#
# External Monte Carlo move integrator
#
class ExternalMCIntegrator(Dynamics.Integrator):
  def __init__(self, universe, molecule, step_size, **options):
    """
    confs - configurations to dart to
    extended - whether or not to use external coordinates
    """
    Dynamics.Integrator.__init__(self, universe, options)
    # Supported features: none for the moment, to keep it simple
    self.features = []

    self.molecule = molecule
    self.step_size = step_size

  def __call__(self, **options):
    # Process the keyword arguments
    self.setCallOptions(options)
    # Check if the universe has features not supported by the integrator
    Features.checkFeatures(self, self.universe)
  
    RT = R*self.getOption('T')
    ntrials = self.getOption('ntrials')

    acc = 0
    xo = np.copy(self.universe.configuration().array)
    eo = self.universe.energy()
    com = self.universe.centerOfMass().array
    
    for c in range(ntrials):
      step = np.random.randn(3)*self.step_size
      if c%2==0:
        # Random translation and full rotation
        xn = np.dot((xo - com), random_rotate()) + com + step
      else:
        # Random translation
        xn = xo + step
      self.universe.setConfiguration(Configuration(self.universe,xn))
      en = self.universe.energy()
      if ((en<eo) or (np.random.random()<np.exp(-(en-eo)/RT))):
        acc += 1
        xo = xn
        eo = en
        com += step
      else:
        self.universe.setConfiguration(Configuration(self.universe,xo))

    return ([np.copy(xo)], [en], acc, ntrials, 0.0)