# This module wraps the velocity verlet integrator,
# requiring that the final energy is a number.
# It is NOT recommended for AlGDock because it can easily attain high energies.

from MMTK import Dynamics, Environment, Features, Trajectory, Units
import MMTK_dynamics
from MMTK.ParticleProperties import Configuration

from Scientific import N
import numpy as np

import random

R = 8.3144621*Units.J/Units.mol/Units.K

#
# Velocity Verlet integrator
#
class VelocityVerletIntegrator(Dynamics.Integrator):

    def __init__(self, universe, **options):
        Dynamics.Integrator.__init__(self, universe, options)
        # Supported features: none for the moment, to keep it simple
        self.features = []

    def __call__(self, **options):
        # Process the keyword arguments
        self.setCallOptions(options)
        # Check if the universe has features not supported by the integrator
        Features.checkFeatures(self, self.universe)
      
        RT = R*self.getOption('T')
        delta_t = self.getOption('delta_t')
        
        if 'steps_per_trial' in self.call_options.keys():
          steps_per_trial = self.getOption('steps_per_trial')
          ntrials = self.getOption('steps')/steps_per_trial
        else:
          steps_per_trial = self.getOption('steps')
          ntrials = 1
  
        if 'normalize' in self.call_options.keys():
          normalize = self.getOption('normalize')
        else:
          normalize = False          

        # Seed the random number generator
        if 'random_seed' in self.call_options.keys():
          np.random.seed(self.getOption('random_seed'))
      
        # Get the universe variables needed by the integrator
        masses = self.universe.masses()
        fixed = self.universe.getAtomBooleanArray('fixed')
        nt = self.getOption('threads')
        comm = self.getOption('mpi_communicator')
        evaluator = self.universe.energyEvaluator(threads=nt,
                                                  mpi_communicator=comm)
        evaluator = evaluator.CEvaluator()

        late_args = (
                masses.array, fixed.array, evaluator,
                N.zeros((0, 2), N.Int), N.zeros((0, ), N.Float),
                N.zeros((1,), N.Int),
                N.zeros((0,), N.Float), N.zeros((2,), N.Float),
                N.zeros((0,), N.Float), N.zeros((1,), N.Float),
                delta_t, self.getOption('first_step'),
                steps_per_trial, self.getActions(),
                'Velocity Verlet step')
      
        xs = []
        energies = []

        # Store initial configuration and potential energy
        xo = np.copy(self.universe.configuration().array)
        pe_o = self.universe.energy()
      
        acc = 0
        for t in range(ntrials):
          # Initialize the velocity
          self.universe.initializeVelocitiesToTemperature(self.getOption('T'))
          # Store total energy
          eo = pe_o + self.universe.kineticEnergy()

          # Run the velocity verlet integrator
          self.run(MMTK_dynamics.integrateVV,
            (self.universe,
             self.universe.configuration().array,
             self.universe.velocities().array) + late_args)

          # Decide whether to accept the move
          pe_n = self.universe.energy()
          en = pe_n + self.universe.kineticEnergy()
          if not np.isnan(en):
            xo = np.copy(self.universe.configuration().array)
            pe_o = pe_n
            acc += 1
            if normalize:
              self.universe.normalizePosition()
          else:
            self.universe.setConfiguration(Configuration(self.universe,xo))
          
          xs.append(np.copy(self.universe.configuration().array))
          energies.append(pe_o)
  
        return (xs, energies, float(acc)/float(ntrials), delta_t)