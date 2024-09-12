# This module implements a Hamiltonian Monte Carlo "integrator"
# with dual Hamiltonians: one for guidance and one for sampling.
# It is a stripped-down version of the velocity verlet integrator
# with a Metropolis acceptance criterion.
# It requires the option 'T' in addition to velocity verlet options.

from MMTK import Dynamics, Environment, Features, Trajectory, Units
import MMTK_dynamics
from MMTK.ParticleProperties import Configuration

from Scientific import N
import numpy as np

R = 8.3144621*Units.J/Units.mol/Units.K

#
# Hamiltonian Monte Carlo integrator
#
class HamiltonianMonteCarloIntegrator(Dynamics.Integrator):

    def __init__(self, universe, sampling_universe, **options):
        Dynamics.Integrator.__init__(self, universe, options)
        self.sampling_universe = sampling_universe
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

        self.universe.initializeVelocitiesToTemperature(self.getOption('T'))
        
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
                'Hamiltonian Monte Carlo step')

        # Variables for velocity assignment
        m3 = np.repeat(np.expand_dims(masses.array,1),3,axis=1)
        sigma_MB = np.sqrt((self.getOption('T')*Units.k_B)/m3)
        natoms = self.universe.numberOfAtoms()

        xs = []
        energies = []

        # Store initial configuration and potential energy
        xo = np.copy(self.universe.configuration().array)
        self.sampling_universe.configuration().array[-natoms:,:] = xo
        pe_o = self.sampling_universe.energy() # <- Using sampling Hamiltonian
                                            # rather than guidance Hamiltonian

        acc = 0
        for t in range(ntrials):
          # Initialize the velocity
          v = self.universe.velocities()
          v.array = np.multiply(sigma_MB,np.random.randn(natoms,3))
    
          # Store total energy
          eo = pe_o + 0.5*np.sum(np.multiply(m3,np.square(v.array)))

          # Run the velocity verlet integrator
          self.run(MMTK_dynamics.integrateVV,
            (self.universe,
             self.universe.configuration().array,
             self.universe.velocities().array) + late_args)

          # Decide whether to accept the move
          self.sampling_universe.configuration().array[-natoms:,:] = \
            self.universe.configuration().array
          pe_n = self.sampling_universe.energy() # <- Using sampling Hamiltonian
                                              # rather than guidance Hamiltonian
          en = pe_n + 0.5*np.sum(np.multiply(m3,np.square(v.array)))

          # TODO: Confirm acceptance criterion
          if ((en<eo) or (np.random.random()<N.exp(-(en-eo)/RT))) and \
             ((abs(pe_o-pe_n)/RT<250.) or (abs(eo-en)/RT<250.)):
            xo = np.copy(self.universe.configuration().array)
            pe_o = pe_n
            acc += 1
            if normalize:
              self.universe.normalizePosition()
          else:
            self.universe.setConfiguration(Configuration(self.universe,xo))
          
          xs.append(np.copy(self.universe.configuration().array))
          energies.append(pe_o)
  
        return (xs, energies, acc, ntrials, delta_t)
