# This module allows one to combine torsional dynamics
# and Hamiltonian Monte Carlo integrator

# It requires the option 'T' in addition to velocity verlet options.

from MMTK import Dynamics, Environment, Features, Trajectory, Units
import MMTK_dynamics
from MMTK.ParticleProperties import Configuration

from Scientific import N
import numpy as np

R = 8.3144621*Units.J/Units.mol/Units.K

#
# Mixed HMC integrator
#
class MixedHMCIntegrator(Dynamics.Integrator):
  def __init__(self, universe, TDintegrator, **options):
    Dynamics.Integrator.__init__(self, universe, options)
    # Supported features: none for the moment, to keep it simple
    self.features = []
  
    self.TDintegrator = TDintegrator
  
  def __call__(self, **options):
    """
    options include:
    
    steps, 
    steps_per_trial, 
    T, 
    delta_t, 
    random_seed, 
    normalize=False, 
    adapt=False
    """
  
    # Process the keyword arguments
    self.setCallOptions(options)
    # Check if the universe has features not supported by the integrator
    Features.checkFeatures(self, self.universe)

    # Required arguments
    steps = self.getOption('steps')
    delta_t = self.getOption('delta_t')
    T = self.getOption('T')
    
    RT = R*T
    
    # Arguments with default values
    try:
      steps_per_cycle = self.getOption('steps_per_trial')
      ncycles = steps/steps_per_cycle
    except ValueError:
      steps_per_cycle = steps
      ncycles = 1

    try:
      fraction_TD = self.getOption('fraction_TD')
    except ValueError:
      fraction_TD = 0.5
      
    try:
      TD_steps_per_trial = self.getOption('TD_steps_per_trial')
    except ValueError:
      TD_steps_per_trial = 5
    
    try:
      delta_t_TD = self.getOption('delta_t_TD')
    except ValueError:
      delta_t_TD = 4.0

    # Allocate steps per cycle to TD and MD
    TD_steps_per_cycle = int(fraction_TD*steps_per_cycle)
    MD_steps_per_cycle = steps_per_cycle - TD_steps_per_cycle*TD_steps_per_trial

    if 'random_seed' in self.call_options.keys():
      random_seed = self.getOption('random_seed')
      np.random.seed(random_seed)
    else:
      import time
      random_seed = np.int(time.time()*1E5)
    
    if 'normalize' in self.call_options.keys():
      normalize = self.getOption('normalize')
    else:
      normalize = False          
    
    if self.universe.velocities() is None:
      self.universe.initializeVelocitiesToTemperature(T)
    
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
            MD_steps_per_cycle, self.getActions(),
            'Mixed Hamiltonian Monte Carlo step')

    # Variables for velocity assignment
    m3 = np.repeat(np.expand_dims(masses.array,1),3,axis=1)
    sigma_MB = np.sqrt((self.getOption('T')*Units.k_B)/m3)
    natoms = self.universe.numberOfAtoms()

    xs = []
    energies = []

    TD_acc = 0
    TD_ntrials = 0
    MD_acc = 0
    for t in range(ncycles):
      # Do the torsional dynamics steps
      (TD_xs_c, TD_energies_c, TD_acc_c, TD_ntrials_c, TD_dts) = \
        self.TDintegrator.Call(TD_steps_per_cycle, TD_steps_per_trial, \
          T, 0.0040, (random_seed+t)%32767, 1, 1, 0.5)
      xs.extend(TD_xs_c)
      energies.extend(TD_energies_c)
      TD_acc += TD_acc_c
      TD_ntrials += TD_ntrials_c

      # Store initial configuration and potential energy
      xo = np.copy(self.universe.configuration().array)
      pe_o = self.universe.energy()

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
      pe_n = self.universe.energy()
      en = pe_n + 0.5*np.sum(np.multiply(m3,np.square(v.array)))
      
      if ((en<eo) or (np.random.random()<N.exp(-(en-eo)/RT))) and \
         ((abs(pe_o-pe_n)/RT<250.) or (abs(eo-en)/RT<250.)):
        xo = np.copy(self.universe.configuration().array)
        pe_o = pe_n
        MD_acc += 1
        if normalize:
          self.universe.normalizePosition()
      else:
        self.universe.setConfiguration(Configuration(self.universe,xo))
      
      xs.append(np.copy(self.universe.configuration().array))
      energies.append(pe_o)
      
    # print 'TD: %d/%d, MD: %d/%d'%(TD_acc,TD_ntrials,MD_acc,ncycles)
    return (xs, energies, TD_acc + MD_acc, TD_ntrials + ncycles, delta_t)
