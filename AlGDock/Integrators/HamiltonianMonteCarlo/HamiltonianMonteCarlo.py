# This module implements a Hamiltonian Monte Carlo "integrator"
# It is a stripped-down version of the velocity verlet integrator
# with a Metropolis acceptance criterion.
# It requires the option 'T' in addition to velocity verlet options.

import numpy as np

try:
  from MMTK import Dynamics, Environment, Features, Units
  import MMTK_dynamics
  from MMTK.ParticleProperties import Configuration
  from Scientific import N
except ImportError:
  MMTK = None

try:
  import openmm
  import openmm.unit as unit
  from openmm.app import AmberPrmtopFile, AmberInpcrdFile, Simulation, NoCutoff
  from openmm import *
  from openmmtools.integrators import VelocityVerletIntegrator
except ImportError:
  OpenMM = None

R = 8.3144621*Units.J/Units.mol/Units.K

#
# Hamiltonian Monte Carlo integrator
#
class HamiltonianMonteCarloIntegrator(Dynamics.Integrator):

    def __init__(self, universe, **options):
        Dynamics.Integrator.__init__(self, universe, options)
        # Supported features: none for the moment, to keep it simple
        self.features = []

    def __call__(self, **options):
        # Process the keyword arguments
        self.setCallOptions(options)
        # Check if the universe has features not supported by the integrator
        Features.checkFeatures(self, self.universe)
        print('options:', options)

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
        pe_o = self.universe.energy()

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
          pe_n = self.universe.energy()
          en = pe_n + 0.5*np.sum(np.multiply(m3,np.square(v.array)))
          
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


class HamiltonianMonteCarloIntegratorUsingOpenMM:
    def __init__(self, molecule, top, OMM_system):
        self.options = {'first_step': 0, 'steps': 100, 'delta_t': 1. * Units.fs,
                           'background': False, 'threads': None,
                           'mpi_communicator': None, 'actions': []}
        self.top = top
        self.molecule = molecule
        self.OMM_system = OMM_system
        self.main()
    def main(self):
        RT = R *  self.options['T']
        delta_t = self.options['delta_t']

        if 'steps_per_trial' in self.options.keys():
            steps_per_trial = self.options['steps_per_trial']
            ntrials = self.options['steps'] / steps_per_trial
        else:
            steps_per_trial = self.options['steps']
            ntrials = 1

        if 'normalize' in self.options.keys():
            normalize = self.options['normalize']
        else:
            normalize = False

            # Seed the random number generator
        if 'random_seed' in self.options.keys():
            np.random.seed(self.options['random_seed'])

        context = self.top.OMM_simulaiton.context
        context.setVelocitiesToTemperature(self.options['T'] * unit.kelvin)
        #self.universe.initializeVelocitiesToTemperature(self.options('T'))

        # Get the universe variables needed by the integrator
        masses = []
        for atom in self.molecule.atoms():
            mass = self.OMM_system.getParticleMass(atom.index)  # Get mass from system
            masses.append(mass.value_in_unit(unit.dalton))  # Convert to Daltons
        masses = np.array(masses)

        # fixed = self.universe.getAtomBooleanArray('fixed')
        # What is fixed? It's [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

        nt = self.options['threads']
        comm = self.options['mpi_communicator']
        # evaluator = self.universe.energyEvaluator(threads=nt,
        #                                           mpi_communicator=comm)
        # evaluator = evaluator.CEvaluator()

        # late_args = (
        #     masses.array, #fixed.array, evaluator,
        #     N.zeros((0, 2), N.Int), N.zeros((0,), N.Float),
        #     N.zeros((1,), N.Int),
        #     N.zeros((0,), N.Float), N.zeros((2,), N.Float),
        #     N.zeros((0,), N.Float), N.zeros((1,), N.Float),
        #     delta_t, options('first_step'),
        #     steps_per_trial, getActions(),
        #     'Hamiltonian Monte Carlo step')

        # Variables for velocity assignment
        m3 = np.repeat(np.expand_dims(masses.array, 1), 3, axis=1)
        k_B = unit.BOLTZMANN_CONSTANT_kB * unit.kilojoule_per_mole / unit.kelvin
        sigma_MB = np.sqrt((self.options['T'] * k_B) / m3)
        natoms = self.OMM_system.getNumParticles()

        xs = []
        energies = []

        # Store initial configuration and potential energy
        state = self.top.OMM_simulation.context.getState(getPositions=True, getEnergy=True, getVelocities=True)
        xo = np.array(state.getPositions(asNumpy=True))
        pe_o = state.getPotentialEnergy().value_in_unit(openmm.unit.kilojoule_per_mole)

        acc = 0
        for t in range(ntrials):
            # Initialize the velocity
            v = state.getVelocities()

            v.array = np.multiply(sigma_MB, np.random.randn(natoms, 3))

            # Store total energy
            eo = pe_o + 0.5 * np.sum(np.multiply(m3, np.square(v.array)))

            # Run the velocity verlet integrator
            # self.run(MMTK_dynamics.integrateVV,
            #          (self.universe,
            #           self.universe.configuration().array,
            #           self.universe.velocities().array) + late_args)

            integrator = VelocityVerletIntegrator(delta_t)
            simulation = Simulation(self.molecule, self.OMM_system, integrator)

            # Set initial positions & velocities
            simulation.context.setPositions(xo)
            simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)

            # Decide whether to accept the move
            state = simulation.context.getState(getEnergy=True, getPositions=True)
            pe_n = state.getPotentialEnergy().value_in_unit(openmm.unit.kilojoule_per_mole)
            en = pe_n + 0.5 * np.sum(np.multiply(m3, np.square(v.array)))

            if ((en < eo) or (np.random.random() < N.exp(-(en - eo) / RT))) and \
                    ((abs(pe_o - pe_n) / RT < 250.) or (abs(eo - en) / RT < 250.)):
                positions = np.array(state.getPositions(asNumpy=True))
                xo = np.copy(positions)
                pe_o = pe_n
                acc += 1
                if normalize:
                    state = self.top.OMM_simulation.context.getState(getPositions=True,
                                                                     enforcePeriodicBox=True)
                    positions = state.getPositions(asNumpy=True)
                    simulation.context.setPositions(positions)
            else:
                simulation.context.setPositions(xo)

            xs.append(xo)
            energies.append(pe_o)

        return (xs, energies, acc, ntrials, delta_t)