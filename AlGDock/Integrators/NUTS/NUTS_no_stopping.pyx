# This module implements the No-U-Turn Sampler (NUTS)
#
# Written by David Minh
#

import numpy as np
cimport numpy as np
import cython

cimport MMTK_trajectory_generator
from MMTK import Units
from MMTK import Features

import MMTK_trajectory
import MMTK_forcefield

cdef extern from "stdlib.h":

    ctypedef long size_t
    cdef void *malloc(size_t size)
    cdef void free(void *ptr)

from MMTK.ParticleProperties import Configuration, ParticleVector

include "MMTK/python.pxi"
include "MMTK/numeric.pxi"
include "MMTK/core.pxi"
include "MMTK/universe.pxi"
include "MMTK/trajectory.pxi"
include "MMTK/forcefield.pxi"

R = 8.3144621*Units.J/Units.mol/Units.K


#
# NUTS integrator
#
cdef class NUTSIntegrator(MMTK_trajectory_generator.EnergyBasedTrajectoryGenerator):

  """
  NUTS integrator
  The integrator is fully thread-safe.
  The integration is started by calling the integrator object.
  All the keyword options (see documnentation of __init__) can be
  specified either when creating the integrator or when calling it.
  The following data categories and variables are available for
  output:
   - category "time": time
   - category "configuration": configuration and box size (for
     periodic universes)
   - category "velocities": atomic velocities
   - category "gradients": energy gradients for each atom
   - category "energy": potential and kinetic energy, plus
     extended-system energy terms if a thermostat and/or barostat
     are used
  """

  cdef double RT
  cdef np.ndarray x, v, g, m
  cdef energy_data energy

  def __init__(self, universe, **options):
    """
    @param universe: the universe on which the integrator acts
    @type universe: L{MMTK.Universe}
    @keyword steps: the number of integration steps (default is 100)
    @type steps: C{int}
    @keyword delta_t: the time step (default is 1 fs)
    @type delta_t: C{float}
    @keyword actions: a list of actions to be executed periodically
                      (default is none)
    @type actions: C{list}
    @keyword threads: the number of threads to use in energy evaluation
                      (default set by MMTK_ENERGY_THREADS)
    @type threads: C{int}
    @keyword background: if True, the integration is executed as a
                         separate thread (default: False)
    @type background: C{bool}
    """
    MMTK_trajectory_generator.EnergyBasedTrajectoryGenerator.__init__(
        self, universe, options, "NUTS integrator")
    # Supported features: none for the moment, to keep it simple
    self.features = []

  default_options = {'first_step': 0, 'steps': 100, 'delta_t': 1.*Units.fs,
                     'background': False, 'threads': None,
                     'actions': []}

  available_data = ['configuration', 'velocities', 'gradients',
                    'energy', 'time']

  restart_data = ['configuration', 'velocities', 'energy']

  def __call__(self, **options):
    self.setCallOptions(options)
    try:
        self.actions = self.getOption('actions')
    except ValueError:
        self.actions = []
    try:
        if self.getOption('background'):
            import MMTK_state_accessor
            self.state_accessor = MMTK_state_accessor.StateAccessor()
            self.actions.append(self.state_accessor)
    except ValueError:
        pass
    Features.checkFeatures(self, self.universe)
    if self.tvars != NULL:
        free(self.tvars)
        self.tvars = NULL
    configuration = self.universe.configuration()
    self.conf_array = configuration.array
    self.declareTrajectoryVariable_array(self.conf_array,
                                         "configuration",
                                         "Configuration:\n",
                                         length_unit_name,
                                         PyTrajectory_Configuration)
    self.universe_spec = <PyUniverseSpecObject *>self.universe._spec
    if self.universe_spec.geometry_data_length > 0:
        self.declareTrajectoryVariable_box(
            self.universe_spec.geometry_data,
            self.universe_spec.geometry_data_length)
    masses = self.universe.masses()
    self.declareTrajectoryVariable_array(masses.array,
                                         "masses",
                                         "Masses:\n",
                                         mass_unit_name,
                                         PyTrajectory_Internal)
    self.natoms = self.universe.numberOfAtoms()
    self.df = self.universe.degreesOfFreedom()
    self.declareTrajectoryVariable_int(&self.df,
                                       "degrees_of_freedom",
                                       "Degrees of freedom: %d\n",
                                       "", PyTrajectory_Internal)
    if self.getOption('background'):
        from MMTK import ThreadManager
        return ThreadManager.TrajectoryGeneratorThread(
            self.universe, self.start_py, (), self.state_accessor)
    else:
        return self.start()

  # Cython compiler directives set for efficiency:
  # - No bound checks on index operations
  # - No support for negative indices
  # - Division uses C semantics
  @cython.boundscheck(False)
  @cython.wraparound(False)
  @cython.cdivision(True)
  cdef start(self):

    cdef double time, delta_t, ke
    cdef int natoms, nsteps, elapsed_steps
    cdef Py_ssize_t i_atm, i_dim

    cdef double joint, logu, e_m
    cdef np.ndarray[double, ndim=2] xminus, xplus, vminus, vplus, x_m
    cdef int j, n, steps_m
        
    # For dual averaging
    cdef double alpha, delta, gamma, kappa, mu, delta_t_bar, Hbar, eta
    cdef int nalpha, t0, m

    # Initialize the velocity
    self.universe.initializeVelocitiesToTemperature(self.getOption('T'))

    # Gather state variables and parameters
    configuration = self.universe.configuration()
    velocities = self.universe.velocities()
    gradients = ParticleVector(self.universe)
    masses = self.universe.masses()
    delta_t = self.getOption('delta_t')
    nsteps = self.getOption('steps')
    natoms = self.universe.numberOfAtoms()

    self.RT = R*self.getOption('T')

    if 'normalize' in self.call_options.keys():
      normalize = self.getOption('normalize')
    else:
      normalize = False
      
    if 'delta' in self.call_options.keys():
      delta = self.getOption('delta')
    else:
      delta = 0.6

    if 'adapt' in self.call_options.keys():
      adapt = self.getOption('adapt')
    else:
      adapt = False

    # Parameters for the dual averaging algorithm.
    if adapt:
      gamma = 0.05
      t0 = 10
      kappa = 0.75
      mu = np.log(10*delta_t)
      # Initialize dual averaging algorithm.
      delta_t_bar = 1.
      Hbar = 0.
    else:
      delta_t_bar = delta_t

    # Seed the random number generator
    if 'seed' in self.call_options.keys():
      np.random.seed(self.getOption('seed'))
    else:
      np.random.seed()

    # For efficiency, the Cython code works at the array
    # level rather than at the ParticleProperty level.
    self.x = configuration.array
    self.v = velocities.array
    self.g = gradients.array
    self.m = np.repeat(np.expand_dims(masses.array,1),3,axis=1)

    # Weight matrix for velocity assignment
    sigma_MB = np.sqrt((self.getOption('T')*Units.k_B)/self.m)
    
    # Ask for energy gradients to be calculated and stored in
    # the array g. Force constants are not requested.
    self.energy.gradients = <void *>self.g
    self.energy.gradient_fn = NULL
    self.energy.force_constants = NULL
    self.energy.fc_fn = NULL
    
    # Declare the variables accessible to trajectory actions.
    self.declareTrajectoryVariable_double(
        &time, "time", "Time: %lf\n", time_unit_name, PyTrajectory_Time)
    self.declareTrajectoryVariable_array(
        self.v, "velocities", "Velocities:\n", velocity_unit_name,
        PyTrajectory_Velocities)
    self.declareTrajectoryVariable_array(
        self.g, "gradients", "Energy gradients:\n", energy_gradient_unit_name,
        PyTrajectory_Gradients)
    self.declareTrajectoryVariable_double(
        &self.energy.energy,"potential_energy", "Potential energy: %lf\n",
        energy_unit_name, PyTrajectory_Energy)
    self.declareTrajectoryVariable_double(
        &ke, "kinetic_energy", "Kinetic energy: %lf\n",
        energy_unit_name, PyTrajectory_Energy)
    self.initializeTrajectoryActions()

    # Acquire the write lock of the universe. This is necessary to
    # make sure that the integrator's modifications to positions
    # and velocities are synchronized with other threads that
    # attempt to use or modify these same values.
    #
    # Note that the write lock will be released temporarily
    # for trajectory actions. It will also be converted to
    # a read lock temporarily for energy evaluation. This
    # is taken care of automatically by the respective methods
    # of class EnergyBasedTrajectoryGenerator.
    self.acquireWriteLock()

    # Initialize the NUTS algorithm
    self.calculateEnergies(self.x, &self.energy, 0)

    xs = [np.copy(self.x)]
    energies = [1.*self.energy.energy]
    
    # Initialize the next sample.
    # If all else fails, the next sample is the previous sample.
    x_m = np.copy(self.x)
    e_m = 1.*self.energy.energy

    # Main integration loop
    elapsed_steps = 0
    m = 1
    while elapsed_steps < nsteps:
      # Resample velocities
      self.v = np.multiply(sigma_MB,np.random.randn(natoms,3))
      ke = 0.5*np.sum(np.multiply(self.m,np.multiply(self.v,self.v)))

      # Joint log-probabiity of positions and velocities
      joint = -(e_m + ke)/self.RT

      # Resample u ~ uniform([0, np.exp(joint)]).
      # Equivalent to (log(u) - joint) ~ exponential(1).
      logu = joint - np.random.exponential(1)
      
      # Initialize tree.
      xminus = np.copy(x_m)
      xplus = np.copy(x_m)
      vminus = np.copy(self.v)
      vplus = np.copy(self.v)

      # Initial height j = 0.
      j = 0
      # Initially the only valid point is the initial point.
      n = 1
      # Initially there have been no integrator steps
      steps_m = 0

      # Main loop
      s = True
      while (s==1):
        # Double the size of the tree
        if np.random.rand()<0.5:
          # Backwards
          self.x = xminus
          self.v = vminus
          self.calculateEnergies(self.x, &self.energy, 0)
          (xminus, vminus, NA0, NA1, 
           xprime, eprime, nprime, sprime, steps_m,
           alpha, nalpha) = \
            self.build_tree(logu, j, -1*delta_t, steps_m, joint)
        else:
          # Forward
          self.x = xplus
          self.v = vplus
          self.calculateEnergies(self.x, &self.energy, 0)
          (NA0, NA1, xplus, vplus, 
           xprime, eprime, nprime, sprime, steps_m,
           alpha, nalpha) = \
            self.build_tree(logu, j, delta_t, steps_m, joint)
        # Use Metropolis-Hastings to decide whether or not to move to a
        # point from the half-tree we just generated
        if (sprime and (np.random.rand() < float(nprime)/n)):
          x_m = xprime
          e_m = eprime
        # Update number of valid points we've seen.
        n += nprime
        # Increment depth
        j += 1
        # Decide if it's time to stop
        s = sprime and \
            ((elapsed_steps + steps_m*2) < nsteps)

      # Keep track of acceptance statistics
      eta = 1./float(m+t0)
      Hbar = (1-eta)*Hbar + eta*(delta-alpha/nalpha)

      # Adapt the time step
      if adapt:
        delta_t = np.exp(mu - np.sqrt(m)/gamma*Hbar)
        eta = m**-kappa
        delta_t_bar = np.exp((1-eta)*np.log(delta_t_bar) + eta*np.log(delta_t))

      xs.append(np.copy(x_m))
      energies.append(1.*e_m)
            
      time += steps_m*delta_t
      self.trajectoryActions(m)
      
      m += 1
      elapsed_steps += steps_m
    
    self.universe.setConfiguration(Configuration(self.universe, x_m), block=False)
    if normalize:
      self.universe.normalizePosition()

    # Release the write lock.
    self.releaseWriteLock()

    # Finalize all trajectory actions (close files etc.)
    self.finalizeTrajectoryActions(nsteps)
    
    return (xs, energies, Hbar, delta_t_bar)

  # The main recursion
  def build_tree(NUTSIntegrator self, double logu, int j, double delta_t, int steps, double joint_o):
    cdef double ke, joint, e_o, eprime, eprime2
    cdef np.ndarray[double, ndim=2] xminus, xplus, vminus, vplus, xprime, xprime2, NA0, NA1
    
    cdef double alphaprime, alphaprime2
    cdef int nalphaprime, nalphaprime2

    if (j==0):
      # Base case: Take a single leapfrog step
      # First half-step
      e_o = 1.*self.energy.energy
      self.v += -0.5*delta_t*np.divide(self.g,self.m)
      self.x += delta_t*self.v      
      # Mid-step energy calculation
      self.foldCoordinatesIntoBox()
      self.calculateEnergies(self.x, &self.energy, 1)
      # Second half-step
      self.v += -0.5*delta_t*np.divide(self.g,self.m)
      ke = 0.5*np.sum(np.multiply(self.m,np.multiply(self.v,self.v)))
      steps += 1

      eprime = 1.*self.energy.energy
      joint = -(eprime+ke)/self.RT
      # Is the new point in the slice?
      nprime = (logu < joint)
      # Is the simulation wildly inaccurate
      sprime = (np.abs(joint_o - joint) < 50.) and (np.abs((e_o - eprime)/self.RT) < 50.)
      # Compute the acceptance probability
      alphaprime = min(1, np.exp(joint - joint_o)) if sprime else 0.
      # Set the return values---minus=plus for all things here, since the
      # "tree" is of depth 0.
      return (np.copy(self.x), np.copy(self.v),
              np.copy(self.x), np.copy(self.v), 
              np.copy(self.x),
              eprime, nprime, sprime, steps,
              alphaprime, 1)
    else:
      (xminus, vminus, xplus, vplus, 
       xprime, eprime, nprime, sprime, steps,
       alphaprime, nalphaprime) = \
        self.build_tree(logu, j-1, delta_t, steps, joint_o)
      # No need to keep going if the stopping criteria were met in the first
      # subtree.
      if sprime:
        if (delta_t < 0):
          (xminus, vminus, NA0, NA1, 
           xprime2, eprime2, nprime2, sprime2, steps,
           alphaprime2, nalphaprime2) = \
              self.build_tree(logu, j-1, delta_t, steps, joint_o)
        else:
          (NA0, NA1, xplus, vplus, 
           xprime2, eprime2, nprime2, sprime2, steps,
           alphaprime2, nalphaprime2) = \
              self.build_tree(logu, j-1, delta_t, steps, joint_o)
        # Choose which subtree to propagate a sample up from.
        if ((nprime + nprime2) > 0) and \
           (np.random.rand() < float(nprime2) / (nprime + nprime2)):
            xprime = xprime2
            eprime = eprime2
        # Update the number of valid points.
        nprime = nprime + nprime2
        # Update the stopping criterion.
        sprime = sprime and sprime2
        # Update the acceptance probability statistics
        alphaprime = alphaprime + alphaprime2
        nalphaprime = nalphaprime + nalphaprime2
      return (xminus, vminus, xplus, vplus, xprime,
              eprime, nprime, sprime, steps,
              alphaprime, nalphaprime)