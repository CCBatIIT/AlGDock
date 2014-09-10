# cython: profile=False

# TO DO: Find out what causes overflow

# This module implements a Velocity Verlet integrator.
#
# This is a stripped down version, without:
#   trajectory output, trajectory actions, or threads

import_MMTK_universe()
import_MMTK_forcefield()

from MMTK import Features
import numpy as N
cimport numpy as N
import cython

cdef extern from "stdlib.h":

    ctypedef long size_t
    cdef void *malloc(size_t size)
    cdef void free(void *ptr)

from MMTK.ParticleProperties import Configuration, ParticleVector
from MMTK import Units
import MMTK_trajectory
import MMTK_forcefield

include "MMTK/python.pxi"
include "MMTK/numeric.pxi"
include "MMTK/core.pxi"
include "MMTK/universe.pxi"
include "MMTK/forcefield.pxi"

R = 8.3144621*Units.J/Units.mol/Units.K

#
# Base class for trajectory generators
#
cdef class TrajectoryGenerator(object):

    """
    Trajectory generator base class

    This base class implements the common aspects of everything that
    generates trajectories: integrators, minimizers, etc.
    """
    
    cdef readonly universe, options, name
    cdef readonly call_options
    cdef readonly features
    cdef readonly actions
    cdef public state_accessor
    cdef PyUniverseSpecObject *universe_spec
    cdef int natoms, df
    cdef N.ndarray conf_array
    cdef int lock_state

    def __init__(self, universe, options, name):
        self.universe = universe
        self.options = options
        self.name = name
        self.call_options = {}
        self.features = []
        
    def setCallOptions(self, options):
        self.call_options = options

    def getActions(self):
        try:
            steps = self.getOption('steps')
        except ValueError:
            steps = None
        return [a.getSpecificationList(self, steps) for a in self.actions]

    def cleanupActions(self):
        for a in self.actions:
            a.cleanup()

    def getOption(self, option):
        try:
            value = self.call_options[option]
        except KeyError:
            try:
                value = self.options[option]
            except KeyError:
                try:
                    value = self.default_options[option]
                except KeyError:
                    raise ValueError('undefined option: ' + option)
        return value

    def optionString(self, options):
        return sum((o + '=' + `self.getOption(o)` + ', ' for o in options),
                   '')[:-2]

    def __call__(self, **options):
        self.setCallOptions(options)
        Features.checkFeatures(self, self.universe)
        configuration = self.universe.configuration()
        self.conf_array = configuration.array
        self.universe_spec = <PyUniverseSpecObject *>self.universe._spec
        masses = self.universe.masses()
        self.natoms = self.universe.numberOfAtoms()
        self.df = self.universe.degreesOfFreedom()
        return self.start()

    def start_py(self):
        self.start()

    cdef start(self):
        pass

    cdef void initializeTrajectoryActions(self) except *:
        pass

    cdef void finalizeTrajectoryActions(self, int last_step,
                                        int error=False) except *:
        pass
        
    cdef void foldCoordinatesIntoBox(self) nogil:
        self.universe_spec.correction_function(<vector3 *>self.conf_array.data,
                                               self.natoms,
                                               self.universe_spec.geometry_data)

    cdef void acquireReadLock(self) nogil:
        PyUniverseSpec_StateLock(self.universe_spec, 1)
        self.lock_state = 1

    cdef void releaseReadLock(self) nogil:
        PyUniverseSpec_StateLock(self.universe_spec, 2)
        self.lock_state = 0

    cdef void acquireWriteLock(self) nogil:
        PyUniverseSpec_StateLock(self.universe_spec, -1)
        self.lock_state = -1

    cdef void releaseWriteLock(self) nogil:
        PyUniverseSpec_StateLock(self.universe_spec, -2)
        self.lock_state = 0

#
# Base class for trajectory generators that call the C-level
# energy evaluators. It implements a mechanism that makes such
# generators thread-safe.
#
cdef class EnergyBasedTrajectoryGenerator(TrajectoryGenerator):

    cdef evaluator_object
    cdef c_evaluator_object
    cdef PyFFEvaluatorObject *evaluator

    def __init__(self, universe, options, name):
        TrajectoryGenerator.__init__(self, universe, options, name)
        self.evaluator_object = None
        self.c_evaluator_object = None
        self.evaluator = NULL
        
    cdef void initializeTrajectoryActions(self) except *:
        TrajectoryGenerator.initializeTrajectoryActions(self)
        # Construct a C evaluator object for the force field, using
        # the specified number of threads or the default value
        nt = self.getOption('threads')
        self.evaluator_object = self.universe.energyEvaluator(threads=nt)
        self.c_evaluator_object = self.evaluator_object.CEvaluator()
        self.evaluator = <PyFFEvaluatorObject*>self.c_evaluator_object        

    cdef void finalizeTrajectoryActions(self, int last_step,
                                        int error=False) except *:
        TrajectoryGenerator.finalizeTrajectoryActions(self, last_step, error)
        self.evaluator = NULL
        self.c_evaluator_object = None
        self.evaluator_object = None

    cdef void calculateEnergies(self, N.ndarray[double, ndim=2] conf_array,
                                energy_data *energy, int small_change=0) \
                            except *:
        if self.lock_state == -1:
            PyUniverseSpec_StateLock(self.universe_spec, -2)
        if self.lock_state != 1:
            PyUniverseSpec_StateLock(self.universe_spec, 1)
            
        self.evaluator.eval_func(self.evaluator, energy, <PyArrayObject *>conf_array,
                                 small_change)
        if self.lock_state != 1:
            PyUniverseSpec_StateLock(self.universe_spec, 2)
        if self.lock_state == -1:
            PyUniverseSpec_StateLock(self.universe_spec, -1)

#
# No U-Turn Sampler "integrator"
#
cdef class NUTSIntegrator(EnergyBasedTrajectoryGenerator):

  """
  NUTS molecular dynamics integrator
  """
  
  cdef float RT
  cdef N.ndarray x, v, g, m
  cdef energy_data energy

  def __init__(NUTSIntegrator self, universe, **options):
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
      EnergyBasedTrajectoryGenerator.__init__(
          self, universe, options, "NUTS integrator")
      # Supported features: none for the moment, to keep it simple
      self.features = []

  default_options = {'first_step': 0, 'steps': 100, 'delta_t': 1.*Units.fs,
                     'background': False, 'threads': None,
                     'actions': []}

#  available_data = ['configuration', 'velocities', 'gradients',
#                    'energy', 'time']
#
#  restart_data = ['configuration', 'velocities', 'energy']

  # Cython compiler directives set for efficiency:
  # - No bound checks on index operations
  # - No support for negative indices
  # - Division uses C semantics
  @cython.boundscheck(False)
  @cython.wraparound(False)
  @cython.cdivision(True)
  cdef start(NUTSIntegrator self):

    cdef float delta_t, ke
    cdef int natoms, max_steps, elapsed_steps
    cdef Py_ssize_t i_atm, i_dim

    cdef float joint, logu, e_m
    cdef N.ndarray[double, ndim=2] xminus, xplus, vminus, vplus, x_m
    cdef int j, n, steps_m
    
    # For dual averaging
    cdef float alpha, delta, gamma, kappa, mu, delta_t_bar, Hbar, eta
    cdef int nalpha, t0, m

    # Initialize the velocity
    self.universe.initializeVelocitiesToTemperature(self.getOption('T'))

    # Gather state variables and parameters
    configuration = self.universe.configuration()
    velocities = self.universe.velocities()
    gradients = ParticleVector(self.universe)
    masses = self.universe.masses()
    delta_t = self.getOption('delta_t')
    max_steps = self.getOption('steps')
    natoms = self.universe.numberOfAtoms()
    
    self.RT = R*self.getOption('T')

    if 'normalize' in self.call_options.keys():
      normalize = self.getOption('normalize')
    else:
      normalize = False

    if 'adapt' in self.call_options.keys():
      adapt = self.getOption('adapt')
    else:
      adapt = False
      
    if 'delta' in self.call_options.keys():
      delta = self.getOption('delta')
    else:
      delta = 0.6

    # For efficiency, the Cython code works at the array
    # level rather than at the ParticleProperty level.
    self.x = configuration.array
    self.v = velocities.array
    self.g = gradients.array
    self.m = N.repeat(N.expand_dims(masses.array,1),3,axis=1)

    # Weight matrix for velocity assignment
    sigma_MB = N.sqrt((self.getOption('T')*Units.k_B)/self.m)

    # Ask for energy gradients to be calculated and stored in
    # the array g. Force constants are not requested.
    self.energy.gradients = <void *>self.g
    self.energy.gradient_fn = NULL
    self.energy.force_constants = NULL
    self.energy.fc_fn = NULL
    
    self.initializeTrajectoryActions()
    
    # Parameters for the dual averaging algorithm.
    if adapt:
      gamma = 0.05
      t0 = 10
      kappa = 0.75
      mu = N.log(10*delta_t)
      # Initialize dual averaging algorithm.
      delta_t_bar = 1.
      Hbar = 0.
    else:
      delta_t_bar = delta_t
      
    self.calculateEnergies(self.x, &self.energy, 0)

    # Initialize the NUTS algorithm
    xs = [N.copy(self.x)]
    energies = [self.energy.energy]
    
    # Initialize the next sample.
    # If all else fails, the next sample is the previous sample.
    x_m = N.copy(self.x)
    e_m = self.energy.energy

    # Main integration loop
    elapsed_steps = 0
    m = 1
    while elapsed_steps < max_steps:
      # Resample velocities
      self.v = N.multiply(sigma_MB,N.random.randn(natoms,3))
      ke = 0.5*N.sum(N.multiply(self.m,N.multiply(self.v,self.v)))

      # Joint log-probabiity of positions and velocities
      joint = -(e_m + ke)/self.RT

      # Resample u ~ uniform([0, N.exp(joint)]).
      # Equivalent to (log(u) - joint) ~ exponential(1).
      logu = joint - N.random.exponential(1)
      
      # Initialize tree.
      xminus = N.copy(x_m)
      xplus = N.copy(x_m)
      vminus = N.copy(self.v)
      vplus = N.copy(self.v)

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
        if N.random.rand()<0.5:
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
        if (sprime and (N.random.rand() < float(nprime)/n)):
          x_m = xprime
          e_m = eprime
        # Update number of valid points we've seen.
        n += nprime
        # Increment depth
        j += 1
        # Decide if it's time to stop
        s = sprime and \
            self.stop_criterion(xminus, xplus, vminus, vplus) and \
            ((elapsed_steps + steps_m*2) < max_steps)

      # Keep track of acceptance statistics
      eta = 1./float(m+t0)
      Hbar = (1-eta)*Hbar + eta*(delta-alpha/nalpha)

      # Adapt the time step
      if adapt:
        delta_t = N.exp(mu - N.sqrt(m)/gamma*Hbar)
        eta = m**-kappa
        delta_t_bar = N.exp((1-eta)*N.log(delta_t_bar) + eta*N.log(delta_t))

      if normalize:
        self.universe.setConfiguration(Configuration(self.universe, x_m))
        self.universe.normalizePosition()
        x_m = self.universe.configuration().array

      xs.append(N.copy(x_m))
      energies.append(e_m)
      m += 1
      elapsed_steps += steps_m

    self.universe.setConfiguration(Configuration(self.universe, x_m))
    return (xs, energies, Hbar, delta_t_bar)

  # The main recursion
  def build_tree(NUTSIntegrator self, float logu, int j, float delta_t, int steps, float joint_o):
    cdef float ke, joint, eprime, eprime2
    cdef N.ndarray[double, ndim=2] xminus, xplus, vminus, vplus, xprime, xprime2, NA0, NA1
    
    cdef float alphaprime, alphaprime2
    cdef int nalphaprime, nalphaprime2

    if (j==0):
      # Base case: Take a single leapfrog step
      # First half-step
      self.v += -0.5*delta_t*N.divide(self.g,self.m)
      self.x += delta_t*self.v
      # Mid-step energy calculation
      self.foldCoordinatesIntoBox()
      self.calculateEnergies(self.x, &self.energy, 1)
      # Second half-step
      self.v += -0.5*delta_t*N.divide(self.g,self.m)
      ke = 0.5*N.sum(N.multiply(self.m,N.multiply(self.v,self.v)))
      steps += 1

      eprime = self.energy.energy
      joint = -(eprime+ke)/self.RT
      # Is the new point in the slice?
      nprime = (logu < joint)
      # Is the simulation wildly inaccurate
      sprime = N.abs(joint_o - joint) < 1000
      # Compute the acceptance probability
      alphaprime = min(1, N.exp(joint - joint_o))
      # Set the return values---minus=plus for all things here, since the
      # "tree" is of depth 0.

#      if N.isnan(joint):
#        print '*** JOINT PROBABILITY IS NAN***'

      return (N.copy(self.x), N.copy(self.v),
              N.copy(self.x), N.copy(self.v), 
              N.copy(self.x),
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
           (N.random.rand() < float(nprime2) / (nprime + nprime2)):
            xprime = xprime2
            eprime = eprime2
        # Update the number of valid points.
        nprime = nprime + nprime2
        # Update the stopping criterion.
        sprime = sprime and sprime2 and \
          self.stop_criterion(xminus, xplus, vminus, vplus)
        # Update the acceptance probability statistics
        alphaprime = alphaprime + alphaprime2
        nalphaprime = nalphaprime + nalphaprime2
      return (xminus, vminus, xplus, vplus, xprime,
              eprime, nprime, sprime, steps,
              alphaprime, nalphaprime)

  def stop_criterion(NUTSIntegrator self, xminus, xplus, vminus, vplus):
    cdef N.ndarray[double] thetavec
    thetavec = N.ravel(xplus-xminus)
    return (N.dot(thetavec,N.ravel(vminus))>0) and \
           (N.dot(thetavec,N.ravel(vplus))>0)

