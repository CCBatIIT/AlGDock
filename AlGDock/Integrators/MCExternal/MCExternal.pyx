# This module implements Markov chain Monte Carlo moves
# for translation and rotation
#
# Written by David Minh
#

import numpy as N
cimport numpy as N
import cython

cimport MMTK_trajectory_generator
from MMTK import Units, ParticleProperties
import MMTK_trajectory
import MMTK_forcefield

include "MMTK/python.pxi"
include "MMTK/numeric.pxi"
include "MMTK/core.pxi"
include "MMTK/universe.pxi"
include "MMTK/trajectory.pxi"
include "MMTK/forcefield.pxi"


#
# MCExternal integrator
#
cdef class MCExternalIntegrator(MMTK_trajectory_generator.EnergyBasedTrajectoryGenerator):

  """
  MCExternal integrator

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
        self, universe, options, "MCExternal integrator")
    # Supported features: none for the moment, to keep it simple
    self.features = []

  default_options = {'first_step': 0, 'steps': 100, 'delta_t': 1.*Units.fs,
                     'background': False, 'threads': None,
                     'actions': []}

  available_data = ['configuration', 'velocities', 'gradients',
                    'energy', 'time']

  restart_data = ['configuration', 'velocities', 'energy']

  # Cython compiler directives set for efficiency:
  # - No bound checks on index operations
  # - No support for negative indices
  # - Division uses C semantics
  @cython.boundscheck(False)
  @cython.wraparound(False)
  @cython.cdivision(True)
  cdef start(self):
    cdef N.ndarray[double, ndim=2] x, v, g, dv
    cdef N.ndarray[double, ndim=1] m
    cdef energy_data energy
    cdef double time, delta_t, ke
    cdef int natoms, nsteps, step
    cdef Py_ssize_t i, j, k

    # Check if velocities have been initialized
    if self.universe.velocities() is None:
        raise ValueError("no velocities")

    # Gather state variables and parameters
    configuration = self.universe.configuration()
    velocities = self.universe.velocities()
    gradients = ParticleProperties.ParticleVector(self.universe)
    masses = self.universe.masses()
    delta_t = self.getOption('delta_t')
    nsteps = self.getOption('steps')
    natoms = self.universe.numberOfAtoms()

    self.RT = 8.3144621*Units.J/Units.mol/Units.K*self.getOption('T')

    # Seed the random number generator
    if 'trials' in self.call_options.keys():
      trials = self.getOption('trials')
    else:
      trials = 20

    if 'step_size' in self.call_options.keys():
      step_size = self.getOption('step_size')
    else:
      step_size = 0.1

    # Seed the random number generator
    if 'seed' in self.call_options.keys():
      N.random.seed(self.getOption('seed'))
    else:
      N.random.seed()

    # For efficiency, the Cython code works at the array
    # level rather than at the ParticleProperty level.
    x = configuration.array
    v = velocities.array
    g = gradients.array
    m = masses.array
    dv = N.zeros((natoms, 3), N.float)

    # Ask for energy gradients to be calculated and stored in
    # the array g. Force constants are not requested.
    energy.gradients = <void *>g
    energy.gradient_fn = NULL
    energy.force_constants = NULL
    energy.fc_fn = NULL

    # Declare the variables accessible to trajectory actions.
    self.declareTrajectoryVariable_double(
        &time, "time", "Time: %lf\n", time_unit_name, PyTrajectory_Time)
    self.declareTrajectoryVariable_array(
        v, "velocities", "Velocities:\n", velocity_unit_name,
        PyTrajectory_Velocities)
    self.declareTrajectoryVariable_array(
        g, "gradients", "Energy gradients:\n", energy_gradient_unit_name,
        PyTrajectory_Gradients)
    self.declareTrajectoryVariable_double(
        &energy.energy,"potential_energy", "Potential energy: %lf\n",
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

    # Preparation: copy initial state energy and configuration
    acc = 0.
    self.calculateEnergies(x, &energy, 0)
    eo = 1.*energy.energy
    com = self.universe.centerOfMass().array

    for c in range(trials):
      step = N.random.randn(size=3)*step_size
      xn = N.dot((x - com), self.random_rotate()) + com + step
      self.calculateEnergies(xn, &energy, 0)
      if ((energy.energy<eo) or \
          (N.random.random()<N.exp(-(energy.energy-eo)/self.RT))):
        acc += 1
        x = N.copy(xn)
        eo = 1.*energy.energy
        com += step

    # Release the write lock.
    self.releaseWriteLock()

    # Finalize all trajectory actions (close files etc.)
    self.finalizeTrajectoryActions(nsteps)

    # TO DO: Work on returns
#    return (xs, energies)

  def random_rotate(self):
    """
    Return a random rotation matrix
    """
    u = N.random.uniform(size=3)

    # Random quaternion
    q = N.array([N.sqrt(1-u[0])*N.sin(2*N.pi*u[1]),
                 N.sqrt(1-u[0])*N.cos(2*N.pi*u[1]),
                 N.sqrt(u[0])*N.sin(2*N.pi*u[2]),
                 N.sqrt(u[0])*N.cos(2*N.pi*u[2])])
    
    # Convert the quaternion into a rotation matrix 
    rotMat = N.array([[q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3],
                       2*q[1]*q[2] - 2*q[0]*q[3],
                       2*q[1]*q[3] + 2*q[0]*q[2]],
                      [2*q[1]*q[2] + 2*q[0]*q[3],
                       q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3],
                       2*q[2]*q[3] - 2*q[0]*q[1]],
                      [2*q[1]*q[3] - 2*q[0]*q[2],
                       2*q[2]*q[3] + 2*q[0]*q[1],
                       q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3]]])
    return rotMat
