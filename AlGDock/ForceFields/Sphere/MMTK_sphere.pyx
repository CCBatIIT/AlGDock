#
# Get all the required declarations
#
include "MMTK/python.pxi"
include "MMTK/numeric.pxi"
include "MMTK/core.pxi"
include "MMTK/universe.pxi"
include 'MMTK/forcefield.pxi'

cdef extern from "math.h":
    double sqrt(double)
    
import numpy as N
cimport numpy as N

#
# The force field term implementation.
# The rules:
#
# - The class must inherit from EnergyTerm.
#
# - EnergyTerm.__init__() must be called with the arguments
#   shown here. The third argument is the name of the EnergyTerm
#   object, the fourth a tuple of the names of all the terms it
#   implements (one object can implement several terms).
#   The assignment to self.eval_func is essential, without it
#   any energy evaluation will crash.
#
# - The function "evaluate" must have exactly the parameter
#   list given in this example.
#
cdef class SphereTerm(EnergyTerm):

    # Parameters
    cdef vector3 center
    cdef double max_R
    cdef char* name
    # Useful variables
    cdef double k
    cdef N.ndarray masses
    cdef double totalMass, max_R2
    cdef int natoms

    def __init__(self, universe,
                 N.ndarray[N.float64_t] center, max_R, name):
                 
        EnergyTerm.__init__(self, universe,
                            name, (name,))

        self.eval_func = <void *>SphereTerm.evaluate

        # Parameters
        cdef N.ndarray[N.float64_t] center_array
        center_array = center
        for i in range(3):
          self.center[i] = (<double *>center_array.data)[i]
        self.max_R = max_R
        self.name = name
        # Useful variables
        self.k = 10000 # kJ/mol nm**2
        self.masses = universe.masses().array
        self.totalMass = sum(self.masses)
        self.natoms = len(self.masses)
        self.max_R2 = self.max_R*self.max_R

    # The function evaluate is called for every single energy
    # evaluation and should therefore be optimized for speed.
    # Its first argument is the global energy evaluator object,
    # which is needed only for parallelized energy terms.
    # The second argument is a C structure that contains all the
    # input data, in particular the particle configuration.
    # The third argument is a C structure that contains the
    # energy term fields and gradient arrays for storing the results.
    # For details, see MMTK_forcefield.pxi.
    cdef void evaluate(self, PyFFEvaluatorObject *eval,
                       energy_spec *input, energy_data *energy):
        # Input
        cdef vector3 *coordinates
        cdef N.float_t *masses
        # Output
        cdef vector3 *gradients
        # Processing
        cdef vector3 com, p # Center of mass & position in sphere
        cdef double kw, kw0, kw1, kw2, overMax, r2
        cdef int atom_index, i, j # Loop variables
        
        coordinates = <vector3 *>input.coordinates.data
        masses = <N.float_t *>self.masses.data

        # Initialize variables
        energy.energy_terms[self.index] = 0
        if energy.gradients != NULL:
          gradients = <vector3 *>(<PyArrayObject *> energy.gradients).data

        # First get the center of mass
        com[0] = com[1] = com[2] = 0.
        for atom_index in range(self.natoms):
          for i in range(3):
            com[i] += masses[atom_index]*coordinates[atom_index][i]
        for i in range(3):
          com[i] /= self.totalMass

        # Position in the sphere
        for i in range(3):
          p[i] = com[i] - self.center[i]

        # Energy and gradients orthogonal to the principal axis
        r2 = p[0]*p[0] + p[1]*p[1] + p[2]*p[2]
        if (r2>self.max_R2):
          r = sqrt(r2)
          overMax = (r-self.max_R)
          energy.energy_terms[self.index] += self.k*overMax**2/2
          if energy.gradients != NULL:
            kw = self.k*overMax/r/self.totalMass
            kw0 = kw*p[0]
            kw1 = kw*p[1]
            kw2 = kw*p[2]
            for atom_index in range(self.natoms):
              gradients[atom_index][0] += kw0*masses[atom_index]
              gradients[atom_index][1] += kw1*masses[atom_index]
              gradients[atom_index][2] += kw2*masses[atom_index]
