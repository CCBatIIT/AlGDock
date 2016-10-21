# Cython force field implementation for OBC

### distutils: language = c++
### distutils: sources = ReferenceObc.cpp, ObcParameters.cpp, ReferenceForce.cpp

#
# Get all the required declarations
#
include "MMTK/python.pxi"
include "MMTK/numeric.pxi"
include "MMTK/core.pxi"
include "MMTK/universe.pxi"
include 'MMTK/forcefield.pxi'

import numpy as np
cimport numpy as np

ctypedef np.float_t float_t
ctypedef np.int_t int_t

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
cdef class OBCTerm(EnergyTerm):
    cdef np.ndarray charges, atomicRadii, scaleFactors
    cdef int natoms

    # The __init__ method remembers parameters and loads the potential
    # file. Note that EnergyTerm.__init__ takes care of storing the
    # name and the universe object.
    def __init__(self, universe, charges, atomicRadii, scaleFactors, name):

        EnergyTerm.__init__(self, universe,
                            name, (name,))
        self.eval_func = <void *>OBCTerm.evaluate

        self.charges = np.array(charges, dtype=float)
        self.atomicRadii = np.array(atomicRadii, dtype=float)
        self.scaleFactors = np.array(scaleFactors, dtype=float)
        self.natoms = len(self.scaling_factor)

    # This method is called for every single energy evaluation, so make
    # it as efficient as possible. The parameters do_gradients and
    # do_force_constants are flags that indicate if gradients and/or
    # force constants are requested.
    cdef void evaluate(self, PyFFEvaluatorObject *eval,
                       energy_spec *input, energy_data *energy):

        # Input
        cdef vector3 *coordinates
        cdef float_t *charges
        cdef float_t *atomicRadii
        cdef float_t *scaleFactors
        # Output
        cdef float_t gridEnergy
        cdef vector3 *gradients

        gridEnergy = 0
        coordinates = <vector3 *>input.coordinates.data

        # Pointers to numpy arrays for faster indexing
        charges = <float_t *>self.charges.data
        atomicRadii = <float_t *>self.atomicRadii.data
        scaleFactors = <float_t *>self.scaleFactors.data

        # Initialize variables
        if energy.gradients != NULL:
          gradients = <vector3 *>(<PyArrayObject *> energy.gradients).data
      
