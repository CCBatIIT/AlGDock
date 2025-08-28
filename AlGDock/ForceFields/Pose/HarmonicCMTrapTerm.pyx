# Center-of-mass restraints
#
# Written by Konrad Hinsen

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
# Trap potential (harmonic restraint to a fixed point in space)
#
cdef class HarmonicCMTrapTerm(EnergyTerm):

    cdef N.ndarray atom_indices, masses
    cdef vector3 ref
    cdef double k, r0

    def __init__(self, universe,
                 N.ndarray[N.int_t] atom_indices,
                 N.ndarray[N.float64_t] masses,
                 reference, force_constant, r0=0.4):
        cdef N.ndarray[N.float64_t] ref_array
        EnergyTerm.__init__(self, universe,
                            "harmonic_cm_trap", ("harmonic_cm_trap",))
        self.eval_func = <void *>HarmonicCMTrapTerm.evaluate
        self.atom_indices = atom_indices
        self.masses = masses
        ref_array = reference.array
        for i in range(3):
            self.ref[i] = (<double *>ref_array.data)[i]
        self.k = force_constant
        self.r0 = r0

    cdef void evaluate(self, PyFFEvaluatorObject *eval,
                       energy_spec *input, energy_data *energy):
        cdef vector3 *coordinates
        cdef vector3 *gradients
        cdef N.ndarray[N.int_t] atom_indices
        cdef N.ndarray[N.float_t] masses
        cdef N.ndarray[N.float_t, ndim=4] fc
        cdef vector3 cm, d
        cdef double m, lsq, kw
        cdef int i, j, n, offset

        coordinates = <vector3 *>input.coordinates.data
        atom_indices = self.atom_indices
        masses = self.masses

        # This code will never be run for a periodic universe, so
        # it's safe to ignore periodic boundary conditions.
        m = 0.
        cm[0] = cm[1] = cm[2] = 0.
        for i in atom_indices:
            m += masses[i]
            for j in range(3):
                cm[j] += masses[i]*coordinates[i][j]
        for j in range(3):
            cm[j] /= m

        for j in range(3):
            d[j] = cm[j]-self.ref[j]
        lsq = d[0]*d[0]+d[1]*d[1]+d[2]*d[2]
        dist = sqrt(lsq)

        # Apply flat bottom
        if dist > self.r0:

            energy.energy_terms[self.index] = 0.5*self.k*lsq
            energy.energy_terms[self.virial_index] -= 2.*self.k*lsq
    
            if energy.gradients != NULL:
                gradients = <vector3 *>(<PyArrayObject *> energy.gradients).data
                kw = 2.*self.k/m
                for i in atom_indices:
                    for j in range(3):
                        gradients[i][j] += kw*d[j]*masses[i]
    
            if energy.force_constants != NULL:
                fc = <object>energy.force_constants
                kw = 2.*self.k/(m*m)
                for i in atom_indices:
                    for j in atom_indices:
                        for n in range(3):
                            fc [i, n, j, n] += kw*masses[i]*masses[j]
        else:
            energy.energy_terms[self.index] = 0                   