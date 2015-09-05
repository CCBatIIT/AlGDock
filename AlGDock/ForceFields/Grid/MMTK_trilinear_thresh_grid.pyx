# Cython force field implementation for trilinear grid

#
# Get all the required declarations
#
include "MMTK/python.pxi"
include "MMTK/numeric.pxi"
include "MMTK/core.pxi"
include "MMTK/universe.pxi"
include 'MMTK/forcefield.pxi'

import gzip

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
cdef class TrilinearThreshGridTerm(EnergyTerm):
    cdef char* grid_name
    cdef np.ndarray scaling_factor, vals, counts, spacing, hCorner
    cdef int npts, nyz, natoms
    cdef float_t energy_thresh, strength, k

    # The __init__ method remembers parameters and loads the potential
    # file. Note that EnergyTerm.__init__ takes care of storing the
    # name and the universe object.
    def __init__(self, universe, spacing, counts, vals, strength,
                 scaling_factor, grid_name, energy_thresh):

        EnergyTerm.__init__(self, universe,
                            grid_name, (grid_name,))
        self.eval_func = <void *>TrilinearThreshGridTerm.evaluate

        self.strength = strength
        self.scaling_factor = np.array(scaling_factor, dtype=float)
        self.natoms = len(self.scaling_factor)
        self.grid_name = grid_name
        self.energy_thresh = energy_thresh
        
        self.spacing = spacing
        self.counts = counts
        self.vals = vals
        self.nyz = self.counts[1]*self.counts[2]
        self.hCorner = np.array((self.spacing[0]*(self.counts[0]-1),
                        self.spacing[1]*(self.counts[1]-1),
                        self.spacing[2]*(self.counts[2]-1)), dtype=float)
        # To keep atoms within the grid
        self.k = 10000. # kJ/mol nm**2

    # This method is called for every single energy evaluation, so make
    # it as efficient as possible. The parameters do_gradients and
    # do_force_constants are flags that indicate if gradients and/or
    # force constants are requested.
    cdef void evaluate(self, PyFFEvaluatorObject *eval,
                       energy_spec *input, energy_data *energy):

        # Input
        cdef vector3 *coordinates
        cdef float_t *scaling_factor
        cdef float_t *vals
        cdef int_t *counts
        cdef float_t *spacing
        cdef float_t *hCorner
        # Output
        cdef float_t gridEnergy
        cdef vector3 *gradients
        # Processing
        cdef int i, ix, iy, iz, atom_index
        cdef float_t vmmm, vmmp, vmpm, vmpp, vpmm, vpmp, vppm, vppp
        cdef float_t vmm, vmp, vpm, vpp, vm, vp
        cdef float_t fx, fy, fz, ax, ay, az
        cdef float_t dvdx, dvdy, dvdz
        cdef float_t Eo, denergy_thresh, coshEo, sinhEo

        gridEnergy = 0
        coordinates = <vector3 *>input.coordinates.data

        # Pointers to numpy arrays for faster indexing
        scaling_factor = <float_t *>self.scaling_factor.data
        vals = <float_t *>self.vals.data
        counts = <int_t *>self.counts.data
        spacing = <float_t *>self.spacing.data
        hCorner = <float_t *>self.hCorner.data

        # Initialize variables
        if energy.gradients != NULL:
          gradients = <vector3 *>(<PyArrayObject *> energy.gradients).data
      
        indicies = [ai for ai in range(self.natoms) if self.scaling_factor[ai]!=0]
        for atom_index in indicies:
          # Check to make sure coordinate is in grid
          if (coordinates[atom_index][0]>0 and 
              coordinates[atom_index][1]>0 and 
              coordinates[atom_index][2]>0 and
              coordinates[atom_index][0]<hCorner[0] and
              coordinates[atom_index][1]<hCorner[1] and
              coordinates[atom_index][2]<hCorner[2]):

            # Index within the grid
            ix = int(coordinates[atom_index][0]/spacing[0])
            iy = int(coordinates[atom_index][1]/spacing[1])
            iz = int(coordinates[atom_index][2]/spacing[2])
            
            i = ix*self.nyz + iy*counts[2] + iz

            # Corners of the box surrounding the point
            vmmm = vals[i]
            vmmp = vals[i+1]
            vmpm = vals[i+counts[2]]
            vmpp = vals[i+counts[2]+1]

            vpmm = vals[i+self.nyz]
            vpmp = vals[i+self.nyz+1]
            vppm = vals[i+self.nyz+counts[2]]
            vppp = vals[i+self.nyz+counts[2]+1]
            
            # Fraction within the box
            fx = (coordinates[atom_index][0] - (ix*spacing[0]))/spacing[0]
            fy = (coordinates[atom_index][1] - (iy*spacing[1]))/spacing[1]
            fz = (coordinates[atom_index][2] - (iz*spacing[2]))/spacing[2]
            
            # Fraction ahead
            ax = 1 - fx
            ay = 1 - fy
            az = 1 - fz
      
            # Trilinear interpolation for energy
            vmm = az*vmmm + fz*vmmp
            vmp = az*vmpm + fz*vmpp
            vpm = az*vpmm + fz*vpmp
            vpp = az*vppm + fz*vppp
            
            vm = ay*vmm + fy*vmp
            vp = ay*vpm + fy*vpp
            
            Eo = (ax*vm + fx*vp)

            # This is slow
            gridEnergy += scaling_factor[atom_index]*self.energy_thresh*np.tanh(Eo/self.energy_thresh)

            # These overflow
            # en2x = np.exp(-2.*Eo)
            # gridEnergy += scaling_factor[atom_index]*self.energy_thresh*(1.-en2x)/(1.+en2x)
            
            # sinhEo = np.sinh(Eo)
            # coshEo = np.cosh(Eo)
            # gridEnergy += scaling_factor[atom_index]*self.energy_thresh*sinhEo/coshEo
          
            if energy.gradients != NULL:
              # x coordinate
              dvdx = -vm + vp
              # y coordinate
              dvdy = (-vmm + vmp)*ax + (-vpm + vpp)*fx
              # z coordinate
              dvdz = ((-vmmm + vmmp)*ay + (-vmpm + vmpp)*fy)*ax + ((-vpmm + vpmp)*ay + (-vppm + vppp)*fy)*fx
              
              # This is slow
              denergy_thresh = (1./np.cosh(Eo/self.energy_thresh))
              denergy_thresh = 1./denergy_thresh

              # These overflow
              # denergy_thresh = 4.*en2x/(1.+en2x)
              # denergy_thresh = 1./coshEo/coshEo
            
              gradients[atom_index][0] += self.strength*scaling_factor[atom_index]*denergy_thresh*dvdx/spacing[0]
              gradients[atom_index][1] += self.strength*scaling_factor[atom_index]*denergy_thresh*dvdy/spacing[1]
              gradients[atom_index][2] += self.strength*scaling_factor[atom_index]*denergy_thresh*dvdz/spacing[2]
          else:
            for i in range(3):
              if (coordinates[atom_index][i]<0):
                gridEnergy += self.k*coordinates[atom_index][i]**2/2.
                if energy.gradients != NULL:
                  gradients[atom_index][i] += self.k*coordinates[atom_index][i]
              elif (coordinates[atom_index][i]>hCorner[i]):
                gridEnergy += self.k*(coordinates[atom_index][i]-hCorner[i])**2/2.
                if energy.gradients != NULL:
                  gradients[atom_index][i] += self.k*(coordinates[atom_index][i]-hCorner[i])

        energy.energy_terms[self.index] = gridEnergy*self.strength
