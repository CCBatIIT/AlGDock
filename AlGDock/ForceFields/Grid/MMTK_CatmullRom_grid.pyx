# Cython force field implementation for Catmull-Rom grid

# TODO: The energy as a function of position is not smooth; it has a bug.

#
# Get all the required declarations
#
include "MMTK/python.pxi"
include "MMTK/numeric.pxi"
include "MMTK/core.pxi"
include "MMTK/universe.pxi"
include 'MMTK/forcefield.pxi'
from Scientific.Geometry import Vector, ex, ey, ez

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
cdef class CatmullRomGridTerm(EnergyTerm):
    cdef char* grid_name
    cdef np.ndarray scaling_factor, vals, counts, spacing, hCorner
    cdef int npts, nyz, natoms
    cdef float_t strength, k
    # The __init__ method remembers parameters and loads the potential
    # file. Note that EnergyTerm.__init__ takes care of storing the
    # name and the universe object.

    
    cdef float_t splineInterpolate(self,float_t p[4],float_t x):
        return p[0]+.5*x*(-p[0]+p[2]+x*(-4.*p[0]+7.*p[1]-2.*p[2]-p[3]+x*(3.*p[0]-5.*p[1]+p[2]+p[3])))
    cdef float_t bisplineInterpolate(self,float_t p[4][4],float_t x,float_t y):
        cdef float_t arr[4]
        arr[0] = self.splineInterpolate(p[0], y)
        arr[1] = self.splineInterpolate(p[1], y)
        arr[2] = self.splineInterpolate(p[2], y)
        arr[3] = self.splineInterpolate(p[3], y)
        return self.splineInterpolate(arr, x)
    cdef float_t trisplineInterpolate(self,float_t p[4][4][4],float_t x,float_t y,float_t z):
        cdef float_t arr[4]
        arr[0] = self.bisplineInterpolate(p[0], y, z)
        arr[1] = self.bisplineInterpolate(p[1], y, z)
        arr[2] = self.bisplineInterpolate(p[2], y, z)
        arr[3] = self.bisplineInterpolate(p[3], y, z)
        return self.splineInterpolate(arr, x)


    cdef float_t derivateOfIntp(self,float_t p[4],float_t x):
        return -.5*p[0]+.5*p[2]+x*(-4.*p[0]+7.*p[1]-2.*p[2]-p[3]+1.5*x*(3.*p[0]-5.*p[1]+p[2]+p[3]))

# the following functions are used to realize the gradients(first dirivative)
    cdef float_t derivateOfIntp_X(self,float_t p[4][4][4],float_t x,float_t y,float_t z):
        cdef float_t arr[4]
        arr[0] = self.bisplineInterpolate(p[0], y, z)
        arr[1] = self.bisplineInterpolate(p[1], y, z)
        arr[2] = self.bisplineInterpolate(p[2], y, z)
        arr[3] = self.bisplineInterpolate(p[3], y, z)
        return self.derivateOfIntp(arr, x)

    cdef float_t derivateOfIntp_Y_2(self,float_t p[4][4],float_t x,float_t y):
        cdef float_t arr[4]
        arr[0] = self.splineInterpolate(p[0], y)
        arr[1] = self.splineInterpolate(p[1], y)
        arr[2] = self.splineInterpolate(p[2], y)
        arr[3] = self.splineInterpolate(p[3], y)
        return self.derivateOfIntp(arr, x)

    cdef float_t derivateOfIntp_Y(self,float_t p[4][4][4],float_t x,float_t y,float_t z):
        cdef float_t arr[4]
        arr[0] = self.derivateOfIntp_Y_2(p[0], y, z)
        arr[1] = self.derivateOfIntp_Y_2(p[1], y, z)
        arr[2] = self.derivateOfIntp_Y_2(p[2], y, z)
        arr[3] = self.derivateOfIntp_Y_2(p[3], y, z)
        return self.splineInterpolate(arr, x)
    cdef float_t derivateOfIntp_Z_2(self,float_t p[4][4],float_t x,float_t y):
        cdef float_t arr[4]
        arr[0] = self.derivateOfIntp(p[0], y)
        arr[1] = self.derivateOfIntp(p[1], y)
        arr[2] = self.derivateOfIntp(p[2], y)
        arr[3] = self.derivateOfIntp(p[3], y)
        return self.splineInterpolate(arr, x)
    cdef float_t derivateOfIntp_Z(self,float_t p[4][4][4],float_t x,float_t y,float_t z):
        cdef float_t arr[4]
        arr[0] = self.derivateOfIntp_Z_2(p[0], y, z)
        arr[1] = self.derivateOfIntp_Z_2(p[1], y, z)
        arr[2] = self.derivateOfIntp_Z_2(p[2], y, z)
        arr[3] = self.derivateOfIntp_Z_2(p[3], y, z)
        return self.splineInterpolate(arr, x)

# the following functions are used to realize the hessian functions(second derivative)

    cdef float_t derivateOfIntp_mm(self,float_t p[4],float_t x):
        return -4.*p[0]+7.*p[1]-2.*p[2]-p[3]+3*x*(3.*p[0]-5.*p[1]+p[2]+p[3])

# calculate the dvdxdx
    cdef float_t derivateOfIntp_XX(self,float_t p[4][4][4],float_t x,float_t y,float_t z):
        cdef float_t arr[4]
        arr[0] = self.bisplineInterpolate(p[0], y, z)
        arr[1] = self.bisplineInterpolate(p[1], y, z)
        arr[2] = self.bisplineInterpolate(p[2], y, z)
        arr[3] = self.bisplineInterpolate(p[3], y, z)
        return self.derivateOfIntp_mm(arr, x)

# calculate the dvdxdy
    cdef float_t derivateOfIntp_XY_2(self,float_t p[4][4],float_t x,float_t y):
        cdef float_t arr[4]
        arr[0] = self.splineInterpolate(p[0], y)
        arr[1] = self.splineInterpolate(p[1], y)
        arr[2] = self.splineInterpolate(p[2], y)
        arr[3] = self.splineInterpolate(p[3], y)
        return self.derivateOfIntp(arr, x)

    cdef float_t derivateOfIntp_XY(self,float_t p[4][4][4],float_t x,float_t y,float_t z):
        cdef float_t arr[4]
        arr[0] = self.derivateOfIntp_XY_2(p[0], y, z)
        arr[1] = self.derivateOfIntp_XY_2(p[1], y, z)
        arr[2] = self.derivateOfIntp_XY_2(p[2], y, z)
        arr[3] = self.derivateOfIntp_XY_2(p[3], y, z)
        return self.derivateOfIntp(arr, x)

# calculate the dvdxdz
    cdef float_t derivateOfIntp_XZ_2(self,float_t p[4][4],float_t x,float_t y):
        cdef float_t arr[4]
        arr[0] = self.derivateOfIntp(p[0], y)
        arr[1] = self.derivateOfIntp(p[1], y)
        arr[2] = self.derivateOfIntp(p[2], y)
        arr[3] = self.derivateOfIntp(p[3], y)
        return self.splineInterpolate(arr, x)

    cdef float_t derivateOfIntp_XZ(self,float_t p[4][4][4],float_t x,float_t y,float_t z):
        cdef float_t arr[4]
        arr[0] = self.derivateOfIntp_XZ_2(p[0], y, z)
        arr[1] = self.derivateOfIntp_XZ_2(p[1], y, z)
        arr[2] = self.derivateOfIntp_XZ_2(p[2], y, z)
        arr[3] = self.derivateOfIntp_XZ_2(p[3], y, z)
        return self.derivateOfIntp(arr, x)

# calculate the dvdydy
    cdef float_t derivateOfIntp_YY_2(self,float_t p[4][4],float_t x,float_t y):
        cdef float_t arr[4]
        arr[0] = self.splineInterpolate(p[0], y)
        arr[1] = self.splineInterpolate(p[1], y)
        arr[2] = self.splineInterpolate(p[2], y)
        arr[3] = self.splineInterpolate(p[3], y)
        return self.derivateOfIntp_mm(arr, x)

    cdef float_t derivateOfIntp_YY(self,float_t p[4][4][4],float_t x,float_t y,float_t z):
        cdef float_t arr[4]
        arr[0] = self.derivateOfIntp_YY_2(p[0], y, z)
        arr[1] = self.derivateOfIntp_YY_2(p[1], y, z)
        arr[2] = self.derivateOfIntp_YY_2(p[2], y, z)
        arr[3] = self.derivateOfIntp_YY_2(p[3], y, z)
        return self.splineInterpolate(arr, x)

# calculate the dvdydz
    cdef float_t derivateOfIntp_YZ_2(self,float_t p[4][4],float_t x,float_t y):
        cdef float_t arr[4]
        arr[0] = self.derivateOfIntp(p[0], y)
        arr[1] = self.derivateOfIntp(p[1], y)
        arr[2] = self.derivateOfIntp(p[2], y)
        arr[3] = self.derivateOfIntp(p[3], y)
        return self.derivateOfIntp(arr, x)

    cdef float_t derivateOfIntp_YZ(self,float_t p[4][4][4],float_t x,float_t y,float_t z):
        cdef float_t arr[4]
        arr[0] = self.derivateOfIntp_YZ_2(p[0], y, z)
        arr[1] = self.derivateOfIntp_YZ_2(p[1], y, z)
        arr[2] = self.derivateOfIntp_YZ_2(p[2], y, z)
        arr[3] = self.derivateOfIntp_YZ_2(p[3], y, z)
        return self.splineInterpolate(arr, x)

# calculate the dvdzdz
    cdef float_t derivateOfIntp_ZZ_2(self,float_t p[4][4],float_t x,float_t y):
        cdef float_t arr[4]
        arr[0] = self.derivateOfIntp_mm(p[0], y)
        arr[1] = self.derivateOfIntp_mm(p[1], y)
        arr[2] = self.derivateOfIntp_mm(p[2], y)
        arr[3] = self.derivateOfIntp_mm(p[3], y)
        return self.splineInterpolate(arr, x)

    cdef float_t derivateOfIntp_ZZ(self,float_t p[4][4][4],float_t x,float_t y,float_t z):
        cdef float_t arr[4]
        arr[0] = self.derivateOfIntp_ZZ_2(p[0], y, z)
        arr[1] = self.derivateOfIntp_ZZ_2(p[1], y, z)
        arr[2] = self.derivateOfIntp_ZZ_2(p[2], y, z)
        arr[3] = self.derivateOfIntp_ZZ_2(p[3], y, z)
        return self.splineInterpolate(arr, x)



    def __init__(self, universe, spacing, counts, vals, strength,
                 scaling_factor, grid_name):
        print "------------test start---------------"
        EnergyTerm.__init__(self, universe,
                            grid_name, (grid_name,))
        self.eval_func = <void *>CatmullRomGridTerm.evaluate

        self.strength = strength
        self.scaling_factor = np.array(scaling_factor, dtype=float)
        self.natoms = len(self.scaling_factor)
        self.grid_name = grid_name

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
        cdef np.ndarray[float_t, ndim=4] force_constants
        # Processing

        cdef int i, ix, iy, iz, atom_index

        #Initialize the 4*4*4 matrix to store the energy value
        cdef float_t vertex[4][4][4]

        #Initialize a hessian matrix
        #cdef float_t hessian[3][3]

        cdef float_t dvdx, dvdy, dvdz
        cdef float_t dvdxdx, dvdxdy, dvdxdz,dvdydy,dvdydz,dvdzdz

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

        # Initialize variables
        if energy.force_constants != NULL:
            force_constants = <np.ndarray[float_t, ndim=4] >(<PyArrayObject *> energy.force_constants)

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
            ix = int(coordinates[atom_index][0]/spacing[0]-1)
            iy = int(coordinates[atom_index][1]/spacing[1]-1)
            iz = int(coordinates[atom_index][2]/spacing[2]-1)
            
            i = ix*self.nyz + iy*counts[2] + iz


            for x in range(0,4):
                for y in range(0,4):
                    for z in range(0,4):
                        vertex[x][y][z]=vals[i+x*self.nyz+y*counts[2]+z]

           # Fraction within the box
            fx = (coordinates[atom_index][0] - (ix*spacing[0]))/spacing[0]
            fy = (coordinates[atom_index][1] - (iy*spacing[1]))/spacing[1]
            fz = (coordinates[atom_index][2] - (iz*spacing[2]))/spacing[2]

            gridEnergy += scaling_factor[atom_index]*self.trisplineInterpolate(vertex,fx,fy,fz)
                # hessian funciton          ************************
            if energy.force_constants !=NULL:
               # x direction
              dvdxdx = self.derivateOfIntp_XX(vertex,fx,fy,fz)
              dvdxdy = self.derivateOfIntp_XY(vertex,fx,fy,fz)
              dvdxdz = self.derivateOfIntp_XZ(vertex,fx,fy,fz)
              # y direction
              dvdydy = self.derivateOfIntp_YY(vertex,fx,fy,fz)
              dvdydz = self.derivateOfIntp_YZ(vertex,fx,fy,fz)
              # z direction
              dvdzdz = self.derivateOfIntp_ZZ(vertex,fx,fy,fz)
              force_constants[atom_index][0][atom_index][0] +=  self.strength*scaling_factor[atom_index]*dvdxdx/spacing[0]/spacing[0]
              force_constants[atom_index][1][atom_index][1] +=  self.strength*scaling_factor[atom_index]*dvdydy/spacing[1]/spacing[1]
              force_constants[atom_index][2][atom_index][2] +=  self.strength*scaling_factor[atom_index]*dvdzdz/spacing[2]/spacing[2]
              force_constants[atom_index][1][atom_index][0] +=  self.strength*scaling_factor[atom_index]*dvdxdy/spacing[0]/spacing[1]
              force_constants[atom_index][2][atom_index][0] +=  self.strength*scaling_factor[atom_index]*dvdxdz/spacing[0]/spacing[2]
              force_constants[atom_index][2][atom_index][1] +=  self.strength*scaling_factor[atom_index]*dvdydz/spacing[1]/spacing[2]
              force_constants[atom_index][0][atom_index][1] = force_constants[atom_index][1][atom_index][0]
              force_constants[atom_index][0][atom_index][2] = force_constants[atom_index][2][atom_index][0]
              force_constants[atom_index][1][atom_index][2] = force_constants[atom_index][2][atom_index][1]

              #*****************************************
            if energy.gradients != NULL:
              # x coordinate
              dvdx = self.derivateOfIntp_X(vertex,fx,fy,fz)
              # y self.coordinate
              dvdy = self.derivateOfIntp_Y(vertex,fx,fy,fz)
              # z coordinate
              dvdz = self.derivateOfIntp_Z(vertex,fx,fy,fz)
            
              gradients[atom_index][0] += self.strength*scaling_factor[atom_index]*dvdx/spacing[0]
              gradients[atom_index][1] += self.strength*scaling_factor[atom_index]*dvdy/spacing[1]
              gradients[atom_index][2] += self.strength*scaling_factor[atom_index]*dvdz/spacing[2]
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
