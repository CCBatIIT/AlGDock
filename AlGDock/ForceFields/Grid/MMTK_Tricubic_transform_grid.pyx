# Implementation of the Tricubic Interpolation Algorith for the AGD Docking
# TODO: Does not seem to compile properly

# Author: Luis Antonio Leite Francisco da Costa

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
cdef class TricubicTransformGridTerm(EnergyTerm):
    cdef char* grid_name
    cdef np.ndarray scaling_factor, vals, counts, spacing, hCorner
    cdef int npts, nyz, natoms
    cdef float_t strength, inv_power, inv_power_m1, k
    # The __init__ method remembers parameters and loads the potential
    # file. Note that EnergyTerm.__init__ takes care of storing the
    # name and the universe object.

    cdef float_t tricubicInterpolate(self, list data, list nkpoints, float_t x):
        
    cdef float_t dfdx[8]
    cdef float_t dfdy[8]
    cdef float_t dfdz[8]    
    cdef float_t d2fdxdy[8]
    cdef float_t d2fdxdz[8]
    cdef float_t d2fdydz[8]
    cdef float_t d3fdxdydz[8]

    cdef int n1 = extract<int>(nkpoints[0])
    cdef int n2 = extract<int>(nkpoints[1])
    cdef int n3 = extract<int>(nkpoints[2])

    cdef fptype x = extract<fptype>(data[0])
    cdef fptype y = extract<fptype>(data[1])
    cdef fptype z = extract<fptype>(data[2])

    # Determine the relative position in the box enclosed by nearest data points
    cdef fptype dx = fmod(x/spacing, n1)
    cdef fptype dy = fmod(y/spacing, n2)
    cdef fptype dz = fmod(z/spacing, n3)        

     # Periodicity
    if(dx < 0): dx += n1
    if(dy < 0): dy += n2
    if(dz < 0): dz += n3
  
    # Calculate lower-bound grid indices
    cdef int xi = (int)floor(dx); 
    cdef int yi = (int)floor(dy);
    cdef int zi = (int)floor(dz);

    # Calculate the derivatives estimatives
    
    # Values of df/dx in each corner in a array
    cdef float_t dfdx(self, float_t p[4][4][4]):
        cdef float_t arr[8]
        arr[0] = p[index(xi,yi,zi)]
        arr[1] = p[index(xi+1,yi,zi)]
        arr[2] = p[index(xi,yi+1,zi)]
        arr[3] = p[index(xi+1,yi+1,zi)]
        arr[4] = p[index(xi,yi,zi+1)]
        arr[5] = p[index(xi+1,yi,zi+1)]
        arr[6] = p[index(xi,yi+1,zi+1)]
        arr[7] = p[index(xi+1,yi+1,zi+1)]
        return self.dfdx(arr)

    # Values of df/dy in each corner
    cdef float_t dfdy(self, float_t p[4][4][4]):
        cdef float_t arr[8]
        arr[0] = 0.5*(p[index(xi+1,yi,zi)]-p[index(xi-1,yi,zi)])
        arr[1] = 0.5*(p[index(xi+2,yi,zi)]-p[index(xi,yi,zi)])
        arr[2] = 0.5*(p[index(xi+1,yi+1,zi)]-p[index(xi-1,yi+1,zi)])
        arr[3] = 0.5*(p[index(xi+2,yi+1,zi)]-p[index((xi,yi+1,zi)])
        arr[4] = 0.5*(p[index(xi+1,yi,zi+1)]-p[index(xi-1,yi,zi+1)])
        arr[5] = 0.5*(p[index(xi+2,yi,zi+1)]-p[index(xi,yi,zi+1)])
        arr[6] = 0.5*(p[index(xi+1,yi+1,zi+1)]-p[index(xi-1,yi+1,zi+1)])
        arr[7] = 0.5*(p[index(xi+2,yi+1,zi+1)]-p[index(xi,yi+1,zi+1)])
        return self.dfdy(arr)

    # Values of df/dz in each corner
    cdef float_t dfdy(self, float_t p[4][4][4]):
        cdef float_t arr[8]
        arr[0] = 0.5*(p[index(xi,yi,zi+1)]-p[index(xi,yi,zi-1)])
        arr[1] = 0.5*(p[index(xi+1,yi,zi+1)]-p[index(xi+1,yi,zi-1)])
        arr[2] = 0.5*(p[index(xi,yi+1,zi+1)]-p[index(xi,yi+1,zi-1)])
        arr[3] = 0.5*(p[index(xi+1,yi+1,zi+1)]-p[index(xi+1,yi+1,zi-1)])
        arr[4] = 0.5*(p[index(xi,yi,zi+2)]-p[index(xi,yi,zi)])
        arr[5] = 0.5*(p[index(xi+1,yi,zi+2)]-p[index(xi+1,yi,zi)])
        arr[6] = 0.5*(p[index(xi,yi+1,zi+2)]-p[index(xi,yi+1,zi)])
        arr[7] = 0.5*(p[index(xi+1,yi+1,zi+2)]-p[index(xi+1,yi+1,zi)])
        return self.dfdz(arr)

    # Values of d2f/dxdy in each corner
    cdef float_t d2fdxdy(self, float_t p[4][4][4]):
        cdef float_t arr[8]
        arr[0] = 0.25*(p[index(xi+1,yi+1,zi)]-p[index(xi-1,yi+1,zi)]-p[index(xi+1,yi-1,zi)]+p[index(xi-1,yi-1,zi)])
        arr[1] = 0.25*(p[index(xi+2,yi+1,zi)]-p[index(xi,yi+1,zi)]-p[index(xi+2,yi-1,zi)]+p[index(xi,yi-1,zi)])
        arr[2] =  0.25*(p[index(xi+1,yi+2,zi)]-p[index(xi-1,yi+2,zi)]-p[index(xi+1,yi,zi)]+p[index(xi-1,yi,zi)])
        arr[3] = 0.25*(p[index(xi+2,yi+2,zi)]-p[index(xi,yi+2,zi)]-p[index(xi+2,yi,zi)]+p[index(xi,yi,zi)])
        arr[4] = 0.25*(p[index(xi+1,yi+1,zi+1)]-p[index(xi-1,yi+1,zi+1)]-p[index(xi+1,yi-1,zi+1)]+p[index(xi-1,yi-1,zi+1)])
        arr[5] = 0.25*(p[index(xi+2,yi+1,zi+1)]-p[index(xi,yi+1,zi+1)]-p[index(xi+2,yi-1,zi+1)]+p[index(xi,yi-1,zi+1)])
        arr[6] = 0.25*(p[index(xi+1,yi+2,zi+1)]-p[index(xi-1,yi+2,zi+1)]-p[index(xi+1,yi,zi+1)]+p[index(xi-1,yi,zi+1)])
        arr[7] = 0.25*(p[index(xi+2,yi+2,zi+1)]-p[index(xi,yi+2,zi+1)]-p[index(xi+2,yi,zi+1)]+p[index(xi,yi,zi+1)])
        return self.d2fdxdy(arr)

    # Values of d2f/dxdz in each corner
    cdef float_t d2fdxdz(self, float_t p[4][4][4]):
        cdef float_t arr[8]
        arr[0] = 0.25*(p[index(xi+1,yi,zi+1)]-p[index(xi-1,yi,zi+1)]-p[index(xi+1,yi,zi-1)]+p[index(xi-1,yi,zi-1)])
        arr[1] = 0.25*(p[index(xi+2,yi,zi+1)]-p[index(xi,yi,zi+1)]-p[index(xi+2,yi,zi-1)]+p[index(xi,yi,zi-1)])
        arr[2] = 0.25*(p[index(xi+1,yi+1,zi+1)]-p[index(xi-1,yi+1,zi+1)]-p[index(xi+1,yi+1,zi-1)]+p[index(xi-1,yi+1,zi-1)])
        arr[3] = 0.25*(p[index(xi+2,yi+1,zi+1)]-p[index(xi,yi+1,zi+1)]-p[index(xi+2,yi+1,zi-1)]+p[index(xi,yi+1,zi-1)])
        arr[4] = 0.25*(p[index(xi+1,yi,zi+2)]-p[index(xi-1,yi,zi+2)]-p[index(xi+1,yi,zi)]+p[index(xi-1,yi,zi)])
        arr[5] = 0.25*(p[index(xi+2,yi,zi+2)]-p[index(xi,yi,zi+2)]-p[index(xi+2,yi,zi)]+p[index(xi,yi,zi)])
        arr[6] = 0.25*(p[index(xi+1,yi+1,zi+2)]-p[index(xi-1,yi+1,zi+2)]-p[index(xi+1,yi+1,zi)]+p[index(xi-1,yi+1,zi)])
        arr[7] = 0.25*(p[index(xi+2,yi+1,zi+2)]-p[index(xi,yi+1,zi+2)]-p[index(xi+2,yi+1,zi)]+p[index(xi,yi+1,zi)])
        return self.d2fdxdz(arr)

    # Values of d2f/dydz in each corner
    cdef float_t d2fdydz(self, float_t p[4][4][4]):
        cdef float_t arr[8]
        arr[0] = 0.25*(p[index(xi,yi+1,zi+1)]-p[index(xi,yi-1,zi+1)]-p[index(xi,yi+1,zi-1)]+p[index(xi,yi-1,zi-1)])
        arr[1] = 0.25*(p[index(xi+1,yi+1,zi+1)]-p[index(xi+1,yi-1,zi+1)]-p[index(xi+1,yi+1,zi-1)]+p[index(xi+1,yi-1,zi-1)])
        arr[2] = 0.25*(p[index(xi,yi+2,zi+1)]-p[index(xi,yi,zi+1)]-p[index(xi,yi+2,zi-1)]+p[index(xi,yi,zi-1)])
        arr[3] = 0.25*(p[index(xi+1,yi+2,zi+1)]-d[index(xi+1,yi,zi+1)]-p[index(xi+1,yi+2,zi-1)]+p[index(xi+1,yi,zi-1)])
        arr[4] = 0.25*(p[index(xi,yi+1,zi+2)]-p[index(xi,yi-1,zi+2)]-p[index(xi,yi+1,zi)]+p[index(xi,yi-1,zi)])
        arr[5] = 0.25*(p[index(xi+1,yi+1,zi+2)]-p[index(xi+1,yi-1,zi+2)]-p[index(xi+1,yi+1,zi)]+p[index(xi+1,yi-1,zi)])
        arr[6] = 0.25*(p[index(xi,yi+2,zi+2)]-p[index(xi,yi,zi+2)]-p[index(xi,yi+2,zi)]+p[index(xi,yi,zi)])
        arr[7] = 0.25*(p[index(xi+1,yi+2,zi+2)]-p[index(xi+1,yi,zi+2)]-p[index(xi+1,yi+2,zi)]+p[index(xi+1,yi,zi)])
        return self.d2fdydz(arr)

    # Values of d3f/dxdydz in each corner
    cdef float_t d3fdxdydz(self, float_t p[4][4][4]):
        cdef float_t arr[8]
        arr[0] = 0.125*(p[index(xi+1,yi+1,zi+1)]-p[index(xi-1,yi+1,zi+1)]-p[index(xi+1,yi-1,zi+1)]+p[index(xi-1,yi-1,zi+1)]-p[index(xi+1,yi+1,zi-1)]+p[index(xi-1,yi+1,zi-1)]+p[index(xi+1,yi-1,zi-1)]-p[index(xi-1,yi-1,zi-1)])
        arr[1] = 0.125*(p[index(xi+2,yi+1,zi+1)]-p[index(xi,yi+1,zi+1)]-p[index(xi+2,yi-1,zi+1)]+p[index(xi,yi-1,zi+1)]-p[index(xi+2,yi+1,zi-1)]+p[index(xi,yi+1,zi-1)]+p[index(xi+2,yi-1,zi-1)]-p[index(xi,yi-1,zi-1)])
        arr[2] = 0.125*(p[index(xi+1,yi+2,zi+1)]-p[index(xi-1,yi+2,zi+1)]-p[index(xi+1,yi,zi+1)]+p[index(xi-1,yi,zi+1)]-p[index(xi+1,yi+2,zi-1)]+p[index(xi-1,yi+2,zi-1)]+p[index(xi+1,yi,zi-1)]-p[index(xi-1,yi,zi-1)])
        arr[3] = 0.125*(p[index(xi+2,yi+2,zi+1)]-p[index(xi,yi+2,zi+1)]-p[index(xi+2,yi,zi+1)]+p[index(xi,yi,zi+1)]-p[index(xi+2,yi+2,zi-1)]+p[index(xi,yi+2,zi-1)]+p[index(xi+2,yi,zi-1)]-p[index(xi,yi,zi-1)])
        arr[4] = 0.125*(p[index(xi+1,yi+1,zi+2)]-p[index(xi-1,yi+1,zi+2)]-p[index(xi+1,yi-1,zi+2)]+p[index(xi-1,yi-1,zi+2)]-p[index(xi+1,yi+1,zi)]+p[index(xi-1,yi+1,zi)]+p[index(xi+1,yi-1,zi)]-p[index(xi-1,yi-1,zi)])
        arr[5] = 0.125*(p[index(xi+2,yi+1,zi+2)]-p[index(xi,yi+1,zi+2)]-p[index(xi+2,yi-1,zi+2)]+p[index(xi,yi-1,zi+2)]-p[index(xi+2,yi+1,zi)]+p[index(xi,yi+1,zi)]+p[index(xi+2,yi-1,zi)]-p[index(xi,yi-1,zi)])
        arr[6] = 0.125*(p[index(xi+1,yi+2,zi+2)]-p[index(xi-1,yi+2,zi+2)]-p[index(xi+1,yi,zi+2)]+p[index(xi-1,yi,zi+2)]-p[index(xi+1,yi+2,zi)]+p[index(xi-1,yi+2,zi)]+p[index(xi+1,yi,zi)]-p[index(xi-1,yi,zi)])
        arr[7] = 0.125*(p[index(xi+2,yi+2,zi+2)]-p[index(xi,yi+2,zi+2)]-p[index(xi+2,yi,zi+2)]+p[index(xi,yi,zi+2)]-p[index(xi+2,yi+2,zi)]+p[index(xi,yi+2,zi)]+p[index(xi+2,yi,zi)]-p[index(xi,yi,zi)])
        return self.d3fddxdydz(arr)


    def __init__(self, universe, spacing, counts, vals, strength,
                 scaling_factor, grid_name, max_val):
        print "------------test start---------------"
        EnergyTerm.__init__(self, universe,
                            grid_name, (grid_name,))
        self.eval_func = <void *>BSplineGridTerm.evaluate

        self.strength = strength
        self.scaling_factor = np.array(scaling_factor, dtype=float)
        self.natoms = len(self.scaling_factor)
        self.grid_name = grid_name
        self.inv_power = float(inv_power)
        self.inv_power_m1 = inv_power - 1.

        self.spacing = spacing
        self.counts = counts
        self.vals = vals
        self.nyz = self.counts[1]*self.counts[2]
        self.hCorner = np.array((self.spacing[0]*(self.counts[0]-1),
                        self.spacing[1]*(self.counts[1]-1),
                        self.spacing[2]*(self.counts[2]-1)), dtype=float)
        # To keep atoms within the grid
        self.k = 10000. # kJ/mol nm**2

        # "Cap" the grid values
        if max_val>0.0:
          self.vals = max_val*np.tanh(self.vals/max_val)
        else:
          self.vals = np.copy(vals)

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

            gridEnergy += scaling_factor[atom_index]*self.tricubicInterpolate(vertex,fx,fy,fz)
                # hessian funciton          ************************
            if energy.force_constants !=NULL:
               # x direction
              dvdxdx = self.dfdx(vertex)
              dvdxdy = self.df2dxdy(vertex)
              dvdxdz = self.df2dxdz(vertex)
              # y direction
              dvdydy = self.dfdy(vertex)
              dvdydz = self.df2dydz(vertex)
              # z direction
              dvdzdz = self.dfdz(vertex)
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
              dvdx = self.dfdx(vertex) + self.df2dxdy(vertex) + self.df2dxdz(vertex)
              # y coordinate
              dvdy = self.dfdy(vertex) + self.df2dxdy(vertex) + self.df2dydz(vertex)
              # z coordinate
              dvdz = self.dfdz(vertex) + self.df2dxdz(vertex) + self.df2dydz(vertex)
            
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

