cimport numpy as np
from cpython cimport bool

from libc.math cimport sin
from libc.math cimport cos
from libc.math cimport sqrt
from libc.math cimport acos
from libc.math cimport atan2

cdef inline dotp(double[:] v1, double[:] v2)
cdef inline norm2(double[:] v)
cdef inline normalize(np.ndarray[np.double_t] v1)
cdef inline cross(double[:] v1, double[:] v2)
cdef inline distance(double[:] p1, double[:] p2)
cdef inline angle(np.ndarray[np.double_t] p1, np.ndarray[np.double_t] p2, \
    np.ndarray[np.double_t] p3)
cdef dihedral(np.ndarray[np.double_t] p1, np.ndarray[np.double_t] p2, \
    np.ndarray[np.double_t] p3, np.ndarray[np.double_t] p4)
cdef extended_coordinates(np.ndarray[np.double_t] p1, \
    np.ndarray[np.double_t] p2, np.ndarray[np.double_t] p3)
cdef rotate123(np.ndarray[np.double_t] p1, \
    np.ndarray[np.double_t] p2, np.ndarray[np.double_t] p3, \
    np.ndarray[np.double_t] origin, \
    double phi, double theta, double omega)

cdef class converter:

  cdef object molecule
  cdef readonly int natoms, ntorsions
  cdef np.ndarray rootInd, _torsionIndL, _firstTorsionTInd

  cpdef BAT(self, double[:,:] XYZ, bool extended)
  cpdef Cartesian(self, double[:] BAT)

