#!/usr/bin/env python

import numpy as np
cimport numpy as np
import cython

ctypedef np.float_t float_t
ctypedef np.int_t int_t

@cython.boundscheck(False)
cpdef attempt_swaps(\
    state_inds, inv_state_inds, u_ij, \
    pairs_to_swap, int nattempts):
  cdef int attempt, t1, t2, a, b
  cdef float ddu
  for attempt in range(nattempts):
    for (t1,t2) in pairs_to_swap:
      a = inv_state_inds[t1]
      b = inv_state_inds[t2]
      ddu = -u_ij[a][b]-u_ij[b][a]+u_ij[a][a]+u_ij[b][b]
      if (ddu>0) or (np.random.uniform()<np.exp(ddu)):
        u_ij[a],u_ij[b] = u_ij[b],u_ij[a]
        state_inds[a],state_inds[b] = state_inds[b],state_inds[a]
        inv_state_inds[state_inds[a]],inv_state_inds[state_inds[b]] = \
          inv_state_inds[state_inds[b]],inv_state_inds[state_inds[a]]
  return state_inds, inv_state_inds