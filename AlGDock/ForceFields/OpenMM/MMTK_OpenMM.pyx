# Cython force field implementation for OpenMM

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

import simtk.openmm
import simtk.openmm.app as OpenMM_app

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
cdef class OpenMMTerm(EnergyTerm):
    cdef char* implicitSolvent
    cdef np.ndarray prmtop_atom_order, inv_prmtop_atom_order
    cdef object prmtop, OMM_system, dummy_integrator, sim, context

    # The __init__ method remembers parameters and loads the potential
    # file. Note that EnergyTerm.__init__ takes care of storing the
    # name and the universe object.
    def __init__(self, universe, prmtopFN, \
        prmtop_atom_order, inv_prmtop_atom_order, implicitSolvent):

      EnergyTerm.__init__(self, universe, 'OpenMM', ('OpenMM',))
      self.eval_func = <void *>OpenMMTerm.evaluate

      self.prmtop_atom_order = np.array(prmtop_atom_order, dtype=float)
      self.inv_prmtop_atom_order = np.array(inv_prmtop_atom_order, dtype=float)

      self.prmtop = OpenMM_app.AmberPrmtopFile(prmtopFN)
      self.OMM_system = self.prmtop.createSystem(\
        nonbondedMethod=OpenMM_app.NoCutoff, \
        constraints=None, \
        implicitSolvent={
          'OpenMM_Gas':None,
          'OpenMM_GBn':OpenMM_app.GBn,
          'OpenMM_GBn2':OpenMM_app.GBn2,
          'OpenMM_HCT':OpenMM_app.HCT,
          'OpenMM_OBC1':OpenMM_app.OBC1,
          'OpenMM_OBC2':OpenMM_app.OBC2}[implicitSolvent])
      self.dummy_integrator = simtk.openmm.LangevinIntegrator(300*simtk.unit.kelvin, \
        1/simtk.unit.picosecond, 0.002*simtk.unit.picoseconds)
      # platform = simtk.openmm.Platform.getPlatformByName('CPU')
      self.sim = OpenMM_app.Simulation(self.prmtop.topology, self.OMM_system, self.dummy_integrator)
      self.context = self.sim.context

    # This method is called for every single energy evaluation, so make
    # it as efficient as possible. The parameters do_gradients and
    # do_force_constants are flags that indicate if gradients and/or
    # force constants are requested.
    cdef void evaluate(self, PyFFEvaluatorObject *eval,
                       energy_spec *input, energy_data *energy):

      # Input
      cdef vector3 *coordinates, *gradients
      cdef int natoms, ai, ain
      cdef np.ndarray[double, ndim=2] openmm_coordinates

      coordinates = <vector3 *>input.coordinates.data
      natoms = input.coordinates.dimensions[0]

      openmm_coordinates = np.zeros((natoms,3))
      for ai in range(natoms):
        ain = self.prmtop_atom_order[ai]
        openmm_coordinates[ai,0] = coordinates[ain][0]
        openmm_coordinates[ai,1] = coordinates[ain][1]
        openmm_coordinates[ai,2] = coordinates[ain][2]
      self.context.setPositions(openmm_coordinates)

      state = self.context.getState(getEnergy=True, getForces=(energy.gradients!=NULL))
      energy.energy_terms[self.index] = state.getPotentialEnergy()/simtk.unit.kilojoule*simtk.unit.mole

      if energy.gradients != NULL:
        openmm_gradients = state.getForces()/simtk.unit.kilojoule*simtk.unit.mole*simtk.unit.nanometer
        gradients = <vector3 *>(<PyArrayObject *> energy.gradients).data
        for ai in range(natoms):
          ain = self.inv_prmtop_atom_order[ai]
          gradients[ai][0] -= openmm_gradients[ain][0]
          gradients[ai][1] -= openmm_gradients[ain][1]
          gradients[ai][2] -= openmm_gradients[ain][2]
