# OpenMM-based potential energy term.
# It can be added to any other force field.

from MMTK.ForceFields.ForceField import ForceField, EnergyTerm
from MMTK import ParticleScalar, ParticleVector, SymmetricPairTensor

import simtk.openmm
import simtk.openmm.app as OpenMM_app
import numpy as np

try:
  from Scientific._vector import Vector
except:
  from Scientific.Geometry.VectorModule import Vector

class OpenMMForceField(ForceField):
  """
  Force field based on OpenMM molecular mechanics
  """

  def __init__(self, prmtopFN, \
      prmtop_atom_order, inv_prmtop_atom_order, implicitSolvent='OpenMM_OBC2'):
    """
    @name: a name for the grid
    @implicitSolvent: the type of implicit solvent, which can be ['OpenMM_Gas','OpenMM_GBn', 'OpenMM_GBn2', 'OpenMM_HCT', 'OpenMM_OBC1', 'OpenMM_OBC2'].
    """
    if not implicitSolvent in \
        ['OpenMM_Gas','OpenMM_GBn', 'OpenMM_GBn2', \
         'OpenMM_HCT', 'OpenMM_OBC1', 'OpenMM_OBC2']:
      raise Exception('Implicit solvent not recognized')

    ForceField.__init__(self, implicitSolvent) # Initialize the ForceField class

    # Store arguments that recreate the force field from a pickled
    # universe or from a trajectory.
    self.arguments = (prmtopFN, prmtop_atom_order, inv_prmtop_atom_order, implicitSolvent)
    
  # The following method is called by the energy evaluation engine
  # to inquire if this force field term has all the parameters it
  # requires. This is necessary for interdependent force field
  # terms.
  def ready(self, global_data):
    return True

  # The following method is returns a dictionary of parameters for
  # the force field
  def evaluatorParameters(self, universe, subset1, subset2, global_data):
    return {'prmtopFN':self.arguments[0], \
            'prmtop_atom_order':self.arguments[1], \
            'inv_prmtop_atom_order':self.arguments[2], \
            'implicitSolvent':self.arguments[3]}

  # The following method is called by the energy evaluation engine
  # to obtain a list of the evaluator objects
  # that handle the calculations.
  def evaluatorTerms(self, universe, subset1, subset2, global_data):
    # The energy for subsets is defined as consisting only
    # of interactions within that subset, so the contribution
    # of an external field is zero. Therefore we just return
    # an empty list of energy terms.
    if subset1 is not None or subset2 is not None:
      return []
    from MMTK_OpenMM import OpenMMTerm
    return [OpenMMTerm(universe, self.arguments[0], \
      self.arguments[1], self.arguments[2], \
      self.arguments[3])]

#  class OpenMMTerm(EnergyTerm):
#      # The __init__ method only remembers parameters. Note that
#      # EnergyTerm.__init__ takes care of storing the name and the
#      # universe object.
#      def __init__(self, universe, prmtopFN, \
#          prmtop_atom_order, inv_prmtop_atom_order, implicitSolvent):
#        
#        EnergyTerm.__init__(self, 'OpenMM', universe)
#        self.prmtop_atom_order = prmtop_atom_order
#        self.inv_prmtop_atom_order = inv_prmtop_atom_order
#
#        self.prmtop = OpenMM_app.AmberPrmtopFile(prmtopFN)
#        self.OMM_system = self.prmtop.createSystem(\
#          nonbondedMethod=OpenMM_app.NoCutoff, \
#          constraints=None, \
#          implicitSolvent={
#            'OpenMM_Gas':None,
#            'OpenMM_GBn':OpenMM_app.GBn,
#            'OpenMM_GBn2':OpenMM_app.GBn2,
#            'OpenMM_HCT':OpenMM_app.HCT,
#            'OpenMM_OBC1':OpenMM_app.OBC1,
#            'OpenMM_OBC2':OpenMM_app.OBC2}[implicitSolvent])
#        self.dummy_integrator = simtk.openmm.LangevinIntegrator(300*simtk.unit.kelvin, \
#          1/simtk.unit.picosecond, 0.002*simtk.unit.picoseconds)
#        # platform = simtk.openmm.Platform.getPlatformByName('CPU')
#        self.sim = OpenMM_app.Simulation(self.prmtop.topology, self.OMM_system, self.dummy_integrator)
#        self.context = self.sim.context
#
#      # This method is called for every single energy evaluation, so make
#      # it as efficient as possible. The parameters do_gradients and
#      # do_force_constants are flags that indicate if gradients and/or
#      # force constants are requested.
#      def evaluate(self, configuration, do_gradients, do_force_constants):
#          results = {}
#
#          self.context.setPositions(configuration.array[self.prmtop_atom_order,:])
#          state = self.context.getState(getEnergy=True, getForces=do_gradients)
#
#          results['energy'] = state.getPotentialEnergy()/simtk.unit.kilojoule*simtk.unit.mole
#          if do_gradients:
#              gradients = ParticleVector(self.universe)
#              gradients.array = -np.array(state.getForces()/simtk.unit.kilojoule*simtk.unit.mole*simtk.unit.nanometer)[self.inv_prmtop_atom_order,:]
#              results['gradients'] = gradients
#          if do_force_constants:
#              pass
#          return results