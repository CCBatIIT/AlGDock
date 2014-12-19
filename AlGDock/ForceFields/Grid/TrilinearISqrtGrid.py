# Grid-based potential energy term.
# It can be added to any other force field.

from MMTK.ForceFields.ForceField import ForceField, EnergyTerm
from MMTK import ParticleScalar, ParticleVector, SymmetricPairTensor

try:
  from Scientific._vector import Vector
except:
  from Scientific.Geometry.VectorModule import Vector

from MMTK_trilinear_isqrt_grid import TrilinearISqrtGridTerm

class TrilinearISqrtGridForceField(ForceField):

    """
    Force field based on trilinear interpolation between points on a transformed 3D grid
    """

    def __init__(self, FN, strength, scaling_property, scaling_prefactor=None, grid_name='trilinear transformed grid', max_val=-1.0):
        """
        @param strength: the electric field vector
        @type strength: L{float}
        @scaling_property: the name of the atomic property in the database
                           that is used to retrieve the atomic scaling_factor.
                           The default is 'amber_charge', the charge property
                           for Amber94 and Amber99.
        @type scaling_property: C{str}
        @grid_name: a name for the grid
        @max_val: the maximum allowed value for a point on the grid.
                  A negative value means that there is no max.
        """
        # Store arguments that recreate the force field from a pickled
        # universe or from a trajectory.
        self.arguments = (FN, strength,
          scaling_property, scaling_prefactor, grid_name, max_val)
        # Initialize the ForceField class, giving a name to this one.
        ForceField.__init__(self, grid_name)
        # Store the parameters for later use.
        self.FN = FN
        self.strength = strength
        self.scaling_property = scaling_property
        self.grid_name = grid_name
        self.max_val = max_val
  
        import AlGDock.IO
        IO_Grid = AlGDock.IO.Grid()
        self.grid_data = IO_Grid.read(self.FN, multiplier=0.1)
        if not (self.grid_data['origin']==0.0).all():
          raise Exception('Trilinear grid origin in %s not at (0, 0, 0)!'%self.FN)

        # Make sure all grid values are positive
        if (self.grid_data['vals']>0).any():
          if (self.grid_data['vals']<0).any():
            raise Exception('All of the grid points do not have the same sign')
          else:
            neg_vals = False
        else:
          neg_vals = True
          self.grid_data['vals'] = -1*self.grid_data['vals']

        if scaling_prefactor is not None:
          self.scaling_prefactor = scaling_prefactor
        else:
          self.scaling_prefactor = -1. if neg_vals else 1.

    # The following method is called by the energy evaluation engine
    # to inquire if this force field term has all the parameters it
    # requires. This is necessary for interdependent force field
    # terms.
    def ready(self, global_data):
        return True
  
    # The following method is returns a dictionary of parameters for
    # the force field
    def evaluatorParameters(self, universe, subset1, subset2, global_data):
        return {self.grid_name+' file name': self.FN,
                self.grid_name+' strength': self.strength,
                self.grid_name+' scaling property': self.scaling_property,
                self.grid_name+' max val': self.max_val}

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
        # Collect the scaling_factor into an array
        scaling_factor = ParticleScalar(universe)
        for o in universe:
            for a in o.atomList():
                scaling_factor[a] = o.getAtomProperty(a, self.scaling_property)
        scaling_factor.scaleBy(self.scaling_prefactor)
        # Here we pass all the parameters to
        # the energy term code that handles energy calculations.
        return [TrilinearISqrtGridTerm(universe, \
          self.grid_data['spacing'], self.grid_data['counts'], \
          self.grid_data['vals'], \
          self.strength, scaling_factor, self.grid_name, self.max_val)]