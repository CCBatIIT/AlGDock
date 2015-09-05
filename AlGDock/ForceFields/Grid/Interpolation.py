# Grid-based potential energy term.
# It can be added to any other force field.

from MMTK.ForceFields.ForceField import ForceField, EnergyTerm
from MMTK import ParticleScalar, ParticleVector, SymmetricPairTensor
from collections import OrderedDict

try:
  from Scientific._vector import Vector
except:
  from Scientific.Geometry.VectorModule import Vector

class InterpolationForceField(ForceField):
  """
  Force fields that interpolate between points on the 3D grid
  """

  def __init__(self, FN,
    name='Interpolation',
    interpolation_type='Trilinear',
    strength=1.0,
    scaling_property='amber_charge',
    scaling_prefactor=None,
    inv_power=None,
    grid_thresh=-1.0,
    energy_thresh=-1.0):
    """
    @FN: the file name.
    @name: a name for the grid
    @interpolation_type: the type of interpolation, which can be
      ['Trilinear','BSpline','CatmullRom', or 'Tricubic'].
    @param strength: scaling factor for all the energies and gradients.
    @type strength: L{float}
    @scaling_property: the name of the atomic property in the database.
      that is used to retrieve the atomic scaling_factor. 
      The default is 'amber_charge'.
    @scaling_prefactor: the atomic scaling_factor is scaled by this value.
    @inv_power: the inverse of the power by which grid points are transformed.
    @grid_thresh: the maximum allowed value for a point on the grid.
      A negative value means that there is no max.
    @energy_thresh: the maximum allowed value for the energy at any point.
      A negative value means that there is no max.
    @type scaling_property: C{str}
    """
    if not interpolation_type in \
        ['Trilinear','BSpline', 'CatmullRom', 'Tricubic']:
      raise Exception('Interpolation type not recognized')

    ForceField.__init__(self, name) # Initialize the ForceField class

    # Store arguments that recreate the force field from a pickled
    # universe or from a trajectory.
    self.arguments = (FN, name, interpolation_type, strength, \
      scaling_property, scaling_prefactor, \
      inv_power, grid_thresh, energy_thresh)
    
    self.params = OrderedDict()
    for key in ['FN','name','interpolation_type','strength','scaling_property',\
        'scaling_prefactor','inv_power','grid_thresh','energy_thresh']:
      self.params[key] = locals()[key]
    
    # Load the grid
    import AlGDock.IO
    IO_Grid = AlGDock.IO.Grid()
    self.grid_data = IO_Grid.read(FN, multiplier=0.1)
    if not (self.grid_data['origin']==0.0).all():
      raise Exception('Trilinear grid origin in %s not at (0, 0, 0)!'%FN)

    # Transform the grid
    neg_vals = False
    if inv_power is not None:
      # Make sure all grid values are positive
      if (self.grid_data['vals']>0).any():
        if (self.grid_data['vals']<0).any():
          raise Exception('All of the grid points do not have the same sign')
      else:
        neg_vals = True
        self.grid_data['vals'] = -1*self.grid_data['vals']

      # Transform all nonzero elements
      nonzero = self.grid_data['vals']!=0
      self.grid_data['vals'][nonzero] = self.grid_data['vals'][nonzero]**(1./inv_power)

    import numpy as np
    # "Cap" the grid values
    if grid_thresh>0.0:
      self.grid_data['vals'] = grid_thresh*np.tanh(self.grid_data['vals']/grid_thresh)

    if scaling_prefactor is None:
      self.params['scaling_prefactor'] = -1. if neg_vals else 1.

  # The following method is called by the energy evaluation engine
  # to inquire if this force field term has all the parameters it
  # requires. This is necessary for interdependent force field
  # terms.
  def ready(self, global_data):
    return True

  # The following method is returns a dictionary of parameters for
  # the force field
  def evaluatorParameters(self, universe, subset1, subset2, global_data):
    return self.params

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
        scaling_factor[a] = o.getAtomProperty(a, self.params['scaling_property'])
    scaling_factor.scaleBy(self.params['scaling_prefactor'])

    # Here we pass all the parameters to
    # the energy term code that handles energy calculations.
    if self.params['interpolation_type']=='Trilinear':
      if self.params['energy_thresh']>0:
        if self.params['inv_power'] is not None:
          raise NotImplementedError
        else:
          from MMTK_trilinear_thresh_grid import TrilinearThreshGridTerm
          return [TrilinearThreshGridTerm(universe, \
            self.grid_data['spacing'], self.grid_data['counts'], \
            self.grid_data['vals'], self.params['strength'], scaling_factor, \
            self.params['name'], self.params['energy_thresh'])]
      elif self.params['inv_power'] is not None:
        if self.params['inv_power']==-2:
          from MMTK_trilinear_isqrt_grid import TrilinearISqrtGridTerm
          return [TrilinearISqrtGridTerm(universe, \
            self.grid_data['spacing'], self.grid_data['counts'], \
            self.grid_data['vals'], self.params['strength'], scaling_factor, \
            self.params['name'])]
        else:
          from MMTK_trilinear_transform_grid import TrilinearTransformGridTerm
          return [TrilinearTransformGridTerm(universe, \
            self.grid_data['spacing'], self.grid_data['counts'], \
            self.grid_data['vals'], self.params['strength'], scaling_factor,
            self.params['name'], self.params['inv_power'])]
      else:
        from MMTK_trilinear_grid import TrilinearGridTerm
        return [TrilinearGridTerm(universe, \
          self.grid_data['spacing'], self.grid_data['counts'], \
          self.grid_data['vals'], self.params['strength'], scaling_factor, \
          self.params['name'])]
    elif self.params['interpolation_type']=='BSpline':
      if self.params['inv_power'] is not None:
        from MMTK_BSpline_transform_grid import BSplineTransformGridTerm
        return [BSplineTransformGridTerm(universe, \
          self.grid_data['spacing'], self.grid_data['counts'], \
          self.grid_data['vals'], self.params['strength'], scaling_factor, \
          self.params['name'], self.params['inv_power'])]
      else:
        from MMTK_BSpline_grid import BSplineGridTerm
        return [BSplineGridTerm(universe, \
          self.grid_data['spacing'], self.grid_data['counts'], \
          self.grid_data['vals'], self.params['strength'], scaling_factor, \
          self.params['name'])]
    elif self.params['interpolation_type']=='CatmullRom':
      if self.params['inv_power'] is not None:
        from MMTK_CatmullRom_transform_grid import CatmullRomTransformGridTerm
        return [CatmullRomTransformGridTerm(universe, \
          self.grid_data['spacing'], self.grid_data['counts'], \
          self.grid_data['vals'], self.params['strength'], scaling_factor, \
          self.params['name'], self.params['inv_power'])]
      else:
        from MMTK_CatmullRom_grid import CatmullRomGridTerm
        return [CatmullRomGridTerm(universe, \
          self.grid_data['spacing'], self.grid_data['counts'], \
          self.grid_data['vals'], self.params['strength'], scaling_factor, \
          self.params['name'])]
    elif self.params['interpolation_type']=='Tricubic':
      raise NotImplementedError
