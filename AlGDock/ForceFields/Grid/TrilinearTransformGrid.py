# Grid-based potential energy term.
# It can be added to any other force field.

from MMTK.ForceFields.ForceField import ForceField, EnergyTerm
from MMTK import ParticleScalar, ParticleVector, SymmetricPairTensor

from Scientific._vector import Vector
from MMTK_trilinear_transform_grid import TrilinearTransformGridTerm

#  class TrilinearTransformGridTerm(EnergyTerm):
#      # The __init__ method remembers parameters and loads the potential
#      # file. Note that EnergyTerm.__init__ takes care of storing the
#      # name and the universe object.
#      def __init__(self, universe, spacing, counts, vals, strength,
#                   scaling_factor, inv_power, grid_name, max_val):
#        
#          EnergyTerm.__init__(self, grid_name, universe)
#          self.strength = strength
#          self.scaling_factor = scaling_factor
#          self.inv_power = float(inv_power)
#          self.inv_power_m1 = inv_power - 1.
#          self.grid_name = grid_name
#          self.max_val = max_val
#    
#          self.counts = counts
#          self.nyz = self.counts[1]*self.counts[2]
#          self.spacing = spacing
#          self.hCorner = [self.spacing[0]*(self.counts[0]-1),
#                          self.spacing[1]*(self.counts[1]-1),
#                          self.spacing[2]*(self.counts[2]-1)]
#
#          import numpy as np
#          # Transform the grid, except where the points are zero
#          self.grid_data['vals'][self.grid_data['vals']!=0] = self.grid_data['vals'][self.grid_data['vals']!=0]**(1./self.inv_power)
#        
#          # "Cap" the grid values
#          if max_val>0.0:
#            self.grid_data['vals'] = max_val*np.tanh(self.grid_data['vals']/max_val)
#    
#      # This method is called for every single energy evaluation, so make
#      # it as efficient as possible. The parameters do_gradients and
#      # do_force_constants are flags that indicate if gradients and/or
#      # force constants are requested.
#      def evaluate(self, configuration, do_gradients, do_force_constants):
#          energy = 0
#          gradients = ParticleVector(self.universe)
#        
#          for (atom_index,(x,y,z)) in zip(range(len(self.scaling_factor)),configuration.array):
#            if self.scaling_factor[atom_index]==0:
#              continue
#            
#            # Check to make sure coordinate is in grid
#            if (x>0 and y>0 and z>0 and
#                x<self.hCorner[0] and y<self.hCorner[1] and z<self.hCorner[2]):
#              # Index within the grid
#              ix = int(x/self.spacing[0])
#              iy = int(y/self.spacing[1])
#              iz = int(z/self.spacing[2])
#              
#              i = ix*self.nyz + iy*self.counts[2] + iz
#
#              # Corners of the box surrounding the point
#              vmmm = self.grid_data['vals'][i]
#              vmmp = self.grid_data['vals'][i+1]
#              vmpm = self.grid_data['vals'][i+self.counts[2]]
#              vmpp = self.grid_data['vals'][i+self.counts[2]+1]
#
#              vpmm = self.grid_data['vals'][i+self.nyz]
#              vpmp = self.grid_data['vals'][i+self.nyz+1]
#              vppm = self.grid_data['vals'][i+self.nyz+self.counts[2]]
#              vppp = self.grid_data['vals'][i+self.nyz+self.counts[2]+1]
#        
#              # Fraction within the box
#              fx = (x - (ix*self.spacing[0]))/self.spacing[0]
#              fy = (y - (iy*self.spacing[1]))/self.spacing[1]
#              fz = (z - (iz*self.spacing[2]))/self.spacing[2]
#              
#              # Fraction ahead
#              ax = 1 - fx
#              ay = 1 - fy
#              az = 1 - fz
#        
#              # Trilinear interpolation for energy
#              vmm = az*vmmm + fz*vmmp
#              vmp = az*vmpm + fz*vmpp
#              vpm = az*vpmm + fz*vpmp
#              vpp = az*vppm + fz*vppp
#              
#              vm = ay*vmm + fy*vmp
#              vp = ay*vpm + fy*vpp
#              
#              interpolated = (ax*vm + fx*vp)
#              if interpolated==0.0:
#                continue
#              energy += self.scaling_factor[atom_index]*interpolated**self.inv_power
#
#              if not do_gradients:
#                continue
#            
#              # x coordinate
#              dvdx = -vm + vp
#
#              # y coordinate
#              dvdy = (-vmm + vmp)*ax + (-vpm + vpp)*fx
#              
#              # z coordinate
#              dvdz = ((-vmmm + vmmp)*ay + (-vmpm + vmpp)*fy)*ax + ((-vpmm + vpmp)*ay + (-vppm + vppp)*fy)*fx
#            
#              gradients[atom_index] = self.scaling_factor[atom_index]*self.inv_power*interpolated**self.inv_power_m1*Vector(dvdx/self.spacing[0],dvdy/self.spacing[1],dvdz/self.spacing[2])
#            else:
#              raise Exception('Coordinate (%f,%f,%f) outside of grid!'%(x,y,z))
#
#          results = {}
#          results['energy'] = energy*self.strength
#          results['virial'] = -energy*self.strength # Not sure
#          
#          if do_gradients:
#              results['gradients'] = gradients*self.strength
#
#          if do_force_constants:
#              # The force constants are zero -> nothing to calculate.
#              results['force_constants'] = SymmetricPairTensor(self.universe)
#          else:
#              force_constants = None            
#          return results
#

#
# Each force field term requires two classes. One represents the
# abstract force field (ElectricField in this example), it knows
# nothing about any molecular system to which it might be applied.
# The other one (ElectricFieldTerm) stores all the parameters required
# for energy evaluation, and provides the actual code.
#
# When a force field is evaluated for the first time for a given universe,
# the method evaluatorTerms() of the force field is called. It creates
# the EnergyTerm objects. The method evaluate() of the energy term object
# is called for every energy evaluation. Consequently, as much of the
# computation as possible should be done outside of this method.
#

class TrilinearTransformGridForceField(ForceField):

    """
    Force field based on trilinear interpolation between points on a transformed 3D grid
    """

    def __init__(self, FN, strength, scaling_property, scaling_prefactor=None, inv_power=-2., grid_name='trilinear transformed grid', max_val=-1.0):
        """
        @param strength: the electric field vector
        @type strength: L{float}
        @scaling_property: the name of the atomic property in the database
                           that is used to retrieve the atomic scaling_factor.
                           The default is 'amber_charge', the charge property
                           for Amber94 and Amber99.
        @type scaling_property: C{str}
        @param inv_power: the inverse of the power by which grid points are transformed
        @grid_name: a name for the grid
        @max_val: the maximum allowed value for a point on the grid.
                  A negative value means that there is no max.
        """
        # Store arguments that recreate the force field from a pickled
        # universe or from a trajectory.
        self.arguments = (FN, strength, scaling_property, inv_power, grid_name, max_val)
        # Initialize the ForceField class, giving a name to this one.
        ForceField.__init__(self, grid_name)
        # Store the parameters for later use.
        self.FN = FN
        self.strength = strength
        self.scaling_property = scaling_property
        self.inv_power = float(inv_power)
        self.grid_name = grid_name
        self.max_val = max_val
  
        import AlGDock.IO
        IO_Grid = AlGDock.IO.Grid()
        self.grid_data = IO_Grid.read(self.FN, multiplier=0.1)
        if not (self.grid_data['origin']==0.0).all():
          raise Exception('Trilinear grid origin in %s not at (0, 0, 0)!'%self.FN)

        # Make sure all grid values are positive
        n_positive = sum(self.grid_data['vals']>0)
        n_negative = sum(self.grid_data['vals']<0)
        if n_positive>0 and n_negative>0:
          raise Exception('All of the grid points do not have the same sign')
        if n_negative>0:
          self.grid_data['vals'] = -1*self.grid_data['vals']
        if scaling_prefactor is not None:
          self.scaling_prefactor = scaling_prefactor
        else:
          self.scaling_prefactor = -1. if n_negative>0 else 1.

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
                self.grid_name+' inverse power': self.inv_power,
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
        return [TrilinearTransformGridTerm(universe, \
          self.grid_data['spacing'], self.grid_data['counts'], \
          self.grid_data['vals'], \
          self.strength, scaling_factor, self.inv_power, \
          self.grid_name, self.max_val)]
