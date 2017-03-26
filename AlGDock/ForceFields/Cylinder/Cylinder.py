# A flat-bottom harmonic potential on the center of mass
# that keeps the system within a cylinder
# with principal axis along the direction (0,0,1).

from MMTK.ForceFields.ForceField import ForceField
from MMTK_cylinder import CylinderTerm

import numpy as N

# Exactly the same as the pure Python version
class CylinderForceField(ForceField):

    """
    Flat-bottom harmonic potential for a cylinder
    """

    def __init__(self, origin, direction, max_Z, max_R, name='cylinder'):
        """
        @param origin: the origin of the principal axis of the cylinder
        @type origin: {numpy.array}
        @param direction: the direction of the principal axis of the cylinder
        @type direction: {numpy.array}
        @param max_Z: the maximal value of the projection along the principal axis
        @type max_Z: C{float}
        @param max_R: the maximal value orthogonal to the principal axis
        @type max_R: C{float}
        """
        # Store arguments that recreate the force field from a pickled
        # universe or from a trajectory.
        self.arguments = (origin, direction, max_Z, max_R)
        # Initialize the ForceField class, giving a name to this one.
        ForceField.__init__(self, name)
        # Store the parameters for later use.
        self.origin = origin
        self.direction = direction
        self.max_Z = max_Z
        self.max_R = max_R
        self.name = name
        if not ((self.direction[0] == 0.0) and
                (self.direction[1] == 0.0) and
                (self.direction[2] == 1.0)):
          raise Exception("For efficiency, principal axis must be along (0,0,1)")
        # Calculate the cylinder volume
        self.volume = N.pi*(self.max_R*self.max_R)*(self.max_Z - self.origin[2])

    # The following method is called by the energy evaluation engine
    # to inquire if this force field term has all the parameters it
    # requires. This is necessary for interdependent force field
    # terms. In our case, we just say "yes" immediately.
    def ready(self, global_data):
        return True
  
    # The following method is returns a dictionary of parameters for
    # the force field
    def evaluatorParameters(self, universe, subset1, subset2, global_data):
        return {self.name+' origin': self.center,
                self.name+' max Z': self.max_Z,
                self.name+' max R': self.max_R}

    # The following method is called by the energy evaluation engine
    # to obtain a list of the low-level evaluator objects (the C routines)
    # that handle the calculations.
    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        # The subset evaluation mechanism does not make much sense for
        # this force field, so we just signal an error if someone
        # uses it by accident.
        if subset1 is not None or subset2 is not None:
            raise ValueError("sorry, no subsets here")
        # Here we pass all the parameters to the code
        # that handles energy calculations.
        return [CylinderTerm(universe,
                  self.origin, self.direction, self.max_Z, self.max_R,
                  self.name)]

    def randomPoint(self):
      """
      Returns a random point within the cylinder
      """
      z = N.random.uniform()
      (x,y) = self._randomPointInCircle()
      return (x*self.max_R + self.origin[0],
              y*self.max_R + self.origin[1],
              z*(self.max_Z - self.origin[2]) + self.origin[2])
    
    def _randomPointInCircle(self):
      """
      Returns a random point within a unit circle
      """
      r2 = 2
      while r2 > 1:
        (x,y) = N.random.uniform(size=2)
        r2 = x*x + y*y
      return (x,y)
