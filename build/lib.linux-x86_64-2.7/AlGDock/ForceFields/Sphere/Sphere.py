# A flat-bottom harmonic potential on the center of mass
# that keeps the system within a sphere

import AlGDock
from MMTK.ForceFields.ForceField import ForceField
from MMTK_sphere import SphereTerm

import numpy as N

# Exactly the same as the pure Python version
class SphereForceField(ForceField):

    """
    Flat-bottom harmonic potential for a sphere
    """

    def __init__(self, center, max_R, name='Sphere'):
        """
        @param center: the center of the principal axis of the sphere
        @type center: {numpy.array}
        @param max_R: the maximal value orthogonal to the principal axis
        @type max_R: C{float}
        """
        # Store arguments that recreate the force field from a pickled
        # universe or from a trajectory.
        self.arguments = (center, max_R)
        # Initialize the ForceField class, giving a name to this one.
        ForceField.__init__(self, name)
        # Store the parameters for later use.
        self.center = center
        self.max_R = max_R
        self.name = name
        # Calculate the sphere volume
        self.volume = 4./3.*N.pi*(self.max_R*self.max_R*self.max_R)

    # The following method is called by the energy evaluation engine
    # to inquire if this force field term has all the parameters it
    # requires. This is necessary for interdependent force field
    # terms. In our case, we just say "yes" immediately.
    def ready(self, global_data):
        return True
  
    # The following method is returns a dictionary of parameters for
    # the force field
    def evaluatorParameters(self, universe, subset1, subset2, global_data):
        return {self.name+' center': self.center,
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
        return [SphereTerm(universe, self.center, self.max_R, self.name)]
  
    def randomPoint(self):
      """
      Returns a random point within the sphere
      """
      (x,y,z) = self._randomPointInSphere()
      return (x*self.max_R + self.center[0],
              y*self.max_R + self.center[1],
              z*self.max_R + self.center[2])

    def _randomPointInSphere(self):
      """
      Returns a random point within a unit circle
      """
      r2 = 2
      while r2 > 1:
        (x,y,z) = N.random.uniform(size=3)
        r2 = x*x + y*y + z*z
      return (x,y,z)

