# Onufriev-Bashford-Case Generalized Born term
# can be added to any other force field

import numpy as np

from MMTK import ParticleScalar
from MMTK.ForceFields.ForceField import ForceField

class OBCForceField(ForceField):

    """
    Onufriev-Bashford-Case Generalized Born
    """

    def __init__(self,
          prmtopFN=None,
          inv_prmtop_atom_order=None,
          desolvationGridFN=None,
          strength=1.0):
        """
        @param prmtopFN: an AMBER parameter and topology file
        @type strength:  C{str}
        """
        # Initialize the ForceField class, giving a name to this one.
        ForceField.__init__(self, 'OBC')

        # Store arguments that recreate the force field from a pickled
        # universe or from a trajectory.
        self.arguments = (prmtopFN, inv_prmtop_atom_order, \
          desolvationGridFN, strength)

        # Load the desolvation grid
        if desolvationGridFN is not None:
          import AlGDock.IO
          IO_Grid = AlGDock.IO.Grid()
          self.grid_data = IO_Grid.read(desolvationGridFN, multiplier=0.1)
          if not (self.grid_data['origin']==0.0).all():
            raise Exception('Trilinear grid origin in %s not at (0, 0, 0)!'%FN)
          self.useDesolvationGrid = True
        else:
          self.grid_data = {'spacing':np.array([0., 0., 0.]), \
                            'counts':np.array([0, 0, 0]), \
                            'vals':np.array([])}
          self.useDesolvationGrid = False
        
        # Store arguments as class variables
        self.prmtopFN = prmtopFN
        self.inv_prmtop_atom_order = inv_prmtop_atom_order
        self.desolvationGridFN = desolvationGridFN
        self.strength = strength

    def set_strength(self, strength):
      self.strength = strength

    # The following method is called by the energy evaluation engine
    # to inquire if this force field term has all the parameters it
    # requires. This is necessary for interdependent force field
    # terms. In our case, we just say "yes" immediately.
    def ready(self, global_data):
        return True

    # The following method is called by the energy evaluation engine
    # to obtain a list of the low-level evaluator objects (the C routines)
    # that handle the calculations.
    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        # The energy for subsets is defined as consisting only
        # of interactions within that subset, so the contribution
        # of an external field is zero. Therefore we just return
        # an empty list of energy terms.
        if subset1 is not None or subset2 is not None:
            return []

        import numpy as np
        if (self.prmtopFN is not None) and \
           (self.inv_prmtop_atom_order is not None):
          # Get charges, radii, and scale factors from OpenMM
          import simtk.openmm
          import simtk.openmm.app as OpenMM_app
          
          prmtop = OpenMM_app.AmberPrmtopFile(self.prmtopFN)
          OMM_system = prmtop.createSystem(\
            nonbondedMethod=OpenMM_app.CutoffNonPeriodic, \
            nonbondedCutoff=1.5, \
            constraints=None, \
            implicitSolvent=OpenMM_app.OBC2)
          f = OMM_system.getForces()[-2]

          numParticles = f.getNumParticles()
          charges = np.zeros(numParticles)
          atomicRadii = np.zeros(numParticles)
          scaleFactors = np.zeros(numParticles)
          for n in range(numParticles):
            (charge, radius, scaleFactor) = f.getParticleParameters(n)
            charges[n] = charge/simtk.unit.elementary_charge
            atomicRadii[n] = radius/simtk.unit.nanometer
            scaleFactors[n] = scaleFactor

          charges = charges[self.inv_prmtop_atom_order]
          atomicRadii = atomicRadii[self.inv_prmtop_atom_order]
          scaleFactors = scaleFactors[self.inv_prmtop_atom_order]
        else:
          # Get charges, radii, and scale factors from the ligand database (preferred)
          charges_ps = ParticleScalar(universe)
          atomicRadii_ps = ParticleScalar(universe)
          scaleFactors_ps = ParticleScalar(universe)
          for o in universe:
            for a in o.atomList():
              charges_ps[a] = o.getAtomProperty(a, 'amber_charge')
              atomicRadii_ps[a] = o.getAtomProperty(a, 'scaling_factor_BornRadii')
              scaleFactors_ps[a] = o.getAtomProperty(a, 'scaling_factor_BornScreening')
          charges = charges_ps.array
          atomicRadii = atomicRadii_ps.array
          scaleFactors = scaleFactors_ps.array
          numParticles = charges.shape[0]
          
#        import time
#        import os.path
#        import MMTK_OBC
#        OBCpath = MMTK_OBC.__file__
#        print """
#        in {0}
#        last modified {1}
#            """.format(OBCpath, time.ctime(os.path.getmtime(OBCpath)))

        # Here we pass all the parameters as "simple" data types to
        # the C code that handles energy calculations.
        if self.useDesolvationGrid:
          # With desolvation grid
          from MMTK_OBC_desolv import OBCDesolvTerm
          return [OBCDesolvTerm(universe._spec, numParticles, self.strength, \
            charges, atomicRadii, scaleFactors, \
            self.grid_data['spacing'], self.grid_data['counts'], \
            self.grid_data['vals'])]
        else:
          # No desolvation grid
          from MMTK_OBC import OBCTerm
          return [OBCTerm(universe._spec, numParticles, self.strength, \
            charges, atomicRadii, scaleFactors)]
