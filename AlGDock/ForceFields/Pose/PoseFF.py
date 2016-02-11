# Force field that makes use of soft internal torsion angles and 
# additional 6 external degrees of freedom to keep the molecule in a 
# particular pose. It uses flat bottom harmonic potentials
# Laurentiu Spiridon 2/11/2016

from MMTK.ForceFields.ForceField import ForceField, EnergyTerm
from MMTK_pose import PoseDihedralTerm, PoseExtDistTerm, PoseExtAnglTerm, PoseExtDiheTerm
from MMTK import ParticleVector, SymmetricPairTensor
from Scientific.Geometry import delta
from Scientific import N

class PoseForceField(ForceField):

    """
    Flat bottom harmonic dihedral and harmonic distance restraint
    """

    def __init__(self, input_table):
        """
        First 2 entries for external dofs:
        [index1, index2, index3, X, Y, Z, bbot, abot]
        [b0, kb, theta0, ka, n, gamma, dbot, kd]
        Internal torsional dofs:
        @param input_table: [index1, index2, index3, X, Y, Z, bbot, abot]
                            [b0, kb, theta0, ka, n, gamma, dbot, kd]
                            [4 atom indices, periodicity, gamma, bottom, force K] x internal dofs
        @type atom: L{L{MMTK.ChemicalObjects.Atom.index,
                        MMTK.ChemicalObjects.Atom.index,
                        MMTK.ChemicalObjects.Atom.index,
                        double, double, double, double, double},
                      L{double, double, double, double, 
                        int, double, double, double},
                      L{MMTK.ChemicalObjects.Atom.index,
                        MMTK.ChemicalObjects.Atom.index,
                        MMTK.ChemicalObjects.Atom.index,
                        MMTK.ChemicalObjects.Atom.index,
                        inp, dobule, double, double} x #of internal dofs}
        """

        # Internal dofs arguments
        self.externalArgs = input_table[:2]
        self.internalArgs = input_table[2:]
        # Halve the flat bottoms
        self.externalArgs[0][6] /= 2
        self.externalArgs[0][7] /= 2
        for ia in self.internalArgs:
          ia[6] /= 2
        # Initialize the ForceField class, giving a name to this one.
        ForceField.__init__(self, 'pose')
        # Store the parameters for later use.

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
        # The subset evaluation mechanism does not make much sense for
        # this force field, so we just signal an error if someone
        # uses it by accident.
        if subset1 is not None or subset2 is not None:
            raise ValueError("sorry, no subsets here")
        # Here we pass all the parameters to the code
        # that handles energy calculations.

        # extXYZ[]: cartesian coordinates of an external reference vertex
        # ext*Is[]: indices of the picked atoms
        # ext[Dist,Angl]Params[]: ref value, force K, bottom width
        # extDiheParams[]: same shape as intDiheParams (described below) 
        extXYZ  = N.array([[self.externalArgs[0][3], self.externalArgs[0][4], self.externalArgs[0][5]]]);
 
        extDistIs  = N.array([[-1, self.externalArgs[0][0]]]);
        extDistParams  = N.array([[self.externalArgs[1][0], self.externalArgs[1][1], self.externalArgs[0][6]]]);

        extAnglIs  = N.array([[-1, self.externalArgs[0][0], self.externalArgs[0][1]]]);
        extAnglParams  = N.array([[self.externalArgs[1][2], self.externalArgs[1][3], self.externalArgs[0][7]]]);

        extDiheIs  = N.array([[-1, self.externalArgs[0][0], self.externalArgs[0][1], self.externalArgs[0][2]]]);
        extDiheParams  = N.array([[self.externalArgs[1][4], self.externalArgs[1][5], self.externalArgs[1][6], self.externalArgs[1][7]]]);

        # intDiheIs[] = 4 atom indices
        # intDiheParams[0] = n : periodicity
        # intDiheParams[1] = gamma: reference angle
        # intDiheParams[2] = b: flat bottom range (it gets halfed below)
        # intDiheParams[3] = k: force constant
        intDiheIs    = N.array(map(lambda d: d[:4], self.internalArgs))
        intDiheParams = N.array(map(lambda d: d[4:], self.internalArgs))

        # Debug:
        #print "PoseFF external XYZ:", extXYZ
        #print "PoseFF external distance indices:", extDistIs
        #print "PoseFF external distance parameters:", extDistParams
        #print "PoseFF external angle indices:", extAnglIs
        #print "PoseFF external angle parameters:", extAnglParams
        #print "PoseFF external dihe indices:", extDiheIs
        #print "PoseFF external dihe parameters:", extDiheParams
        #print "PoseFF internal dihedrals indices:", intDiheIs
        #print "PoseFF internal dihedrals parameters: ", intDiheParams

        # -1 index blocks swapi? for force constants in pose.c
        return [PoseExtDistTerm(universe._spec, extDistIs, extDistParams, extXYZ), \
                PoseExtAnglTerm(universe._spec, extAnglIs, extAnglParams, extXYZ), \
                PoseExtDiheTerm(universe._spec, extDiheIs, extDiheParams, extXYZ), \
                PoseDihedralTerm(universe._spec, intDiheIs, intDiheParams)]




