# Force field that makes use of soft internal torsion angles and 
# additional 6 external degrees of freedom to keep the molecule in a 
# particular pose. It uses flat bottom harmonic potentials
# Laurentiu Spiridon 2/11/2016

from MMTK.ForceFields.ForceField import ForceField, EnergyTerm
from MMTK_pose import PoseDihedralTerm, PoseExtDistTerm, PoseExtAnglTerm, PoseExtDiheTerm, PoseExtDihe2Term
from MMTK import ParticleVector, SymmetricPairTensor
from Scientific.Geometry import delta
from Scientific import N

class PoseForceField(ForceField):

    """
    Flat bottom harmonic dihedral and harmonic distance restraint
    """

    def __init__(self, input):
        """
        First 2 entries for external dofs:
        [index1, index2, index3, X, Y, Z, bbot, abot]
        [b0, kb, theta0, ka, n, gamma, dbot, kd]
        Internal torsional dofs:
        @param input: [index1, index2, index3, X, Y, Z, bbot, abot]
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

        self.input = input
        #print "input", self.input
        offset = 0.3
        
        # 0.03490659/2 is 20 degrees (full width)/2 (to make half width)

        self.extXYZ  = N.array([[self.input[0][3], self.input[0][4], self.input[0][5]]])
        self.extDistIs  = N.array([[-1, int(self.input[0][0])]]);
        self.extDistParams  = N.array([[0.0, 1000.0, 0.01/2]]); # Half the flat bot
        print "PoseFF: extXYZ: %.8lf %.8lf %.8lf"%(self.extXYZ[0][0], self.extXYZ[0][1], self.extXYZ[0][2])
        print "PoseFF: extDistIs:", self.extDistIs
        print "PoseFF: extDistParams: %.8lf %.8lf %.8lf\n"%(self.extDistParams[0][0], self.extDistParams[0][1], self.extDistParams[0][2])

        self.PhiDummy1 = N.array([[offset, 0.0, -offset]]) + self.input[0][3:6]
        self.PhiDummy2 = N.array([[0.0, 0.0, -offset]]) + self.input[0][3:6]
        self.extDihe2Is  = N.array([[-2, -1, int(self.input[0][0]), int(self.input[0][1])]]);
        self.extDihe2Params  = N.array([[0, self.input[0][9], 0.03490659/2, 1000.0]]);
        print "PoseFF: PhiDummies: [%.4f %.4f %.4f] [%.4f %.4f %.4f]"\
          %(self.PhiDummy1[0][0], self.PhiDummy1[0][1], self.PhiDummy1[0][2], \
            self.PhiDummy2[0][0], self.PhiDummy2[0][1], self.PhiDummy2[0][2])
        print "PoseFF: extDihe2Is", self.extDihe2Is
        print "PoseFF: extDihe2Params: %.8lf %.8lf %.8lf\n"\
          %(self.extDihe2Params[0][0], self.extDihe2Params[0][1], self.extDihe2Params[0][2])

        self.ThetaDummy = N.array([[0.0, 0.0, offset]]) + self.input[0][3:6]
        self.extAnglIs  = N.array([[-1, int(self.input[0][0]), int(self.input[0][1])]]);
        self.extAnglParams  = N.array([[self.input[0][13], 1000.0 , 0.03490659/2]]); # Half the flat bot
        print "PoseFF: ThetaDummy:", self.ThetaDummy
        print "PoseFF: extAnglIs", self.extAnglIs
        print "PoseFF: extAnglParams: %.8lf %.8lf %.8lf\n"%(self.extAnglParams[0][0], self.extAnglParams[0][1], self.extAnglParams[0][2])

        self.OmegaDummy = N.array([[0.0, 0.0, -1.0]]) + self.input[0][3:6]
        self.extDiheIs  = N.array([[-1, int(self.input[0][0]), int(self.input[0][1]), int(self.input[0][2])]]);
        self.extDiheParams  = N.array([[0, self.input[0][14], 0.03490659/2, 1000.0]]);
        print "PoseFF: OmegaDummy:", self.OmegaDummy
        print "PoseFF: extDiheIs", self.extDiheIs
        print "PoseFF: extDiheParams: %.8lf %.8lf %.8lf"%(self.extDiheParams[0][0], self.extDiheParams[0][1], self.extDiheParams[0][2])

        # Halve the flat bottoms !!!!!!!!!!!!!!!!! 
        #self.input[0][6] /= 2
        #self.input[0][7] /= 2
        #for ia in self.internalArgs:
        #  ia[6] /= 2
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
        if subset1 is not None or subset2 is not None:
            raise ValueError("sorry, no subsets here")
        # Here we pass all the parameters to the code
        # that handles energy calculations.

        # extXYZ[]: cartesian coordinates of an external reference vertex
        # ext*Is[]: indices of the picked atoms
        # ext[Dist,Angl]Params[]: ref value, force K, bottom width
        # extDiheParams[]: same shape as intDiheParams (described below) 



        # intDiheIs[] = 4 atom indices
        # intDiheParams[0] = n : periodicity
        # intDiheParams[1] = gamma: reference angle
        # intDiheParams[2] = b: flat bottom range (it gets halfed below)
        # intDiheParams[3] = k: force constant
        intDiheIs    = N.array(map(lambda d: d[:4], self.input[1]))
        intDiheParams = N.array(map(lambda d: d[4:], self.input[1]))
        print "intDiheIs", intDiheIs
        print "intDiheParams", intDiheParams

        return [PoseExtDistTerm(universe._spec, self.extDistIs, self.extDistParams, self.extXYZ), \
                PoseExtDihe2Term(universe._spec, self.extDihe2Is, self.extDihe2Params, self.PhiDummy1, self.PhiDummy2), \
                PoseExtAnglTerm(universe._spec, self.extAnglIs, self.extAnglParams, self.ThetaDummy), \
                PoseExtDiheTerm(universe._spec, self.extDiheIs, self.extDiheParams, self.OmegaDummy), 
                PoseDihedralTerm(universe._spec, intDiheIs, intDiheParams)
               ]




