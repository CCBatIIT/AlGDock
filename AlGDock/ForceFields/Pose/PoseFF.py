# Force field that makes use of soft internal torsion angles and 
# additional 6 external degrees of freedom to keep the molecule in a 
# particular pose. It uses flat bottom harmonic potentials
# Laurentiu Spiridon 2/11/2016

from MMTK.ForceFields.ForceField import ForceField, EnergyTerm
from MMTK_pose import PoseDihedralTerm, PoseExtDistTerm, PoseExtAnglTerm, PoseExtDiheTerm, PoseExtDihe2Term
from MMTK import ParticleVector, SymmetricPairTensor
from Scientific.Geometry import delta
from Scientific import N

from collections import OrderedDict

class InternalRestraintForceField(ForceField):
  """
  Flat bottom harmonic dihedral for internal degrees of freedom
  """

  def __init__(self, torsions, \
               hwidth = 0.7853981633974483, \
               k = 1000.0):
    """
    Internal torsional dofs:
    @param input: [4 atom indices, gamma] x internal dofs
    @type atom:   L{MMTK.ChemicalObjects.Atom.index,
                    MMTK.ChemicalObjects.Atom.index,
                    MMTK.ChemicalObjects.Atom.index,
                    MMTK.ChemicalObjects.Atom.index,
                    dobule} x # of internal dofs}
    """
    
    # Initialize the ForceField class, giving a name to this one.
    ForceField.__init__(self, 'internal dihedral restraint')

    # Store the parameters for later use
    self.params = OrderedDict()
    for key in [\
      'torsions', 'hwidth', 'k']:
      self.params[key] = locals()[key]

  def set_hwidth(self, hwidth):
    self.params['hwidth'] = hwidth

  def set_k(self, k):
    self.params['k'] = k

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

    # intDiheIs[] = 4 atom indices
    intDiheIs = N.array([t[:4] for t in self.params['torsions']])

    # intDiheParams[0] = n: periodicity
    # intDiheParams[1] = gamma: reference angle
    # intDiheParams[2] = b: flat bottom range (it gets halfed below)
    # intDiheParams[3] = k: force constant
    intDiheParams = N.array([[0, t[4], self.params['hwidth'], self.params['k']] \
      for t in self.params['torsions']])

    return [PoseDihedralTerm(universe._spec, intDiheIs, intDiheParams)]

class ExternalRestraintForceField(ForceField):
  """
  Flat bottom harmonic dihedral and harmonic distance restraint for external degrees of freedom
  """

  def __init__(self, ind1, ind2, ind3, X1, Y1, Z1, phi, theta, omega, \
               hwidth_spatial = 0.05, \
               k_spatial = 1000.0, \
               hwidth_angular = 0.7853981633974483, \
               k_angular = 1000.0):
    """
    @param input: [index1, index2, index3, X1, Y1, Z1, phi, theta, omega]
    ind1 - index of atom 1
    ind2 - index of atom 2
    ind3 - index of atom 3
    The next 6 variables are the center of the flat-bottom harmonic bias
    X1 - x coordinate of atom 1
    X2 - x coordinate of atom 2
    X3 - x coordinate of atom 3
    phi, theta, omega - angular degrees of freedom, in radians
    hwidth_spatial - the half width of the spatial flat-bottom harmonic terms (nm)
    k_spatial - the spring constant on the spatial flat-bottom harmonic terms (kJ/nm)
    hwidth_angular - the half width of the angular flat-bottom harmonic terms (radians)
    k_angular - the spring constant on the angular flat-bottom harmonic terms (kJ/nm)
    """
    
    # Initialize the ForceField class, giving a name to this one.
    ForceField.__init__(self, 'external restraint')

    # Store the parameters for later use.
    self.params = OrderedDict()
    for key in [\
      'ind1', 'ind2', 'ind3', 'X1', 'Y1', 'Z1', 'phi', 'theta', 'omega', \
      'hwidth_spatial', 'k_spatial', 'hwidth_angular', 'k_angular']:
      self.params[key] = locals()[key]

  def set_hwidth_spatial(self, hwidth_spatial):
    self.params['hwidth_spatial'] = hwidth_spatial

  def set_k_spatial(self, k_spatial):
    self.params['k_spatial'] = k_spatial

  def set_hwidth_angular(self, hwidth_angular):
    self.params['hwidth_angular'] = hwidth_angular
  
  def set_k_angular(self, k_angular):
    self.params['k_angular'] = k_angular
  
  def get_reference_external_BAT(self):
    return N.array([self.params[key] \
      for key in ['X1','Y1','Z1','phi','theta','omega']])
  
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

    offset = 0.3

    ind1 = self.params['ind1']
    ind2 = self.params['ind2']
    ind3 = self.params['ind3']
    phi  = self.params['phi']
    theta = self.params['theta']
    omega = self.params['omega']
    hwidth_spatial = self.params['hwidth_spatial']
    k_spatial = self.params['k_spatial']
    hwidth_angular = self.params['hwidth_angular']
    k_angular = self.params['k_angular']
    
    # extXYZ[]: cartesian coordinates of an external reference vertex
    # ext*Is[]: indices of the picked atoms
    # ext[Dist,Angl]Params[]: ref value, force K, bottom width
    # extDiheParams[]: same shape as intDiheParams

    extXYZ  = N.array([[self.params['X1'], self.params['Y1'], self.params['Z1']]])
    extDistIs  = N.array([[-1, ind1]]);
    extDistParams  = N.array([[0.0, k_spatial, hwidth_spatial]]);

    PhiDummy1 = N.array([[offset, 0.0, -offset]]) + extXYZ
    PhiDummy2 = N.array([[0.0, 0.0, -offset]]) + extXYZ
    extDihe2Is  = N.array([[-2, -1, ind1, ind2]]);
    extDihe2Params  = N.array([[0, phi, hwidth_angular, k_angular]]);

    ThetaDummy = N.array([[0.0, 0.0, offset]]) + extXYZ
    extAnglIs  = N.array([[-1, ind1, ind2]]);
    extAnglParams  = N.array([[theta, k_angular, hwidth_angular]]);

    OmegaDummy = N.array([[0.0, 0.0, -1.0]]) + extXYZ
    extDiheIs  = N.array([[-1, ind1, ind2, ind3]]);
    extDiheParams  = N.array([[0, omega, hwidth_angular, k_angular]]);

    # Here we pass all the parameters to the code
    # that handles energy calculations.

    return [PoseExtDistTerm(universe._spec, extDistIs, extDistParams, extXYZ), \
            PoseExtDihe2Term(universe._spec, extDihe2Is, extDihe2Params, PhiDummy1, PhiDummy2), \
            PoseExtAnglTerm(universe._spec, extAnglIs, extAnglParams, ThetaDummy), \
            PoseExtDiheTerm(universe._spec, extDiheIs, extDiheParams, OmegaDummy)]
