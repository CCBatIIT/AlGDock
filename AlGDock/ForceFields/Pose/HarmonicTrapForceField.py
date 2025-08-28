__docformat__ = 'restructuredtext'

from MMTK.ChemicalObjects import isChemicalObject
from MMTK.Collections import isCollection
from MMTK.ForceFields.ForceField import ForceField
from MMTK import Utility
from MMTK_forcefield import HarmonicDistanceTerm, HarmonicAngleTerm, CosineDihedralTerm
from HarmonicCMTrapTerm import HarmonicCMTrapTerm
from Scientific import N

class HarmonicTrapForceField(ForceField):

    """
    Harmonic potential with respect to a fixed point in space
    """

    def __init__(self, obj, center, force_constant):
        """
        :param obj: the object on whose center of mass the force field acts
        :type obj: :class:`~MMTK.Collections.GroupOfAtoms`
        :param center: the point to which the atom is attached by
                        the harmonic potential
        :type center: Scientific.Geometry.Vector
        :param force_constant: the force constant of the harmonic potential
        :type force_constant: float
        """
        if isChemicalObject(obj):
            obj = obj.atomList()
        self.atom_indices = self.getAtomParameterIndices(obj)
        self.arguments = (self.atom_indices, center, force_constant)
        ForceField.__init__(self, 'harmonic_trap')
        self.center = center
        self.force_constant = force_constant

    def ready(self, global_data):
        return True

    def evaluatorParameters(self, universe, subset1, subset2, global_data):
        if universe.is_periodic and len(self.atom_indices) > 1:
            raise ValueError("Center-of-mass restraints not implemented"
                             " for periodic universes")
        for subset in [subset1, subset2]:
            if subset is not None:
                subset = set(a.index for a in subset.atomIterator())
                diff = set(self.atom_indices).difference(subset)
                if diff:
                    if len(diff) == len(self.atom_indices):
                        # The subset doesn't contain the restrained atoms.
                        return {'harmonic_trap_cm': []}
                    else:
                        # The subset contains some but not all of the
                        # restrained atoms.
                        raise ValueError("Restrained atoms partially "
                                         "in a subset")
        return {'harmonic_trap_cm':
                [(self.atom_indices, self.center, self.force_constant)]}

    def evaluatorTerms(self, universe, subset1, subset2, global_data):
        params = self.evaluatorParameters(universe, subset1, subset2,
                                          global_data)['harmonic_trap_cm']
        return [HarmonicCMTrapTerm(universe,
                                   N.array(p[0]),
                                   universe.masses().array,
                                   p[1],
                                   p[2])
                for p in params]