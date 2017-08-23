#!/usr/bin/env python

import MMTK
from AlGDock.BAT import *

def _in_ring(ancestors):
  """
  Returns a list of atoms in the same circular path as ancestors,
  which can be a single atom or list of atoms.
  """
  if isinstance(ancestors,MMTK.ChemicalObjects.Atom):
    ancestors = [ancestors]
  in_ring = []
  for child in [a for a in ancestors[-1].bondedTo() if a not in ancestors[1:]]:
    if child==ancestors[0]:
      if len(ancestors)>2:
        return ancestors
      else:
        return []
    in_ring += _in_ring(ancestors+[child])
    if len(in_ring)>0:
      return in_ring
  return []

def _join_sets(sets):
  """
  Join together sets that have intersecting elements
  :param sets: a list of sets
  """
  while not _sets_are_unique(sets):
    output_sets = []
    for new_set in sets:
      joined = False
      for old_set in output_sets:
        if len(new_set.intersection(old_set))>0:
          output_sets.remove(old_set)
          joined_set = old_set.union(new_set)
          output_sets.append(joined_set)
          joined = True
      if not joined:
        output_sets.append(new_set)
    sets = output_sets
  return sets

def _sets_are_unique(sets):
  """
  Returns True if each set has no intersecting elements with every other set in the list
  :param sets: a list of sets
  """
  nsets = len(sets)
  for a in range(nsets):
    for b in range(a+1,nsets):
      if len(sets[a].intersection(sets[b]))>0:
        return False
  return True

class identifier(converter):
  """
  Identifies rings, rigid bodies, soft torsions
  """
  def __init__(self, universe, molecule):
    self.universe = universe
    self.molecule = molecule
    self.natoms = universe.numberOfAtoms()

    # Identify all unique rings
    unique_paths = []
    for ring in [set(_in_ring(a)) for a in self.molecule.atoms]:
      if ring!=set() and (not ring in unique_paths):
        unique_paths.append(ring)
    self.rings = _join_sets(unique_paths)
    
    # Rigid bodies also include terminal atoms adjacent to ring atoms
    # TO DO: include planar systems
    rigid_bodies = []
    for index in range(len(self.rings)):
      attached_terminal_atoms = set(\
        [a for a in self.molecule.atoms \
          if len(a.bondedTo())==1 and a.bondedTo()[0] in self.rings[index]])
      rigid_bodies.append(self.rings[index].union(attached_terminal_atoms))
    self.rigid_bodies = _join_sets(rigid_bodies)
    
    # Choose initial atom, 
    # preferably the heaviest terminal atom attached to the largest ring
    initial_atom = None
    ring_lists = [sorted(list(ring), key=lambda a:a.fullName()) \
      for ring in self.rings]
    ordered_rings = sorted(sorted(ring_lists, \
      key=lambda r:repr(r)), key=lambda r:len(r), reverse=True)
    for ring in ordered_rings:
      attached_terminal_atoms = sorted(sorted([a for a in self.molecule.atoms \
        if len(a.bondedTo())==1 and a.bondedTo()[0] in ring], \
        key=lambda a:a.fullName()), key=lambda a:a.mass())
      if len(attached_terminal_atoms)>0:
        initial_atom = attached_terminal_atoms[-1]
        break
    self._converter_setup(initial_atom=initial_atom)
    
    # Select soft torsions, primary torsions that are not in the same rigid body
    softTorsionInd = []
    for torsion_ind in set(self._firstTorsionTInd):
      body_id = [-1, -1, -1]
      for atom_ind in [1,2]:
        in_bodies = [self._torsionIndL[torsion_ind][atom_ind] in rigid_body \
          for rigid_body in rigid_bodies]
        if True in in_bodies:
          body_id[atom_ind] = in_bodies.index(True)
        else:
          body_id[atom_ind] = -atom_ind
      if body_id[1] != body_id[2]:
        softTorsionInd.append(torsion_ind)
    self._softTorsionInd = softTorsionInd
    
  def poseInp(self):
    """
    Generates input for the pose force field
    """
    # Restraints on external degrees of freedom
    XYZ = self.universe.configuration().array

    p1 = XYZ[self.rootInd[0]]
    p2 = XYZ[self.rootInd[1]]
    p3 = XYZ[self.rootInd[2]]
    
    ExternalRestraintSpecs = self.rootInd + \
      list(self.extended_coordinates(p1,p2,p3))
    
    # Restraints on internal torsions
    TorsionRestraintSpecs = []
    for ind in self._softTorsionInd:
      t = self._torsionIndL[ind]
      TorsionRestraintSpecs.append(self._torsionIndL[ind]
        + [dihedral(XYZ[t[0]],XYZ[t[1]],XYZ[t[2]],XYZ[t[3]])])

    return [TorsionRestraintSpecs, ExternalRestraintSpecs]

  def setOccupancyTo(self, property='rings'):
    for atom in self.molecule.atoms:
      atom.occupancy = 0
    if property=='rings' or property=='rigid_bodies':
      for index in range(len(self.rings)):
        for atom in getattr(self, property)[index]:
          atom.occupancy = index+1
    elif property=='first_torsions':
      for index in self._firstTorsionTInd:
        for atomInd in self._torsionIndL[index]:
          self.molecule.atoms[atomInd].occupancy = index+1
    elif property=='soft_torsions':
      for index in self._softTorsionInd:
        for atomInd in self._torsionIndL[index]:
          self.molecule.atoms[atomInd].occupancy = index+1

########
# MAIN #
########

if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(
    description='Bond-Angle-Torsion converter')
  parser.add_argument('--database',
    help='MMTK database that describes the molecule of interest',
    default='ligand.db')
  args = parser.parse_args()

  import os.path
  if args.database=='all':
    import glob
    dbFNs = glob.glob('*/*.db')
  elif args.database=='random':
    import glob
    from random import randrange
    dbFNs = glob.glob('*/*.db')
    dbFNs = [dbFNs[randrange(len(dbFNs))]]
  elif os.path.isfile(args.database):
    dbFNs = [args.database]
  else:
    raise Exception('Database file %s not found'%args.database)

  dirname = os.path.dirname(os.path.abspath(dbFNs[0]))
  for FN in dbFNs:
    print 'Loading', FN
    dbFN = os.path.abspath(FN)
    if os.path.dirname(dbFN)!=dirname:
      raise Exception('Cannot change ligand directory in MMTK')
    MMTK.Database.molecule_types.directory = os.path.dirname(dbFN)
    universe = MMTK.Universe.InfiniteUniverse()
    molecule = MMTK.Molecule(os.path.basename(dbFN))
    universe.addObject(molecule)

    original_xyz = np.copy(universe.configuration().array)
    self = identifier(universe, molecule)

    print 'There are %d unique rings'%len(self.rings)
    self.setOccupancyTo('soft_torsions')

    # This tests a conversion to BAT coordinates and back
    BAT = self.BAT(original_xyz, extended=True)
    new_xyz = self.Cartesian(BAT)
    print sum(sum(new_xyz - original_xyz))

    # This rotates a random primary torsion
    from random import randrange
    softTorsionInd = self._softTorsionInd[randrange(len(self._softTorsionInd))]
    BAT_ind = np.array(range(6+5,self.natoms*3,3))[softTorsionInd]
    confs = []
    for torsion_offset in np.linspace(0,2*np.pi):
      BAT_n = [BAT[ind] if ind!=BAT_ind else BAT[ind] + torsion_offset \
        for ind in range(len(BAT))]
      XYZ = self.Cartesian(BAT_n)
      confs.append(XYZ)

    import AlGDock.IO
    IO_dcd = AlGDock.IO.dcd(molecule)
    IO_dcd.write('rotation.dcd', confs)
    self.showMolecule(dcdFN='rotation.dcd', colorBy='Occupancy')
    os.remove('rotation.dcd')


