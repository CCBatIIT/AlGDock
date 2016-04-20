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
    self.rings = sorted(_join_sets(unique_paths), key=lambda r:len(r))

    # Rigid bodies also include terminal atoms adjacent to ring atoms
    # TO DO: include planar systems
    rigid_bodies = []
    for index in range(len(self.rings)):
      attached_terminal_atoms = set(\
        [a for a in self.molecule.atoms \
          if len(a.bondedTo())==1 and a.bondedTo()[0] in self.rings[index]])
      rigid_bodies.append(self.rings[index].union(attached_terminal_atoms))
    self.rigid_bodies = sorted(_join_sets(rigid_bodies), key=lambda b:len(b))

    # Choose initial atom
    if len(self.rings)>0 and len(attached_terminal_atoms)>0:
      # heaviest terminal atom attached to the largest ring
      attached_terminal_atoms = sorted(list(attached_terminal_atoms), \
         key=lambda atom:atom.mass())
      self._converter_setup(
        initial_atom=attached_terminal_atoms[-1])
    else:
      self._converter_setup()
      
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

  def setOccupancyTo(self, property='rings'):
    for atom in self.molecule.atoms:
      atom.occupancy = 0
    if property=='rings' or property=='rigid_bodies':
      for index in range(len(self.rings)):
        for atom in getattr(self, property)[index]:
          atom.occupancy = index+1
    elif property=='soft_torsions':
      for index in self._softTorsionInd:
        for atom in self._torsionIndL[index]:
          atom.occupancy = index+1

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

    self = identifier(universe, molecule)
    print 'There are %d unique rings'%len(self.rings)
    self.setOccupancyTo('soft_torsions')

    # This rotates the last primary torsion
    BAT = self.BAT(extended=True)

    from random import randrange
    torsion_ind = self._softTorsionInd[randrange(len(self._softTorsionInd))]
    BAT_ind = len(BAT)-self.ntorsions+torsion_ind
    confs = []
    for torsion_offset in np.linspace(0,2*np.pi):
      BAT_n = [BAT[ind] if ind!=BAT_ind else BAT[ind] + torsion_offset \
        for ind in range(len(BAT))]
      self.Cartesian(BAT_n)
      confs.append(np.copy(self.universe.configuration().array))

    import AlGDock.IO
    IO_dcd = AlGDock.IO.dcd(molecule)
    IO_dcd.write('rotation.dcd', confs)
    self.showMolecule(dcdFN='rotation.dcd', colorBy='Occupancy')
    os.remove('rotation.dcd')


