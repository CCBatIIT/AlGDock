#!/usr/bin/env python

# TODO: Find a sincos function that simultaneously calculates sine and cosine?

import numpy as np
cimport numpy as np
import cython

from libc.math cimport sin
from libc.math cimport cos
from libc.math cimport sqrt
from libc.math cimport acos
from libc.math cimport atan2

ctypedef np.float_t float_t
ctypedef np.int_t int_t

import MMTK

### Vector functions
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline dotp(np.ndarray[np.float_t] v1, v2):
  return (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline norm2(np.ndarray[np.float_t] v):
  return (v[0]*v[0] + v[1]*v[1] + v[2]*v[2])

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline normalize(np.ndarray[np.float_t] v1):
  return v1/sqrt(norm2(v1))

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline cross(np.ndarray[np.float_t] v1, np.ndarray[np.float_t] v2):
  return np.array([v1[1]*v2[2]-v1[2]*v2[1], \
    v1[2]*v2[0]-v1[0]*v2[2], \
    v1[0]*v2[1]-v1[1]*v2[0]])

### BAT coordinate measurement
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline distance(np.ndarray[np.float_t] p1,np.ndarray[np.float_t] p2):
  cdef np.ndarray[np.float_t] v1 = p2 - p1
  return sqrt(norm2(v1))

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline angle(np.ndarray[np.float_t] p1, np.ndarray[np.float_t] p2, \
    np.ndarray[np.float_t] p3):
  cdef np.ndarray[np.float_t] v1 = p2 - p1
  cdef np.ndarray[np.float_t] v2 = p2 - p3
  return acos(max(-1.,min(1.,dotp(v1,v2)/\
    sqrt(norm2(v1)*norm2(v2)))))

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef dihedral(np.ndarray[np.float_t] p1, np.ndarray[np.float_t] p2, \
    np.ndarray[np.float_t] p3, np.ndarray[np.float_t] p4):
  cdef np.ndarray[np.float_t] v1, v2, v3, a, b
  cdef float c, s
  v1 = p2 - p1
  v2 = p2 - p3
  v3 = p3 - p4
  a = normalize(cross(v1,v2))
  b = normalize(cross(v3,v2))
  c = dotp(a,b)
  s = dotp(cross(b,a),normalize(v2))
  return atan2(s,c)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef BAT4(np.ndarray[np.float_t] p1, np.ndarray[np.float_t] p2, \
    np.ndarray[np.float_t] p3, np.ndarray[np.float_t] p4):
  cdef np.ndarray[np.float_t] v1, v2, v3, a, b
  cdef float c, s, norm_v1_2, norm_v2_2
  v1 = p2 - p1
  v2 = p2 - p3
  v3 = p3 - p4
  a = normalize(cross(v1,v2))
  b = normalize(cross(v3,v2))
  c = dotp(a,b)
  norm_v1_2 = norm2(v1)
  norm_v2_2 = norm2(v2)
  s = dotp(cross(b,a),v2/sqrt(norm_v2_2))
  return (sqrt(norm_v1_2), acos(max(-1.,min(1.,dotp(v1,v2)/\
    sqrt(norm_v1_2*norm_v2_2)))), atan2(s,c))

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef extended_coordinates(np.ndarray[np.float_t] p1, \
    np.ndarray[np.float_t] p2, np.ndarray[np.float_t] p3):
  # The rotation axis is a normalized vector pointing from atom 0 to 1
  # It is described in two degrees of freedom by the polar angle and azimuth
  cdef np.ndarray[np.float_t] e = normalize(p2-p1)
  cdef float phi = atan2(e[1],e[0]) # Polar angle
  cdef float theta = acos(e[2]) # Azimuthal angle
  # Rotation to the z axis
  cdef float cp = cos(phi)
  cdef float sp = sin(phi)
  cdef float ct = cos(theta)
  cdef float st = sin(theta)
  cdef np.ndarray[np.float_t, ndim=2] Rz = np.array([[cp*ct,ct*sp,-st],[-sp,cp,0],[cp*st,sp*st,ct]])
  cdef np.ndarray[np.float_t] pos2 = Rz.dot(np.array(p3-p1))
  # Angle about the rotation axis
  cdef float omega = atan2(pos2[1],pos2[0])
  cdef np.ndarray[np.float_t] coords = np.zeros((6,))
  coords[:3] = p1
  coords[3] = phi
  coords[4] = theta
  coords[5] = omega
  return coords

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef get_XYZ1(np.ndarray[np.float_t] p2, \
    np.ndarray[np.float_t] p3, np.ndarray[np.float_t] p4,
    float bond, float angle, float torsion):

  cdef float s, c
  cdef np.ndarray[np.float_t] n23, v21

  n23 = normalize(p3-p2)

  s = sin(angle)
  c = cos(angle)
  v21 = (bond*c)*n23 - (bond*s)*cross(normalize(cross(p4-p3,n23)),n23)

  s = sin(torsion)
  c = cos(torsion)
  return p2 - cross(n23,v21)*s + dotp(n23,v21)*n23*(1.0-c) + v21*c

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef rotate123(np.ndarray[np.float_t] p1, \
    np.ndarray[np.float_t] p2, np.ndarray[np.float_t] p3, \
    np.ndarray[np.float_t] origin, \
    float phi, float theta, float omega):
  # Rotate the third atom by the appropriate value
  cdef float co = cos(omega)
  cdef float so = sin(omega)
  cdef np.ndarray[np.float_t, ndim=2] Romega = np.array([[co, -so, 0],[so, co, 0], [0, 0, 1]])
  p3 = Romega.dot(p3)
  # Rotate the second two atoms to point in the right direction
  cdef float cp = cos(phi)
  cdef float sp = sin(phi)
  cdef float ct = cos(theta)
  cdef float st = sin(theta)
  cdef np.ndarray[np.float_t, ndim=2] Re = np.array([[cp*ct,-sp,cp*st],[ct*sp,cp,sp*st],[-st,0,ct]])
  p2 = Re.dot(p2)
  p3 = Re.dot(p3)
  # Translate the first three atoms by the origin
  p1 += origin
  p2 += origin
  p3 += origin
  return (p1,p2,p3)

# Main converter class
cdef class converter:
  """
  Interconverts Cartesian and Bond-Angle-Torsion coordinates
  """
  cdef readonly int natoms, ntorsions
  cdef object rootInd, _torsionIndL
  cdef object _firstTorsionTInd

  def __init__(self, universe, molecule, initial_atom=None):
    self.natoms = universe.numberOfAtoms()
    self._converter_setup(molecule.atoms, initial_atom)

  def _converter_setup(self, atoms, initial_atom=None):
    atom_mass = lambda atom:atom.mass()
    terminal_atoms = sorted(\
      [a for a in atoms if len(a.bondedTo())==1], \
      key=atom_mass)
    if (initial_atom is None):
      # Select the heaviest root atoms from the heaviest terminal atom
      root = [terminal_atoms[-1]]
    else:
      if initial_atom in terminal_atoms:
        root = [initial_atom]
      else:
        raise Exception('Initial atom is not a terminal atom')
    root.append(sorted(root[0].bondedTo(),key=atom_mass)[-1])
    root.append(sorted([a for a in root[-1].bondedTo() \
      if (a not in root) and (a not in terminal_atoms)],key=atom_mass)[-1])

    def _find_dihedral(selected):
      """
      Finds a dihedral angle adjacent to the selected atoms that includes a new atom
      :param selected: a list of atoms that have already been selected
      :returns: a list of atoms that includes the new atom and its neighboring selected atoms
      """
      atom_mass = lambda atom:atom.mass()
      # Loop over possible nearest neighbors
      for a2 in selected:
        # Find the new atom
        for a1 in sorted([a for a in a2.bondedTo() \
            if a not in selected],key=atom_mass,reverse=True):
          # Find the third atom
          for a3 in sorted([a for a in a2.bondedTo() \
              if (a in selected) and (a!=a1)],key=atom_mass,reverse=True):
            # Find the last atom
            for a4 in sorted([a for a in a3.bondedTo() \
                if (a in selected) and (a!=a2)],key=atom_mass,reverse=True):
              return (a1,a2,a3,a4)
      print 'Selected atoms:', selected
      raise Exception('No new dihedral angle found!')

    # Construct a list of torsion angles
    torsionL = []
    selected = [a for a in root]
    while len(selected)<self.natoms:
      (a1,a2,a3,a4) = _find_dihedral(selected)
      torsionL.append((a1,a2,a3,a4))
      selected.append(a1)
    # If _firstTorsionTInd is not equal to the list index,
    # then the dihedrals will likely be correlated and it is more appropriate
    # to use a relative phase angle
    prior_atoms = [(a2,a3,a4) for (a1,a2,a3,a4) in torsionL]

    self.rootInd = [r.index for r in root]
    self._torsionIndL = [[a.index for a in tset] for tset in torsionL]
    self._firstTorsionTInd = [prior_atoms.index(prior_atoms[n]) \
      for n in range(len(prior_atoms))]
    self.ntorsions = self.natoms-3

  def getFirstTorsionInds(self, extended):
    offset = 6 if extended else 0
    torsionInds = np.array(range(offset+5,self.natoms*3,3))
    primaryTorsions = sorted(list(set(self._firstTorsionTInd)))
    return list(torsionInds[primaryTorsions])

  @cython.boundscheck(False)
  @cython.wraparound(False)
  @cython.cdivision(True)
  def BAT(self, np.ndarray[np.float_t, ndim=2] XYZ, extended=False):
    """
    Conversion from Cartesian to Bond-Angle-Torsion coordinates
    :param extended: whether to include external coordinates or not
    :param Cartesian: Cartesian coordinates. If None, then the molecules' coordinates will be used
    """
    cdef int offset = 6 if extended else 0
    cdef np.ndarray[np.float_t] bat = np.zeros((self.natoms*3-6+offset,))
    cdef int a1, a2, a3, a4, n, batInd

    bat[offset] = distance(XYZ[self.rootInd[0]],XYZ[self.rootInd[1]])
    bat[offset+1] = distance(XYZ[self.rootInd[1]],XYZ[self.rootInd[2]])
    bat[offset+2] = angle(XYZ[self.rootInd[0]],\
      XYZ[self.rootInd[1]],XYZ[self.rootInd[2]])

    for n in range(len(self._torsionIndL)):
      (a1,a2,a3,a4) = self._torsionIndL[n]
      batInd = offset+3*n+3
      bat[batInd:batInd+3] = BAT4(XYZ[a1],XYZ[a2],XYZ[a3],XYZ[a4])
      if self._firstTorsionTInd[n] != n:
        bat[batInd+2] -= bat[offset+5+3*self._firstTorsionTInd[n]]

    if extended:
      bat[:offset] = extended_coordinates(XYZ[self.rootInd[0]],\
        XYZ[self.rootInd[1]],XYZ[self.rootInd[2]])

    return bat

  @cython.boundscheck(False)
  @cython.wraparound(False)
  @cython.cdivision(True)
  def Cartesian(self, np.ndarray[np.float_t] BAT):
    """
    Conversion from (internal or extended) Bond-Angle-Torsion 
    to Cartesian coordinates
    """
    cdef int offset = 6 if len(BAT)==(3*self.natoms) else 0

    cdef np.ndarray[np.float_t] p1, p2, p3
    p1 = np.array([0.,0.,0.])
    p2 = np.array([0.,0.,BAT[offset]])
    p3 = np.array([BAT[offset+1]*sin(BAT[offset+2]), 0., \
      BAT[offset]-BAT[offset+1]*cos(BAT[offset+2])])

    # If appropriate, rotate and translate the first three atoms
    if offset==6:
      (p1,p2,p3) = rotate123(p1,p2,p3,np.array(BAT[:3]),BAT[3],BAT[4],BAT[5])

    cdef np.ndarray[np.float_t, ndim=2] XYZ = np.zeros((self.natoms,3))
    
    XYZ[self.rootInd[0]] = p1
    XYZ[self.rootInd[1]] = p2
    XYZ[self.rootInd[2]] = p3

    cdef int a1, a2, a3, a4, batInd
    cdef float torsion

    for n in range(len(self._torsionIndL)):
      (a1,a2,a3,a4) = self._torsionIndL[n]
      batInd = offset+3*n+3
      if self._firstTorsionTInd[n] == n:
        torsion = BAT[batInd+2]
      else:
        torsion = BAT[batInd+2] + BAT[offset+5+3*self._firstTorsionTInd[n]]
      XYZ[a1] = get_XYZ1(XYZ[a2],XYZ[a3],XYZ[a4],BAT[batInd],BAT[batInd+1],torsion)

    return XYZ

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
    self = converter(universe, molecule)

    # This tests a conversion to BAT coordinates and back
    BAT = self.BAT(original_xyz, extended=True)
    new_xyz = self.Cartesian(BAT)
    print sum(sum(new_xyz - original_xyz))