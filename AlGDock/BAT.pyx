#!/usr/bin/env python

import numpy as np
cimport numpy as np
import cython

ctypedef np.float_t float_t
ctypedef np.int_t int_t

import MMTK

# Vector functions
@cython.boundscheck(False)
cdef inline normalize(np.ndarray[np.float_t] v1):
  return v1/np.sqrt(np.sum(v1*v1))

@cython.boundscheck(False)
cdef inline cross(np.ndarray[np.float_t] v1, np.ndarray[np.float_t] v2):
  return np.array([v1[1]*v2[2]-v1[2]*v2[1], \
    v1[2]*v2[0]-v1[0]*v2[2], \
    v1[0]*v2[1]-v1[1]*v2[0]])

@cython.boundscheck(False)
cdef inline distance(np.ndarray[np.float_t] p1,np.ndarray[np.float_t] p2):
  cdef np.ndarray[np.float_t] v1 = p2 - p1
  return np.sqrt(np.sum(v1*v1))

@cython.boundscheck(False)
cdef inline angle(np.ndarray[np.float_t] p1, np.ndarray[np.float_t] p2, \
    np.ndarray[np.float_t] p3):
  cdef np.ndarray[np.float_t] v1 = p2 - p1
  cdef np.ndarray[np.float_t] v2 = p2 - p3
  return np.arccos(max(-1.,min(1.,np.sum(v1*v2)/\
    np.sqrt(np.sum(v1*v1)*np.sum(v2*v2)))))

@cython.boundscheck(False)
cdef dihedral(np.ndarray[np.float_t] p1, np.ndarray[np.float_t] p2, \
    np.ndarray[np.float_t] p3, np.ndarray[np.float_t] p4):
  cdef np.ndarray[np.float_t] v1, v2, v3, a, b
  cdef float c, s
  v1 = p2 - p1
  v2 = p2 - p3
  v3 = p3 - p4
  a = normalize(cross(v1,v2))
  b = normalize(cross(v3,v2))
  c = np.sum(a*b)
  s = np.sum(cross(b,a)*normalize(v2))
  return np.arctan2(s,c)

@cython.boundscheck(False)
cdef BAT4(np.ndarray[np.float_t] p1, np.ndarray[np.float_t] p2, \
    np.ndarray[np.float_t] p3, np.ndarray[np.float_t] p4):
  cdef np.ndarray[np.float_t] v1, v2, v3, a, b
  cdef float c, s, norm_v1_2, norm_v2_2
  v1 = p2 - p1
  v2 = p2 - p3
  v3 = p3 - p4
  a = normalize(cross(v1,v2))
  b = normalize(cross(v3,v2))
  c = np.sum(a*b)
  norm_v1_2 = np.sum(v1*v1)
  norm_v2_2 = np.sum(v2*v2)
  s = np.sum(cross(b,a)*v2/np.sqrt(norm_v2_2))
  return (np.sqrt(norm_v1_2), np.arccos(max(-1.,min(1.,np.sum(v1*v2)/\
    np.sqrt(norm_v1_2*norm_v2_2)))), np.arctan2(s,c))

@cython.boundscheck(False)
cdef extended_coordinates(np.ndarray[np.float_t] p1, \
    np.ndarray[np.float_t] p2, np.ndarray[np.float_t] p3):
  # The rotation axis is a normalized vector pointing from atom 0 to 1
  # It is described in two degrees of freedom by the polar angle and azimuth
  cdef np.ndarray[np.float_t] e = normalize(p2-p1)
  cdef float phi = np.arctan2(e[1],e[0]) # Polar angle
  cdef float theta = np.arccos(e[2]) # Azimuthal angle
  # Rotation to the z axis
  cdef float cp = np.cos(phi)
  cdef float sp = np.sin(phi)
  cdef float ct = np.cos(theta)
  cdef float st = np.sin(theta)
  cdef np.ndarray[np.float_t, ndim=2] Rz = np.array([[cp*ct,ct*sp,-st],[-sp,cp,0],[cp*st,sp*st,ct]])
  cdef np.ndarray[np.float_t] pos2 = Rz.dot(np.array(p3-p1))
  # Angle about the rotation axis
  cdef float omega = np.arctan2(pos2[1],pos2[0])
  cdef np.ndarray[np.float_t] coords = np.zeros((6,))
  coords[:3] = p1
  coords[3] = phi
  coords[4] = theta
  coords[5] = omega
  return coords

@cython.boundscheck(False)
cdef get_XYZ1(np.ndarray[np.float_t] p2, \
    np.ndarray[np.float_t] p3, np.ndarray[np.float_t] p4,
    float bond, float angle, float torsion):

  # sphere.intersectWith(cone) is a circle
  cdef float from_center = bond*np.cos(angle)
  cdef float circle_radius = bond*np.sin(angle)
  cdef np.ndarray[np.float_t] circle_normal = normalize(p3-p2)
  cdef np.ndarray[np.float_t] circle_center =  p2 + circle_normal*from_center

  # plane.intersectWith(Plane(circle.center, circle.normal)) is a line
  cdef np.ndarray[np.float_t] plane123_normal = normalize(cross(p3-p4,p2-p3))
  cdef float plane123_distance_from_zero = np.sum(plane123_normal*p3)
  cdef np.ndarray[np.float_t] line_direction = cross(plane123_normal,circle_normal)
  cdef np.ndarray[np.float_t] point_in_plane123 = plane123_distance_from_zero*plane123_normal
  cdef np.ndarray[np.float_t] line_point = point_in_plane123 - \
    (np.sum(point_in_plane123*circle_normal) - \
     np.sum(circle_normal*circle_center))*circle_normal

  cdef np.ndarray[np.float_t] dv = line_point-circle_center
  dv -= np.sum(dv*line_direction)*line_direction
  cdef float x = np.sqrt(np.sum(dv*dv))
  cdef float a = np.arccos(x/circle_radius)
  cdef float along_line = np.sin(a)*circle_radius
  cdef np.ndarray[np.float_t] normal = cross(circle_normal,line_direction)
  dv = line_point-(circle_center+normal)
  dv -= np.sum(dv*line_direction)*line_direction
  cdef float d = np.sqrt(np.sum(dv*dv))
  if d>x:
    normal = -normal
  cdef np.ndarray[np.float_t] base = circle_center+x*normal
  cdef np.ndarray[np.float_t] p = base-line_direction*along_line
  if np.sum(normalize(cross(p2-p3,p-p2))*plane123_normal)<0:
    p = base+line_direction*along_line

  # Rotate the point about the p2-p3 axis by the torsion angle
  cdef np.ndarray[np.float_t] direction = normalize(p2-p3)
  cdef np.ndarray[np.float_t] vector = p-p2
  cdef float s = np.sin(torsion)
  cdef float c = np.cos(torsion)
  cdef float c1 = 1-c
  return p2 + cross(direction,vector)*s + np.sum(direction*vector)*direction*c1 + vector*c

@cython.boundscheck(False)
cdef rotate123(np.ndarray[np.float_t] p1, \
    np.ndarray[np.float_t] p2, np.ndarray[np.float_t] p3, \
    np.ndarray[np.float_t] origin, \
    float phi, float theta, float omega):
  # Rotate the third atom by the appropriate value
  cdef float co = np.cos(omega)
  cdef float so = np.sin(omega)
  cdef np.ndarray[np.float_t, ndim=2] Romega = np.array([[co, -so, 0],[so, co, 0], [0, 0, 1]])
  p3 = Romega.dot(p3)
  # Rotate the second two atoms to point in the right direction
  cdef float cp = np.cos(phi)
  cdef float sp = np.sin(phi)
  cdef float ct = np.cos(theta)
  cdef float st = np.sin(theta)
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
  # TODO: Make this attributes private
  cdef int natoms
  cdef public int ntorsions
  cdef object rootInd, _torsionIndL
  cdef public object _firstTorsionInd

  def __init__(self, universe, molecule, initial_atom=None):
    self.natoms = universe.numberOfAtoms()
    self._converter_setup(molecule.atoms, initial_atom)
    print 'Using cython converter'

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
    # If _firstTorsionInd is not equal to the list index,
    # then the dihedrals will likely be correlated and it is more appropriate
    # to use a relative phase angle
    prior_atoms = [(a2,a3,a4) for (a1,a2,a3,a4) in torsionL]

    self.rootInd = [r.index for r in root]
    self._torsionIndL = [[a.index for a in tset] for tset in torsionL]
    self._firstTorsionInd = [prior_atoms.index(prior_atoms[n]) \
      for n in range(len(prior_atoms))]
    self.ntorsions = self.natoms-3

  def BAT(self, np.ndarray[np.float_t, ndim=2] XYZ, extended=False):
    """
    Conversion from Cartesian to Bond-Angle-Torsion coordinates
    :param extended: whether to include external coordinates or not
    :param Cartesian: Cartesian coordinates. If None, then the molecules' coordinates will be used
    """
    cdef int offset = 6 if extended else 0
    cdef np.ndarray[np.float_t] bat = np.zeros((self.natoms*3-6+offset,))
    cdef int a1, a2, a3, a4, n, lastInd

    bat[offset] = distance(XYZ[self.rootInd[0]],XYZ[self.rootInd[1]])
    bat[offset+1] = distance(XYZ[self.rootInd[1]],XYZ[self.rootInd[2]])
    bat[offset+2] = angle(XYZ[self.rootInd[0]],XYZ[self.rootInd[1]],XYZ[self.rootInd[2]])

    for n in range(len(self._torsionIndL)):
      (a1,a2,a3,a4) = self._torsionIndL[n]
      lastInd = offset+3*n+6
      bat[lastInd-3:lastInd] = BAT4(XYZ[a1],XYZ[a2],XYZ[a3],XYZ[a4])
      if self._firstTorsionInd[n] != n:
        bat[lastInd-1] -= bat[offset+5+3*self._firstTorsionInd[n]]

    if extended:
      bat[:offset] = extended_coordinates(XYZ[self.rootInd[0]],XYZ[self.rootInd[1]],XYZ[self.rootInd[2]])

    return bat

  def Cartesian(self, np.ndarray[np.float_t] BAT):
    """
    Conversion from (internal or extended) Bond-Angle-Torsion 
    to Cartesian coordinates
    """
    # TODO: Don't create these lists
    # Arrange BAT coordinates in convenient lists
    cdef int offset = 6 if len(BAT)==(3*self.natoms) else 0
    cdef np.ndarray[np.float_t] bonds = BAT[offset+3::3]
    cdef np.ndarray[np.float_t] angles = BAT[offset+4::3]
    phase_torsions = BAT[offset+5::3]
    torsions = [(phase_torsions[n] + phase_torsions[self._firstTorsionInd[n]]) \
      if self._firstTorsionInd[n]!=n else phase_torsions[n] \
      for n in range(self.ntorsions)]

    cdef np.ndarray[np.float_t] p1, p2, p3
    p1 = np.array([0.,0.,0.])
    p2 = np.array([0.,0.,BAT[offset]])
    p3 = np.array([BAT[offset+1]*np.sin(BAT[offset+2]), 0., \
      BAT[offset]-BAT[offset+1]*np.cos(BAT[offset+2])])

    # If appropriate, rotate and translate the first three atoms
    if offset==6:
      (p1,p2,p3) = rotate123(p1,p2,p3,np.array(BAT[:3]),BAT[3],BAT[4],BAT[5])

    cdef np.ndarray[np.float_t, ndim=2] XYZ = np.zeros((self.natoms,3))
    
    XYZ[self.rootInd[0]] = p1
    XYZ[self.rootInd[1]] = p2
    XYZ[self.rootInd[2]] = p3

    cdef int a1, a2, a3, a4
    cdef float bond, angle, torsion

    for ((a1,a2,a3,a4), bond, angle, torsion) in \
        zip(self._torsionIndL,bonds,angles,torsions):
      XYZ[a1] = get_XYZ1(XYZ[a2],XYZ[a3],XYZ[a4],bond,angle,torsion)

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