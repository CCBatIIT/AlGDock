#!/usr/bin/env python

# TODO: Find a sincos function that simultaneously calculates sine and cosine?

import cython

import numpy as np
cimport numpy as np

from cpython cimport bool

from libc.math cimport sin
from libc.math cimport cos
from libc.math cimport sqrt
from libc.math cimport acos
from libc.math cimport atan2

import MMTK

### Vector functions
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline dotp(double[:] v1, double[:] v2):
  return (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2])

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline norm2(double[:] v):
  return (v[0]*v[0] + v[1]*v[1] + v[2]*v[2])

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline normalize(np.ndarray[np.double_t] v1):
  return v1/sqrt(norm2(v1))

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline cross(double[:] v1, double[:] v2):
  return np.array([v1[1]*v2[2]-v1[2]*v2[1], \
    v1[2]*v2[0]-v1[0]*v2[2], \
    v1[0]*v2[1]-v1[1]*v2[0]])

### BAT coordinate measurement
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline distance(double[:] p1, double[:] p2):
  cdef double d2, delta_c
  cdef int c
  for c in range(3):
    delta_c = p2[c] - p1[c]
    d2 += delta_c*delta_c
  return sqrt(d2)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline angle(np.ndarray[np.double_t] p1, np.ndarray[np.double_t] p2, \
    np.ndarray[np.double_t] p3):
  cdef np.ndarray[np.double_t] v1 = p2 - p1
  cdef np.ndarray[np.double_t] v2 = p2 - p3
  return acos(max(-1.,min(1.,dotp(v1,v2)/\
    sqrt(norm2(v1)*norm2(v2)))))

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef dihedral(np.ndarray[np.double_t] p1, np.ndarray[np.double_t] p2, \
    np.ndarray[np.double_t] p3, np.ndarray[np.double_t] p4):
  cdef np.ndarray[np.double_t] v1, v2, v3, a, b
  cdef double c, s
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
cdef extended_coordinates(np.ndarray[np.double_t] p1, \
    np.ndarray[np.double_t] p2, np.ndarray[np.double_t] p3):
  # The rotation axis is a normalized vector pointing from atom 0 to 1
  # It is described in two degrees of freedom by the polar angle and azimuth
  cdef np.ndarray[np.double_t] e = normalize(p2-p1)
  cdef double phi = atan2(e[1],e[0]) # Polar angle
  cdef double theta = acos(e[2]) # Azimuthal angle
  # Rotation to the z axis
  cdef double cp = cos(phi)
  cdef double sp = sin(phi)
  cdef double ct = cos(theta)
  cdef double st = sin(theta)
  cdef np.ndarray[np.double_t, ndim=2] Rz = np.array([[cp*ct,ct*sp,-st],[-sp,cp,0],[cp*st,sp*st,ct]])
  cdef np.ndarray[np.double_t] pos2 = Rz.dot(np.array(p3-p1))
  # Angle about the rotation axis
  cdef double omega = atan2(pos2[1],pos2[0])
  cdef np.ndarray[np.double_t] coords = np.zeros((6,))
  coords[:3] = p1
  coords[3] = phi
  coords[4] = theta
  coords[5] = omega
  return coords

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef rotate123(np.ndarray[np.double_t] p1, \
    np.ndarray[np.double_t] p2, np.ndarray[np.double_t] p3, \
    np.ndarray[np.double_t] origin, \
    double phi, double theta, double omega):
  # Rotate the third atom by the appropriate value
  cdef double co = cos(omega)
  cdef double so = sin(omega)
  cdef np.ndarray[np.double_t, ndim=2] Romega = np.array([[co, -so, 0],[so, co, 0], [0, 0, 1]])
  p3 = Romega.dot(p3)
  # Rotate the second two atoms to point in the right direction
  cdef double cp = cos(phi)
  cdef double sp = sin(phi)
  cdef double ct = cos(theta)
  cdef double st = sin(theta)
  cdef np.ndarray[np.double_t, ndim=2] Re = np.array([[cp*ct,-sp,cp*st],[ct*sp,cp,sp*st],[-st,0,ct]])
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

  # These attributes are already declared in the pxd file
  # cdef readonly int natoms, ntorsions
  # cdef np.ndarray rootInd, _torsionIndL, _firstTorsionTInd

  def __init__(self, universe, molecule, initial_atom=None):
    self.molecule = molecule
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

    self.rootInd = np.array([r.index for r in root], dtype=int)
    self._torsionIndL = np.array(\
      [[a.index for a in tset] for tset in torsionL], dtype=int)
    self._firstTorsionTInd = np.array([prior_atoms.index(prior_atoms[n]) \
      for n in range(len(prior_atoms))], dtype=int)
    self.ntorsions = self.natoms-3

  def getFirstTorsionInds(self, extended):
    offset = 6 if extended else 0
    torsionInds = np.array(range(offset+5,self.natoms*3,3))
    primaryTorsions = sorted(list(set(self._firstTorsionTInd)))
    return list(torsionInds[primaryTorsions])

  @cython.boundscheck(False)
  @cython.wraparound(False)
  @cython.cdivision(True)
  cpdef BAT(self, double[:,:] XYZ, bool extended):
    """
    Conversion from Cartesian to Bond-Angle-Torsion coordinates
    :param extended: whether to include external coordinates or not
    :param Cartesian: Cartesian coordinates. If None, then the molecules' coordinates will be used
    """
    cdef int offset = 6 if extended else 0
    cdef np.ndarray[np.double_t] bat
    cdef np.ndarray[np.double_t] p1, p2, p3, p4
    cdef int a1, a2, a3, a4, n, i, batInd

    cdef np.ndarray[np.double_t] v1, v2, v3, a, b
    cdef double c, s, norm_v1_2, norm_v2_2

    bat = np.zeros((self.natoms*3-6+offset,))
    p1 = np.zeros((3,))
    p2 = np.zeros((3,))
    p3 = np.zeros((3,))
    p4 = np.zeros((3,))

    # Memory Views
    cpdef double[:] bat_v = bat
    cdef double[:] p1_v = p1
    cdef double[:] p2_v = p2
    cdef double[:] p3_v = p3
    cdef double[:] p4_v = p4

    for i in range(3):
      p1_v[i] = XYZ[self.rootInd[0],i]
      p2_v[i] = XYZ[self.rootInd[1],i]
      p3_v[i] = XYZ[self.rootInd[2],i]

    if extended:
      bat[:offset] = extended_coordinates(p1,p2,p3)
    bat_v[offset] = distance(p1,p2)
    bat_v[offset+1] = distance(p2,p3)
    bat_v[offset+2] = angle(p1,p2,p3)

    for n in range(self._torsionIndL.shape[0]):
      for i in range(3):
        p1_v[i] = XYZ[self._torsionIndL[n,0],i]
        p2_v[i] = XYZ[self._torsionIndL[n,1],i]
        p3_v[i] = XYZ[self._torsionIndL[n,2],i]
        p4_v[i] = XYZ[self._torsionIndL[n,3],i]

      v1 = p2 - p1
      v2 = p2 - p3
      v3 = p3 - p4
      a = normalize(cross(v1,v2))
      b = normalize(cross(v3,v2))
      c = dotp(a,b)
      norm_v1_2 = norm2(v1)
      norm_v2_2 = norm2(v2)
      s = dotp(cross(b,a),v2/sqrt(norm_v2_2))

      batInd = offset+3*n+3
      bat_v[batInd] = sqrt(norm_v1_2)
      bat_v[batInd+1] = acos(max(-1.,min(1.,dotp(v1,v2)/\
        sqrt(norm_v1_2*norm_v2_2))))
      bat_v[batInd+2] = atan2(s,c)
      if self._firstTorsionTInd[n] != n:
        bat_v[batInd+2] -= bat_v[offset+5+3*self._firstTorsionTInd[n]]

    return bat

  @cython.boundscheck(False)
  @cython.wraparound(False)
  @cython.cdivision(True)
  cpdef Cartesian(self, double[:] BAT):
    """
    Conversion from (internal or extended) Bond-Angle-Torsion 
    to Cartesian coordinates
    """
    cdef int offset, batInd, firstTorsionInd, i
    cdef np.ndarray[np.double_t] p1, p2, p3, p4, n23, v21
    cdef np.ndarray[np.double_t, ndim=2] XYZ
    cdef double s, c, bond, angle, torsion

    offset = 6 if BAT.shape[0]==(3*self.natoms) else 0

    p1 = np.array([0.,0.,0.])
    p2 = np.array([0.,0.,BAT[offset]])
    p3 = np.array([BAT[offset+1]*sin(BAT[offset+2]), 0., \
      BAT[offset]-BAT[offset+1]*cos(BAT[offset+2])])
    p4 = np.array([0.,0.,0.])

    # If appropriate, rotate and translate the first three atoms
    if offset==6:
      (p1,p2,p3) = rotate123(p1,p2,p3,np.array(BAT[:3]),BAT[3],BAT[4],BAT[5])

    XYZ = np.zeros((self.natoms,3))

    # Memoryviews for faster indexing
    cdef double[:,:] XYZ_v = XYZ
    cdef double[:] p1_v = p1
    cdef double[:] p2_v = p2
    cdef double[:] p3_v = p3
    cdef double[:] p4_v = p4

    for i in range(3):
      XYZ_v[self.rootInd[0],i] = p1[i]
      XYZ_v[self.rootInd[1],i] = p2[i]
      XYZ_v[self.rootInd[2],i] = p3[i]

    for n in range(self._torsionIndL.shape[0]):
      for i in range(3):
        p2_v[i] = XYZ_v[self._torsionIndL[n,1],i]
        p3_v[i] = XYZ_v[self._torsionIndL[n,2],i]
        p4_v[i] = XYZ_v[self._torsionIndL[n,3],i]

      batInd = offset+3*n+3
      bond = BAT[batInd]
      angle = BAT[batInd+1]
      if self._firstTorsionTInd[n]==n:
        torsion = BAT[batInd+2]
      else:
        firstTorsionInd = offset+5+3*self._firstTorsionTInd[n]
        torsion = BAT[batInd+2] + BAT[firstTorsionInd]

      n23 = normalize(p3-p2)

      s = sin(angle)
      c = cos(angle)
      v21 = (bond*c)*n23 - (bond*s)*cross(normalize(cross(p4-p3,n23)),n23)

      s = sin(torsion)
      c = cos(torsion)
      p1 = p2 - cross(n23,v21)*s + dotp(n23,v21)*n23*(1.0-c) + v21*c
      XYZ[self._torsionIndL[n,0],:] = p1

    return XYZ

  def showMolecule(self, colorBy=None, label=False, dcdFN=None):
    """
    Opens the molecule in VMD
    :param colorBy: color atoms by 'Occupancy', or 'Beta'. None uses default colors.
    """
    # Write PDB file
    # To set Occupancy, change atom.occupancy
    # To set Beta, change atom.temperature_factor
    import os.path
    pdbFN = os.path.join(MMTK.Database.molecule_types.directory,'showMolecule.pdb')
    outF = MMTK.PDB.PDBOutputFile(pdbFN)
    outF.write(self.molecule)
    outF.close()
    # Write VMD script
    script  = 'set ligand [mol new '+pdbFN+']\n'
    if colorBy is not None:
      script += 'mol modcolor 0 $ligand '+colorBy+'\n'
    script += 'mol modstyle 0 0 CPK 1.000000 0.300000 10.000000 10.000000\n'
    if label:
      script += """
proc label_atoms { molid seltext } {
  set sel [atomselect $molid $seltext] 
  set atomlist [$sel list] 
  foreach {atom} $atomlist { 
    set atomlabel [format "%d/%d" $molid $atom]
    label add Atoms $atomlabel 
  } 
  $sel delete 
} 
label_atoms 0 all
"""
    if dcdFN is not None:
      script += 'animate delete all $ligand\n'
      script += 'mol addfile '+dcdFN+' type dcd waitfor all\n'
    scriptF = open('showMolecule.vmd','w')
    scriptF.write(script)
    scriptF.close()
    # Find and run vmd
    import AlGDock
    vmdCommand = AlGDock.findPath(AlGDock.search_paths['vmd'])
    import subprocess
    subprocess.call([vmdCommand, '-e', 'showMolecule.vmd'])
    # Remove files
    os.remove(pdbFN)
    os.remove('showMolecule.vmd')
