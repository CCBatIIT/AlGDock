#!/usr/bin/env python

import numpy as np
from Scientific.Geometry.Objects3D import Sphere, Cone, Plane, Line, \
                                          rotatePoint
from Scientific.Geometry import Vector
import MMTK

# Vector functions
def normalize(v1):
  return v1/np.sqrt(np.sum(v1*v1))

def cross(v1,v2):
  return np.array([v1[1]*v2[2]-v1[2]*v2[1], \
    v1[2]*v2[0]-v1[0]*v2[2], \
    v1[0]*v2[1]-v1[1]*v2[0]])

def distance(p1,p2):
  v1 = p2 - p1
  return np.sqrt(np.sum(v1*v1))

def angle(p1,p2,p3):
  v1 = p2 - p1
  v2 = p2 - p3
  return np.arccos(max(-1.,min(1.,np.sum(v1*v2)/\
    np.sqrt(np.sum(v1*v1)*np.sum(v2*v2)))))

def dihedral(p1,p2,p3,p4):
  v1 = p2 - p1
  v2 = p2 - p3
  v3 = p3 - p4
  a = normalize(cross(v1,v2))
  b = normalize(cross(v3,v2))
  c = np.sum(a*b)
  s = np.sum(cross(b,a)*normalize(v2))
  return np.arctan2(s,c)

def BAT4(p1,p2,p3,p4):
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

# Main converter class
class converter():
  """
  Interconverts Cartesian and Bond-Angle-Torsion coordinates
  """
  def __init__(self, universe, molecule, initial_atom=None):
    self.universe = universe
    self.molecule = molecule
    self.natoms = universe.numberOfAtoms()
    self._converter_setup(initial_atom)

  def _converter_setup(self, initial_atom=None):
    atom_name = lambda atom:atom.fullName()
    atom_mass = lambda atom:atom.mass()
    terminal_atoms = sorted(\
      [a for a in self.molecule.atoms if len(a.bondedTo())==1], key=atom_name)
    terminal_atoms = sorted(terminal_atoms, key=atom_mass)
    if (initial_atom is None):
      # Select the heaviest root atoms from the heaviest terminal atom
      root = [terminal_atoms[-1]]
    else:
      if initial_atom in terminal_atoms:
        root = [initial_atom]
      else:
        raise Exception('Initial atom is not a terminal atom')
    self.initial_atom = initial_atom
    
    attached_to_zero = sorted(root[0].bondedTo(), key=atom_name)
    attached_to_zero = sorted(attached_to_zero, key=atom_mass)
    root.append(attached_to_zero[-1])
    attached_to_one = sorted([a for a in root[-1].bondedTo() \
      if (a not in root) and (a not in terminal_atoms)], key=atom_name)
    attached_to_one = sorted(attached_to_one, key=atom_mass)
    root.append(attached_to_one[-1])

    def _find_dihedral(selected):
      """
      Finds a dihedral angle adjacent to the selected atoms that includes a new atom
      :param selected: a list of atoms that have already been selected
      :returns: a list of atoms that includes the new atom and its neighboring selected atoms
      """
      atom_name = lambda atom:atom.fullName()
      atom_mass = lambda atom:atom.mass()
      # Loop over possible nearest neighbors
      for a2 in selected:
        # Find the new atom
        attached_to_a2 = sorted([a for a in a2.bondedTo() \
            if a not in selected], key=atom_name)
        for a1 in sorted(attached_to_a2, key=atom_mass, reverse=True):
          # Find the third atom
          attached_to_a3 = sorted([a for a in a2.bondedTo() \
              if (a in selected) and (a!=a1)], key=atom_name)
          for a3 in sorted(attached_to_a3, key=atom_mass, reverse=True):
            # Find the last atom
            attached_to_a4 = sorted([a for a in a3.bondedTo() \
                if (a in selected) and (a!=a2)], key=atom_name)
            for a4 in sorted(attached_to_a4, key=atom_mass, reverse=True):
              return (a1,a2,a3,a4)
      print 'Selected atoms:', selected
      raise Exception('No new dihedral angle found!')

    # Construct a list of torsion angles
    torsionL = []
    selected = [a for a in root]
    while len(selected)<self.universe.numberOfAtoms():
      (a1,a2,a3,a4) = _find_dihedral(selected)
      torsionL.append((a1,a2,a3,a4))
      selected.append(a1)
    # If _firstTorsionTInd is not equal to the list index,
    # then the dihedrals will likely be correlated and it is more appropriate
    # to use a relative phase angle
    prior_atoms = [sorted([a2.index,a3.index]) for (a1,a2,a3,a4) in torsionL]

    self.rootInd = [r.index for r in root]
    self._torsionIndL = [[a.index for a in tset] for tset in torsionL]
    self._firstTorsionTInd = [prior_atoms.index(prior_atoms[n]) \
      for n in range(len(prior_atoms))]
    self.ntorsions = self.natoms-3

  def getFirstTorsionInds(self, extended):
    """
    Indices of the first torsions in the BAT array
    """
    offset = 6 if extended else 0
    torsionInds = np.array(range(offset+5,self.natoms*3,3))
    primaryTorsions = sorted(list(set(self._firstTorsionTInd)))
    return list(torsionInds[primaryTorsions])

  def BAT(self, XYZ, extended=False):
    """
    Conversion from Cartesian to Bond-Angle-Torsion coordinates
    :param extended: whether to include external coordinates or not
    :param Cartesian: Cartesian coordinates. If None, then the molecules' coordinates will be used
    """
    root = [distance(XYZ[self.rootInd[0]],XYZ[self.rootInd[1]]),\
      distance(XYZ[self.rootInd[1]],XYZ[self.rootInd[2]]),\
      angle(XYZ[self.rootInd[0]],XYZ[self.rootInd[1]],XYZ[self.rootInd[2]])]

    import itertools
    internal = root + \
      [val for val in itertools.chain.from_iterable([\
        BAT4(XYZ[a1],XYZ[a2],XYZ[a3],XYZ[a4]) \
          for (a1,a2,a3,a4) in self._torsionIndL])]

    torsions = internal[5::3]
    phase_torsions = [(torsions[n] - torsions[self._firstTorsionTInd[n]]) \
      if self._firstTorsionTInd[n]!=n else torsions[n] \
      for n in range(len(torsions))]
    internal[5::3] = phase_torsions

    if not extended:
      return np.array(internal)

    external = self.extended_coordinates(XYZ[self.rootInd[0]], \
      XYZ[self.rootInd[1]], XYZ[self.rootInd[2]])
    return np.array(list(external)+list(internal))
    
  def extended_coordinates(self,p1,p2,p3):
    # The rotation axis is a normalized vector pointing from atom 0 to 1
    # It is described in two degrees of freedom by the polar angle and azimuth
    e = normalize(p2-p1)
    phi = np.arctan2(e[1],e[0]) # Polar angle
    theta = np.arccos(e[2]) # Azimuthal angle
    # Rotation to the z axis
    cp = np.cos(phi)
    sp = np.sin(phi)
    ct = np.cos(theta)
    st = np.sin(theta)
    Rz = np.array([[cp*ct,ct*sp,-st],[-sp,cp,0],[cp*st,sp*st,ct]])
    pos2 = Rz.dot(np.array(p3-p2))
    # Angle about the rotation axis
    omega = np.arctan2(pos2[1],pos2[0])
    return np.array(list(p1) + [phi, theta, omega])
    
  def Cartesian(self, BAT):
    """
    Conversion from (internal or extended) Bond-Angle-Torsion 
    to Cartesian coordinates
    """
    # Arrange BAT coordinates in convenient arrays
    offset = 6 if len(BAT)==(3*self.natoms) else 0
    bonds = BAT[offset+3::3]
    angles = BAT[offset+4::3]
    phase_torsions = BAT[offset+5::3]
    torsions = [(phase_torsions[n] + phase_torsions[self._firstTorsionTInd[n]]) \
      if self._firstTorsionTInd[n]!=n else phase_torsions[n] \
      for n in range(self.ntorsions)]
    
    p1 = np.array([0.,0.,0.])
    p2 = np.array([0.,0.,BAT[offset]])
    p3 = np.array([BAT[offset+1]*np.sin(BAT[offset+2]), 0., \
      BAT[offset]-BAT[offset+1]*np.cos(BAT[offset+2])])

    # If appropriate, rotate and translate the first three atoms
    if offset==6:
      # Rotate the third atom by the appropriate value
      (phi,theta,omega) = BAT[3:6]
      co = np.cos(omega)
      so = np.sin(omega)
      Romega = np.array([[co, -so, 0],[so, co, 0], [0, 0, 1]])
      p3 = Romega.dot(p3)
      # Rotate the second two atoms to point in the right direction
      cp = np.cos(phi)
      sp = np.sin(phi)
      ct = np.cos(theta)
      st = np.sin(theta)
      Re = np.array([[cp*ct,-sp,cp*st],[ct*sp,cp,sp*st],[-st,0,ct]])
      p2 = Re.dot(p2)
      p3 = Re.dot(p3)
      # Translate the first three atoms by the origin
      origin = np.array(BAT[:3])
      p1 += origin
      p2 += origin
      p3 += origin

    XYZ = np.zeros((self.natoms,3))
    
    XYZ[self.rootInd[0]] = p1
    XYZ[self.rootInd[1]] = p2
    XYZ[self.rootInd[2]] = p3

    for ((a1,a2,a3,a4), bond, angle, torsion) in \
        zip(self._torsionIndL,bonds,angles,torsions):
      sphere = Sphere(Vector(XYZ[a2]), bond)
      cone = Cone(Vector(XYZ[a2]), Vector(XYZ[a3]-XYZ[a2]), angle)
      plane123 = Plane(Vector(XYZ[a4]), Vector(XYZ[a3]), Vector(XYZ[a2]))
      points = sphere.intersectWith(cone).intersectWith(plane123)
      p = points[0] if (Plane(Vector(XYZ[a3]), Vector(XYZ[a2]), points[0]).normal*plane123.normal)>0 else points[1]
      p = rotatePoint(Vector(p), Line(Vector(XYZ[a2]), Vector(XYZ[a2]-XYZ[a3])), torsion)
      XYZ[a1] = p.array

    return XYZ

    for ((a1,a2,a3,a4), bond, angle, torsion) in \
        zip(self._torsionIndL,bonds,angles,torsions):

      p2 = XYZ[a2]
      p3 = XYZ[a3]
      p4 = XYZ[a4]

      # circle = sphere.intersectWith(cone)
      n23 = normalize(p3-p2)

      # points = circle.intersectWith(plane123)
      # plane.intersectWith(Plane(circle.center, circle.normal)) is a line
      # line_direction = cross(normalize(cross(p4-p3,n23)),n23)

      # Rotate the point about the p2-p3 axis by the torsion angle
      v21 = (bond*np.cos(angle))*n23 - (bond*np.sin(angle))*cross(normalize(cross(p4-p3,n23)),n23)
      s = np.sin(torsion)
      c = np.cos(torsion)
      XYZ[a1] = p2 - cross(n23,v21)*s + np.sum(n23*v21)*n23*(1.0-c) + v21*c

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

    # This rotates a random primary torsion
    from random import randrange
    firstTorsionInds = self.getFirstTorsionInds(True)
    BAT_ind = firstTorsionInds[randrange(len(firstTorsionInds))]
    confs = []
    for torsion_offset in np.linspace(0,2*np.pi):
      BAT_n = np.array([BAT[ind] if ind!=BAT_ind else BAT[ind] + torsion_offset \
        for ind in range(len(BAT))])
      XYZ = self.Cartesian(BAT_n)
      confs.append(XYZ)

    import AlGDock.IO
    IO_dcd = AlGDock.IO.dcd(molecule)
    IO_dcd.write('rotation.dcd', confs)
    self.showMolecule(dcdFN='rotation.dcd')
    os.remove('rotation.dcd')

#  [[51, 10, 5, 46],
#   [2, 5, 10, 51],
#   [4, 5, 10, 51],
