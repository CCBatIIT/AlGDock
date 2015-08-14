#!/usr/bin/env python

import numpy as np
from Scientific.Geometry.Objects3D import Sphere, Cone, Plane, Line, \
                                          rotatePoint
from Scientific.Geometry import Vector
import MMTK

class converter():
  """
  Interconverts Cartesian and Bond-Angle-Torsion coordinates
  """
  def __init__(self, universe, molecule, initial_atom=None):
    self.universe = universe
    self.molecule = molecule
    self._converter_setup(initial_atom)

  def _converter_setup(self, initial_atom=None):
    atom_mass = lambda atom:atom.mass()
    terminal_atoms = sorted(\
      [a for a in self.molecule.atoms if len(a.bondedTo())==1], \
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
    self.root = root

    # Construct a list of torsion angles
    self._torsionL = []
    selected = [a for a in root]
    while len(selected)<self.molecule.numberOfAtoms():
      (a1,a2,a3,a4) = self._find_dihedral(selected)
      self._torsionL.append((a1,a2,a3,a4))
      selected.append(a1)
    # If _firstTorsionInd is not equal to the list index,
    # then the dihedrals will likely be correlated and it is more appropriate
    # to use a relative phase angle
    prior_atoms = [(a2,a3,a4) for (a1,a2,a3,a4) in self._torsionL]
    self._firstTorsionInd = [prior_atoms.index(prior_atoms[n]) \
      for n in range(len(prior_atoms))]
    self.ntorsions = self.molecule.numberOfAtoms()-3

  def _find_dihedral(self, selected):
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

  def BAT(self, extended=False, Cartesian=None):
    """
    Conversion from Cartesian to Bond-Angle-Torsion coordinates
    :param extended: whether to include external coordinates or not
    :param Cartesian: Cartesian coordinates. If None, then the molecules' coordinates will be used
    """
    if Cartesian is not None:
      self.universe.setConfiguration(MMTK.Configuration(self.universe,Cartesian))
    root = [self.universe.distance(self.root[0],self.root[1]),\
      self.universe.distance(self.root[1],self.root[2]),\
      self.universe.angle(self.root[0],self.root[1],self.root[2])]
    bonds = [self.universe.distance(a1,a2) for (a1,a2,a3,a4) in self._torsionL]
    angles = [self.universe.angle(a1,a2,a3) for (a1,a2,a3,a4) in self._torsionL]
    torsions = [self.universe.dihedral(a1,a2,a3,a4) \
      for (a1,a2,a3,a4) in self._torsionL]
    phase_torsions = [(torsions[n] - torsions[self._firstTorsionInd[n]]) \
      if self._firstTorsionInd[n]!=n else torsions[n] \
      for n in range(len(torsions))]
    internal = root + bonds + angles + phase_torsions
  
    if not extended:
      return internal

    # The rotation axis is a normalized vector pointing from atom 0 to 1
    # It is described in two degrees of freedom by the polar angle and azimuth
    e = (self.root[1].position()-self.root[0].position()).normal()
    phi = np.arctan2(e[1],e[0]) # Polar angle
    theta = np.arccos(e[2]) # Azimuthal angle
    # Rotation to the z axis
    cp = np.cos(phi)
    sp = np.sin(phi)
    ct = np.cos(theta)
    st = np.sin(theta)
    Rz = np.array([[cp*ct,ct*sp,-st],[-sp,cp,0],[cp*st,sp*st,ct]])
    pos2 = Rz.dot(np.array(self.root[2].position()-self.root[0].position()))
    # Angle about the rotation axis
    omega = np.arctan2(pos2[1],pos2[0])
    external = list(self.root[0].position()) + [phi, theta, omega]
    return external + internal

  def Cartesian(self, BAT):
    """
    Conversion from (internal or extended) Bond-Angle-Torsion 
    to Cartesian coordinates
    """
    # Arrange BAT coordinates in convenient arrays
    offset = 6 if len(BAT)==(3*self.molecule.numberOfAtoms()) else 0
    bonds = BAT[offset+3:self.ntorsions+offset+3]
    angles = BAT[self.ntorsions+offset+3:2*self.ntorsions+offset+3]
    phase_torsions = BAT[2*self.ntorsions+offset+3:]
    torsions = [(phase_torsions[n] + phase_torsions[self._firstTorsionInd[n]]) \
      if self._firstTorsionInd[n]!=n else phase_torsions[n] \
      for n in range(self.ntorsions)]
    
    # Determine the positions of the first three atoms
    p1 = Vector(0,0,0) # First atom at origin
    p2 = Vector(0,0,BAT[offset]) # Second atom along z-axis
    # Third atom in xz-plane
    sphere = Sphere(p2, BAT[offset+1])
    cone = Cone(p2, -p2, BAT[offset+2])
    plane = Plane(p1, Vector(0,1,0))
    p3 = sphere.intersectWith(cone).intersectWith(plane)[0]

    # If appropriate, rotate and translate the first three atoms
    if offset==6:
      p1 = np.array(p1)
      p2 = np.array(p2)
      p3 = np.array(p3)
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
    
    self.root[0].setPosition(Vector(p1))
    self.root[1].setPosition(Vector(p2))
    self.root[2].setPosition(Vector(p3))

    # Add fourth and remaining atoms
    for ((a1,a2,a3,a4), bond, angle, torsion) in zip(self._torsionL,bonds,angles,torsions):
      sphere = Sphere(a2.position(), bond)
      cone = Cone(a2.position(), a3.position()-a2.position(), angle)
      plane123 = Plane(a4.position(), a3.position(), a2.position())
      points = sphere.intersectWith(cone).intersectWith(plane123)
      for p in points:
        if Plane(a3.position(), a2.position(), p).normal * plane123.normal > 0:
          break
      # The line points in the opposite direction to the ZMatrix constructor from
      # MMTK, but it seems to be correct
      p = rotatePoint(p, Line(a2.position(), a2.position()-a3.position()), torsion)
      a1.setPosition(p)
    
    return self.universe.configuration().array

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
    BAT = self.BAT(extended=True)
    self.Cartesian(BAT)
    print sum(sum(self.universe.configuration().array - original_xyz))

    # This rotates the last primary torsion
    torsion_ind = self._firstTorsionInd[-1]
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
    self.showMolecule(dcdFN='rotation.dcd')
    os.remove('rotation.dcd')
