import os, sys
import numpy as np

try:
    import MMTK
    from MMTK.ForceFields import Amber12SBForceField
except ImportError:
    print('MMTK not found.')
try:
    from openmm.app import AmberPrmtopFile, AmberInpcrdFile, Simulation, NoCutoff
    from openmm import *
    from openmm.unit import *
except ImportError:
    print('OpenMM not found.')
    openmm = None

term_map = {
  'cosine dihedral angle': 'MM',
  'electrostatic/pair sum': 'MM',
  'harmonic bond': 'MM',
  'harmonic bond angle': 'MM',
  'Lennard-Jones': 'MM',
  'OpenMM': 'MM',
  'OBC': 'OBC',
  'OBC_desolv': 'OBC',
  'site': 'site',
  'sLJr': 'sLJr',
  'sELE': 'sELE',
  'sLJa': 'sLJa',
  'LJr': 'LJr',
  'LJa': 'LJa',
  'ELE': 'ELE',
  'pose dihedral angle': 'k_angular_int',
  'pose external dihedral': 'k_angular_ext',
  'pose external distance': 'k_spatial_ext',
  'pose external angle': 'k_angular_ext'
}

ligand_database = 'ligand.db'
gaff_forcefield_file = 'gaff2.dat'
ligand_frcmod_file = 'ligand.frcmod'

MMTK.Database.molecule_types.directory = os.path.dirname(os.path.abspath(ligand_database))
molecule = MMTK.Molecule(os.path.basename(ligand_database))
universe = MMTK.Universe.InfiniteUniverse()
universe.addObject(molecule)
forcefields_gaff = Amber12SBForceField(parameter_file = gaff_forcefield_file, mod_files = [ligand_frcmod_file])
universe.setForceField(forcefields_gaff)

eT = universe.energyTerms()
E = {}

from AlGDock.BindingPMF import scalables
for term in (scalables):
  E[term] = np.zeros(1, dtype=float)

for (key, value) in eT.iteritems():
  if key == 'electrostatic':
    pass  # For some reason, MMTK double-counts electrostatic energies
  elif key.startswith('pose'):
    # For pose restraints, the energy is per spring constant unit
    E[term_map[key]][c] += value / params_full[term_map[key]]
  else:
    try:
      E[term_map[key]][c] += value
    except KeyError:
      print(key)
      print('Keys in eT', eT.keys())
      print('Keys in term map', term_map.keys())
      print('Keys in E', E.keys())
      raise Exception('key not found in term map or E')
      
print(E)
