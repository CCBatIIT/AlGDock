import numpy as np
#from typing import Optional
from simtk.openmm import System
from simtk.openmm.app import Topology


def hydrogen_mass_repartitioning(molecule, H_mass):
  H_mass_o = np.sum([atom.mass() for atom in molecule.atoms \
    if atom.type.name=='hydrogen'])
  H_mass_f = np.sum([1 for atom in molecule.atoms \
    if atom.type.name=='hydrogen'])*H_mass
  total_mass_o = np.sum([atom.mass() for atom in molecule.atoms])
  nonH_mass_o = np.sum([atom.mass() for atom in molecule.atoms \
    if atom.type.name!='hydrogen'])
  nonH_mass_f = total_mass_o - H_mass_f
  mass_ratio = nonH_mass_f / nonH_mass_o
  for atom in molecule.atoms:
    if atom.type.name == 'hydrogen':
      atom.setMass(H_mass)
    else:
      atom.setMass(mass_ratio * atom.mass())
  return molecule

def hydrogen_mass_repartitioning_openmm(molecule, omm_system, H_mass):
# def hydrogen_mass_repartitioning_openmm(molecule: Optional[Topology], omm_system: Optional[System], H_mass: int) -> Optional[System]:
# python 2 doesn't support type annotations.
  """
  Adjusts the mass of hydrogen atoms in an OpenMM System.

  Parameters:
      molecule (Topology): The topology of the system.
      omm_system (System): The OpenMM System to modify.
      H_mass (int): The new mass to set for hydrogen atoms.

  Returns:
      System: The modified OpenMM System with updated hydrogen masses.
  """
  import openmm
  H_mass_o = 0
  H_mass_f = 0
  total_mass_o = 0
  nonH_mass_o = 0
  for atom in molecule.atoms():
    atom_mass = omm_system.getParticleMass(atom.index).value_in_unit(openmm.unit.dalton)
    if atom.element == openmm.app.element.hydrogen:
      H_mass_o += atom_mass
      H_mass_f += 1 * H_mass
    else:
      nonH_mass_o += atom_mass
    total_mass_o += atom_mass

  nonH_mass_f = total_mass_o - H_mass_f
  mass_ratio = nonH_mass_f / nonH_mass_o

  for atom in molecule.atoms():
    index = atom.index
    if atom.element == openmm.app.element.hydrogen:
      omm_system.setParticleMass(index, H_mass * openmm.unit.dalton)
    else:
      atom_mass = omm_system.getParticleMass(atom.index)
      omm_system.setParticleMass(index, mass_ratio * atom_mass)
  return omm_system

