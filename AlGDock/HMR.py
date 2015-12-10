import numpy as np

def hydrogen_mass_repartitioning(molecule, H_mass):
  H_mass_o = np.sum([atom.mass() for atom in molecule.atoms \
    if atom.type.name=='hydrogen'])
  H_mass_f = np.sum([1 for atom in molecule.atoms \
    if atom.type.name=='hydrogen'])*H_mass
  total_mass_o = np.sum([atom.mass() for atom in molecule.atoms])
  nonH_mass_o = np.sum([atom.mass() for atom in molecule.atoms \
    if atom.type.name!='hydrogen'])
  nonH_mass_f = total_mass_o - H_mass_f
  mass_ratio = nonH_mass_f/nonH_mass_o
  for atom in molecule.atoms:
    if atom.type.name=='hydrogen':
      atom.setMass(H_mass)
    else:
      atom.setMass(mass_ratio*atom.mass())
  return molecule