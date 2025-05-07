import os
import MMTK
from AlGDock.HMR import hydrogen_mass_repartitioning
H_mass = 4
ligand_database = 'ligand.db'

MMTK.Database.molecule_types.directory = os.path.dirname(os.path.abspath(ligand_database))
molecule = MMTK.Molecule(os.path.basename(ligand_database))
molecule = hydrogen_mass_repartitioning(molecule, H_mass)
MMTK_atom_masses = {}
for atom in molecule.atoms:
    MMTK_atom_masses[atom.name.split('i')[0]] = round(atom.mass(), 2)

#----------------------------------------------------------
import numpy as np
import openmm
from openmm.app import AmberPrmtopFile, NoCutoff
from AlGDock.HMR import hydrogen_mass_repartitioning_openmm
H_mass = 4

prmtopL = AmberPrmtopFile("ligand.prmtop")
OMM_systemL = prmtopL.createSystem(nonbondedMethod=NoCutoff, constraints=None)
molecule = prmtopL.topology
OMM_systemL = hydrogen_mass_repartitioning_openmm(molecule, OMM_systemL, H_mass)
OpenMM_atom_masses = {}
for atom in molecule.atoms():
    atom_mass = OMM_systemL.getParticleMass(atom.index).value_in_unit(openmm.unit.dalton)
    OpenMM_atom_masses[atom.name]= round(atom_mass,2)

assert OpenMM_atom_masses == MMTK_atom_masses

#=============================================================
# validate self.L_first_atom in topology.py

import AlGDock.IO
IO_prmtop = AlGDock.IO.prmtop()
prmtop_R = IO_prmtop.read('receptor.prmtop')
prmtop_RL = IO_prmtop.read('complex.prmtop')
ligand_ind = [ind for ind in range(len(prmtop_RL['RESIDUE_LABEL']))
  if prmtop_RL['RESIDUE_LABEL'][ind] not in prmtop_R['RESIDUE_LABEL']]
L_first_atom = prmtop_RL['RESIDUE_POINTER'][ligand_ind[0]] - 1

#---------------------
import openmm
from openmm.app import AmberPrmtopFile
prmtopR = AmberPrmtopFile("receptor.prmtop")
prmtopRL = AmberPrmtopFile("complex.prmtop")
receptor_atoms = set(atom.index for atom in prmtopR.topology.atoms())
complex_atoms = set(atom.index for atom in prmtopRL.topology.atoms())
ligand_atoms = complex_atoms - receptor_atoms
if ligand_atoms:
  L_first_atom = min(ligand_atoms) #5011

