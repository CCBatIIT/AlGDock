import os, sys
import numpy as np

import MMTK

from AlGDock.logger import NullDevice

class Topology:
  """Describes the system to simulate

  ...

  Attributes
  ----------
  molecule : MMTK.Molecule
    Ligand molecule, like an OpenMM Chain
  molecule_R : MMTK.Molecule
    If includeReceptor, Receptor molecule, like an OpenMM Chain
  universe : MMTK.InfiniteUniverse
    A universe containing all the molecules

  L_first_atom : int
    Index of the first ligand atom in the AMBER prmtop file
  inv_prmtop_atom_order_L : np.array
    Indices to convert from prmtop to MMTK ordering
  prmtop_atom_order_L : np.array
    Indices to convert from prmtop to MMTK ordering

  """
  def __init__(self, args, includeReceptor=False):
    """Initializes the class

    Parameters
    ----------
    args : simulation_arguments.SimulationArguments
      Simulation arguments
    includeReceptor : bool
      Includes the receptor in the topology
    """

    # Set up the system
    original_stderr = sys.stderr
    sys.stderr = NullDevice()
    MMTK.Database.molecule_types.directory = \
      os.path.dirname(args.FNs['ligand_database'])
    self.molecule = MMTK.Molecule(\
      os.path.basename(args.FNs['ligand_database']))
    if includeReceptor:
      self.molecule_R = MMTK.Molecule(\
        os.path.basename(args.FNs['receptor_database']))
    else:
      self.molecule_R = None
    sys.stderr = original_stderr

    # Hydrogen Mass Repartitioning
    # (sets hydrogen mass to H_mass and scales other masses down)
    if args.params['BC']['H_mass'] > 0.:
      from AlGDock.HMR import hydrogen_mass_repartitioning
      self.molecule = hydrogen_mass_repartitioning(self.molecule, \
        args.params['BC']['H_mass'])

    # # Helpful variables for referencing and indexing atoms in the molecule
    # self.molecule.heavy_atoms = [ind for (atm,ind) in \
    #   zip(self.molecule.atoms,range(self.molecule.numberOfAtoms())) \
    #   if atm.type.name!='hydrogen']
    # self.molecule.nhatoms = len(self.molecule.heavy_atoms)

    self.prmtop_atom_order_L = np.array([atom.number \
      for atom in self.molecule.prmtop_order], dtype=int)
    self.inv_prmtop_atom_order_L = \
      np.zeros(shape=len(self.prmtop_atom_order_L), dtype=int)
    for i in range(len(self.prmtop_atom_order_L)):
      self.inv_prmtop_atom_order_L[
          self.prmtop_atom_order_L[i]] = i

    # Create universe and add molecule to universe
    self.universe = MMTK.Universe.InfiniteUniverse()
    self.universe.addObject(self.molecule)
    if includeReceptor:
      self.universe.addObject(self.molecule_R)

    # Define L_first_atom
    if includeReceptor:
      if (args.FNs['prmtop']['R'] is not None) and \
         (args.FNs['prmtop']['RL'] is not None):
        import AlGDock.IO
        IO_prmtop = AlGDock.IO.prmtop()
        prmtop_R = IO_prmtop.read(args.FNs['prmtop']['R'])
        prmtop_RL = IO_prmtop.read(args.FNs['prmtop']['RL'])
        ligand_ind = [
            ind for ind in range(len(prmtop_RL['RESIDUE_LABEL']))
            if prmtop_RL['RESIDUE_LABEL'][ind] not in prmtop_R['RESIDUE_LABEL']
        ]
        if len(ligand_ind) == 0:
          raise Exception('Ligand not found in complex prmtop')
        elif len(ligand_ind) > 1:
          print '  possible ligand residue labels: '+\
            ', '.join([prmtop_RL['RESIDUE_LABEL'][ind] for ind in ligand_ind])
        print 'ligand residue name: ' + \
          prmtop_RL['RESIDUE_LABEL'][ligand_ind[0]].strip()
        self.L_first_atom = prmtop_RL['RESIDUE_POINTER'][ligand_ind[0]] - 1
      else:
        self.L_first_atom = 0
    else:
      self.L_first_atom = 0
