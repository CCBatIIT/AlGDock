import os, sys
import numpy as np

try:
  import MMTK
  from MMTK.ParticleProperties import Configuration
  from MMTK.ForceFields import ForceField
  from MMTK.ForceFields import Amber12SBForceField
except ImportError:
  MMTK = None

try:
  import openmm
  from openmm.app import AmberPrmtopFile, AmberInpcrdFile, Simulation, NoCutoff
  import openmm.unit as unit

except ImportError:
  openmm = None

from AlGDock.logger import NullDevice

class TopologyMMTK:
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
    if includeReceptor and \
        (args.FNs['receptor_database'] is not None) and \
        os.path.isfile(args.FNs['receptor_database']):
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
      self.inv_prmtop_atom_order_L[self.prmtop_atom_order_L[i]] = i

    # Create universe and add molecule to universe
    self.universe = MMTK.Universe.InfiniteUniverse()
    self.universe.addObject(self.molecule)
    if includeReceptor:
      if self.molecule_R is not None:
        self.universe.addObject(self.molecule_R)
      else:
        self.universe = None

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

class TopologyUsingOpenMM:
  """Describes the system to simulate using OpenMM
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
    original_stderr = sys.stderr
    sys.stderr = NullDevice()
    prmtopL = AmberPrmtopFile(self.args.FNs['prmtop']["L"])
    inpcrdL = AmberInpcrdFile(self.args.FNs['inpcrd']["L"])
    self.OMM_systemL = prmtopL.createSystem(nonbondedMethod=NoCutoff, constraints=None)
    self.molecule = prmtopL.topology

    if includeReceptor:
      prmtopR = AmberPrmtopFile(self.args.FNs['prmtop']["R"])
      inpcrdR = AmberInpcrdFile(self.args.FNs['inpcrd']["R"])
      self.OMM_systemR = prmtopR.createSystem(nonbondedMethod=NoCutoff, constraints=None)
      self.molecule_R = prmtopR.topology
    else:
      self.molecule_R = None

    sys.stderr = original_stderr
    dummy_integrator = openmm.LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)

    # Hydrogen Mass Repartitioning
    # (sets hydrogen mass to H_mass and scales other masses down)
    if args.params['BC']['H_mass'] > 0.:
      from AlGDock.HMR import hydrogen_mass_repartitioning_openmm
      self.OMM_systemL = hydrogen_mass_repartitioning_openmm(self.molecule,  self.OMM_systemL,\
        args.params['BC']['H_mass'])
    self.OMM_simulaitonL = openmm.app.Simulation(prmtopL.topology, self.OMM_systemL, dummy_integrator)
    self.OMM_simulaitonL.context.setPositions(inpcrdL.positions)
    self.inv_prmtop_atom_order_L = self.prmtop_atom_order_L = np.array([atom.index \
      for atom in self.molecule.atoms()], dtype=int)

    if includeReceptor:
      self.OMM_simulaitonR = openmm.app.Simulation(prmtopR.topology, self.OMM_systemR, dummy_integrator)
      self.OMM_simulaitonR.context.setPositions(inpcrdR.positions)

      if (args.FNs['prmtop']['R'] is not None) and \
         (args.FNs['prmtop']['RL'] is not None):
        prmtopRL = AmberPrmtopFile(self.args.FNs['prmtop']["RL"])
        receptor_atoms = set(atom.index for atom in prmtopR.topology.atoms())
        complex_atoms = set(atom.index for atom in prmtopRL.topology.atoms())
        ligand_atoms = complex_atoms - receptor_atoms
        if ligand_atoms:
          self.L_first_atom = min(ligand_atoms)
        else:
          print("No ligand atoms found in complex prmtop file.")
          self.L_first_atom = 0
      else:
        self.L_first_atom = 0
    else:
      self.L_first_atom = 0

