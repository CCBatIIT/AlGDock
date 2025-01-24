import os, sys
import numpy as np
import copy

try:
  import MMTK
  from MMTK.ParticleProperties import Configuration
  from MMTK.ForceFields import ForceField
  from MMTK.ForceFields import Amber12SBForceField
except ImportError:
  MMTK = None

try:
  import openmm
  import openmm.unit as unit
  from openmm.app import AmberPrmtopFile, AmberInpcrdFile, Simulation, NoCutoff
  from openmm import *
except ImportError:
  openmm = None

from AlGDock.BindingPMF import scalables
from AlGDock.BindingPMF import HMStime
import AlGDock.IO
prmtop_IO = AlGDock.IO.prmtop()
varnames = ['POINTERS','TITLE','ATOM_NAME','AMBER_ATOM_TYPE','CHARGE','MASS',\
            'NONBONDED_PARM_INDEX','LENNARD_JONES_ACOEF','LENNARD_JONES_BCOEF',\
            'ATOM_TYPE_INDEX','BONDS_INC_HYDROGEN','BONDS_WITHOUT_HYDROGEN',\
            'RADII','SCREEN']

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

class System:
  """Forces and masses

  ...

  Attributes
  ----------
  _evaluators : dict
  _forceFields : dict
  T : float
    Current temperature
  args : AlGDock.simulation_arguments.SimulationArguments
    Simulation arguments
  log : AlGDock.logger.Logger
    Simulation log
  top : AlGDock.topology.Topology
    Topology of the ligand
  top_RL : AlGDock.topology.Topology
    Topology of the complex
  starting_pose : numpy.array
    Starting pose used for simulations with an external restraint
  """
  def __init__(self, args, log, top, top_RL=None, starting_pose=None):
    """Initializes the class

    Parameters
    ----------
    args : AlGDock.simulation_arguments.SimulationArguments
      Simulation arguments
    log : AlGDock.logger.Logger
      Simulation log
    top : AlGDock.topology.Topology
      Topology of the ligand
    top_RL : AlGDock.topology.Topology
      Topology of the complex
    starting_pose : numpy.array
      Starting pose used for simulations with an external restraint
    """
    self.args = args
    self.log = log
    self.top = top   # old-> AlGDock.topology.TopologyMMTK, now: AlGDock.topology.TopologyUsingOpenMM
    self.top_RL = top_RL
    self.starting_pose = starting_pose

    self._evaluators = {}
    self._forceFields = {}

    # Molecular mechanics force fields
    if MMTK:
      self._forceFields['gaff'] = Amber12SBForceField(
      parameter_file=self.args.FNs['forcefield'],
      mod_files=self.args.FNs['frcmodList'])
    else:
      self._forceFields['gaff'] = None # Use ambertools to generate topology files, which include GAFF parameters.

  def setParams(self, params):
    """Sets the universe evaluator to values appropriate for the parameters.

    Parameters
    ----------
    params : dict
      The elements in the dictionary params can be:
        MM - True, to turn on the Generalized AMBER force field
        site - True, to turn on the binding site
        sLJr - scaling of the soft Lennard-Jones repulsive grid
        sLJa - scaling of the soft Lennard-Jones attractive grid
        sELE - scaling of the soft electrostatic grid
        LJr - scaling of the Lennard-Jones repulsive grid
        LJa - scaling of the Lennard-Jones attractive grid
        ELE - scaling of the electrostatic grid
        k_angular_int - spring constant of flat-bottom wells for angular internal degrees of freedom (kJ/nm)
        k_spatial_ext - spring constant of flat-bottom wells for spatial external degrees of freedom (kJ/nm)
        k_angular_ext - spring constant of flat-bottom wells for angular external degrees of freedom (kJ/nm)
        T - the temperature in K
    """

    self.T = params['T']

    # Reuse evaluators that have been stored
    evaluator_key = ','.join(['%s:%s'%(k,params[k]) \
      for k in sorted(params.keys())])
    if MMTK:
      if evaluator_key in self._evaluators.keys():
        self.top.universe._evaluator[(None,None,None)] = \
          self._evaluators[evaluator_key]
    # In MMTK, these evaluators were associated with the universe (equal to OpenMM system)
    return

    # Otherwise create a new evaluator
    fflist = []
    if ('MM' in params.keys()) and params['MM']:
      fflist.append(self._forceFields['gaff'])
    if ('site' in params.keys()) and params['site']:
      # Set up the binding site in the force field
      append_site = True  # Whether to append the binding site to the force field
      if not 'site' in self._forceFields.keys():
        print(self.args.params['CD']['site'], self.args.params['CD']['site_center'], self.args.params['CD']['site_max_R'])
        if (self.args.params['CD']['site']=='Sphere') and \
           (self.args.params['CD']['site_center'] is not None) and \
           (self.args.params['CD']['site_max_R'] is not None):
          from AlGDock.ForceFields.Sphere.Sphere import SphereForceField
          self._forceFields['site'] = SphereForceField(
            center=self.args.params['CD']['site_center'],
            max_R=self.args.params['CD']['site_max_R'],
            name='site')
        elif (self.args.params['CD']['site']=='Cylinder') and \
             (self.args.params['CD']['site_center'] is not None) and \
             (self.args.params['CD']['site_direction'] is not None):
          from AlGDock.ForceFields.Cylinder.Cylinder import CylinderForceField
          self._forceFields['site'] = CylinderForceField(
            origin=self.args.params['CD']['site_center'],
            direction=self.args.params['CD']['site_direction'],
            max_Z=self.args.params['CD']['site_max_Z'],
            max_R=self.args.params['CD']['site_max_R'],
            name='site')
        else:
          # Do not append the site if it is not defined
          print('Binding site not defined!')
          append_site = False
      if append_site:
        fflist.append(self._forceFields['site'])

    # Add scalable terms
    for scalable in scalables:
      if (scalable in params.keys()) and params[scalable] > 0:
        # Load the force field if it has not been loaded
        if not scalable in self._forceFields.keys():
          if scalable == 'OBC':
            from AlGDock.ForceFields.OBC.OBC import OBCForceField
            if self.args.params['CD']['solvation']=='Fractional' and \
                ('ELE' in params.keys()):
              self.log.recordStart('grid_loading')
              self._forceFields['OBC'] = OBCForceField(\
                desolvationGridFN=self.args.FNs['grids']['desolv'])
              self.log.tee('  %s grid loaded from %s in %s'%(scalable, \
                os.path.basename(self.args.FNs['grids']['desolv']), \
                HMStime(self.log.timeSince('grid_loading'))))
            else:
              self._forceFields['OBC'] = OBCForceField()
          else:  # Grids
            self.log.recordStart('grid_loading')
            grid_FN = self.args.FNs['grids'][{
              'sLJr': 'LJr',
              'sLJa': 'LJa',
              'sELE': 'ELE',
              'LJr': 'LJr',
              'LJa': 'LJa',
              'ELE': 'ELE'
            }[scalable]]
            grid_scaling_factor = 'scaling_factor_' + \
              {'sLJr':'LJr','sLJa':'LJa','sELE':'electrostatic', \
               'LJr':'LJr','LJa':'LJa','ELE':'electrostatic'}[scalable]

            # Determine the grid threshold
            if scalable == 'sLJr':
              grid_thresh = 10.0
            elif scalable == 'sELE':
              # The maximum value is set so that the electrostatic energy
              # less than or equal to the Lennard-Jones repulsive energy
              # for every heavy atom at every grid point

              if MMTK:
                scaling_factors_ELE = np.array([ \
                self.top.molecule.getAtomProperty(a, 'scaling_factor_electrostatic') \
                  for a in self.top.molecule.atomList()],dtype=float)
                scaling_factors_LJr = np.array([ \
                self.top.molecule.getAtomProperty(a, 'scaling_factor_LJr') \
                  for a in self.top.molecule.atomList()],dtype=float)
              else:
                """
                get the scaling_factors_ELE and scaling_factors_LJr using openmm instead of MMTK.
                Note that MMTK may provide the values in a different atom order compared to OpenMM.
                """
                scaling_factors_ELE = []
                scaling_factors_LJr = []
                ligand_prmtop = self.args.FNs['prmtop']['L']
                prmtop = AmberPrmtopFile(ligand_prmtop)
                topology = prmtop.topology

                prmtop_alg = prmtop_IO.read(ligand_prmtop, varnames)
                NATOM = prmtop_alg['POINTERS'][0]
                NTYPES = prmtop_alg['POINTERS'][1]
                LJ_radius = np.ndarray(shape=(NTYPES), dtype=float)
                LJ_depth = np.ndarray(shape=(NTYPES), dtype=float)
                for i in range(NTYPES):
                  LJ_index = prmtop_alg['NONBONDED_PARM_INDEX'][NTYPES * i + i] - 1
                  if prmtop_alg['LENNARD_JONES_ACOEF'][LJ_index] < 1.0e-6:
                    LJ_radius[i] = 0
                    LJ_depth[i] = 0
                  else:
                    factor = 2 * prmtop_alg['LENNARD_JONES_ACOEF'][LJ_index] / prmtop_alg['LENNARD_JONES_BCOEF'][
                      LJ_index]
                    LJ_radius[i] = pow(factor, 1.0 / 6.0) * 0.5
                    LJ_depth[i] = prmtop_alg['LENNARD_JONES_BCOEF'][LJ_index] / 2 / factor
                root_LJ_depth = np.sqrt(LJ_depth)
                LJ_diameter = LJ_radius * 2
                atom_type_indicies = [prmtop_alg['ATOM_TYPE_INDEX'][atom_index] - 1 for atom_index in range(NATOM)]
                scaling_factor_LJr_dict = dict()
                for (name, type_index) in zip(prmtop_alg['ATOM_NAME'], atom_type_indicies):
                  scaling_factor_LJr_dict[name.strip()] = round(
                    4.184 * root_LJ_depth[type_index] * (LJ_diameter[type_index] ** 6), 6)

                atoms_name = [a.name for a in topology.atoms()]
                atoms_elements = dict()
                for a in topology.atoms():
                  atoms_elements[a.name] = a.element.symbol
                system = prmtop.createSystem(nonbondedMethod=NoCutoff, nonbondedCutoff=1.0 * nanometer,
                                             constraints=None)
                force = system.getForce(3)
                assert force.getName() == 'NonbondedForce'
                for i in range(system.getNumParticles()):
                  atom_name = atoms_name[i]
                  charge, sigma, epsilon = force.getParticleParameters(i)
                  charge_val = charge.value_in_unit(elementary_charge)
                  epsilon = epsilon.value_in_unit(kilojoule / mole)
                  if epsilon == 0:
                    continue
                  scaling_factors_ELE.append(round(charge_val * 4.184, 6))
                  scaling_factors_LJr.append(scaling_factor_LJr_dict[atom_name])
                scaling_factors_ELE = np.array(scaling_factors_ELE)
                scaling_factors_LJr = np.array(scaling_factors_LJr)

              #---------
              toKeep = np.logical_and(scaling_factors_LJr > 10.,
                                      abs(scaling_factors_ELE) > 0.1)
            
              scaling_factors_ELE = scaling_factors_ELE[toKeep]
              scaling_factors_LJr = scaling_factors_LJr[toKeep]
            
              grid_thresh = min(
                abs(scaling_factors_LJr * 10.0 / scaling_factors_ELE))
            else:
              grid_thresh = -1  # There is no threshold for grid points

            from AlGDock.ForceFields.Grid.Interpolation \
              import InterpolationForceField
            self._forceFields[scalable] = InterpolationForceField(grid_FN, \
              name=scalable, interpolation_type='Trilinear', \
              strength=params[scalable], scaling_property=grid_scaling_factor,
              inv_power=4 if scalable=='LJr' else None, \
              grid_thresh=grid_thresh)
            self.log.tee('  %s grid loaded from %s in %s'%(scalable, \
              os.path.basename(grid_FN), \
              HMStime(self.log.timeSince('grid_loading'))))

        # Set the force field strength to the desired value
        self._forceFields[scalable].set_strength(params[scalable])
        fflist.append(self._forceFields[scalable])

    if ('k_angular_int' in params.keys()) or \
       ('k_spatial_ext' in params.keys()) or \
       ('k_angular_ext' in params.keys()):

      # Load the force field if it has not been loaded
      if MMTK:
        if not ('ExternalRestraint' in self._forceFields.keys()):
          Xo = np.copy(self.top.universe.configuration().array)
          self.top.universe.setConfiguration(
            Configuration(self.top.universe, self.starting_pose))
          import AlGDock.rigid_bodies
          rb = AlGDock.rigid_bodies.identifier(self.top.universe,
                                               self.top.molecule)
          (TorsionRestraintSpecs, ExternalRestraintSpecs) = rb.poseInp()
          self.top.universe.setConfiguration(Configuration(
            self.top.universe, Xo))
      else:
        if not ('ExternalRestraint' in self._forceFields.keys()):
          import AlGDock.rigid_bodies
          # TODO, modify rigid_bodies using openmm
          # rb = AlGDock.rigid_bodies.identifier(self.top.universe,
          #                                      self.top.molecule)
          self.top.OMM_simulation.context.setPositions(self.starting_pose)

        # Create force fields
        from AlGDock.ForceFields.Pose.PoseFF import InternalRestraintForceField
        self._forceFields['InternalRestraint'] = \
          InternalRestraintForceField(TorsionRestraintSpecs)
        from AlGDock.ForceFields.Pose.PoseFF import ExternalRestraintForceField
        self._forceFields['ExternalRestraint'] = \
          ExternalRestraintForceField(*ExternalRestraintSpecs)

      # Set parameter values
      if ('k_angular_int' in params.keys()):
        self._forceFields['InternalRestraint'].set_k(\
          params['k_angular_int'])
        fflist.append(self._forceFields['InternalRestraint'])

      if ('k_spatial_ext' in params.keys()):
        self._forceFields['ExternalRestraint'].set_k_spatial(\
          params['k_spatial_ext'])
        fflist.append(self._forceFields['ExternalRestraint'])

      if ('k_angular_ext' in params.keys()):
        self._forceFields['ExternalRestraint'].set_k_angular(\
          params['k_angular_ext'])

    if MMTK:
      compoundFF = fflist[0]
      for ff in fflist[1:]:
        compoundFF += ff
      self.top.universe.setForceField(compoundFF)

      if self.top_RL.universe is not None:
        if 'OBC_RL' in params.keys():
          if not 'OBC_RL' in self._forceFields.keys():
            from AlGDock.ForceFields.OBC.OBC import OBCForceField
            self._forceFields['OBC_RL'] = OBCForceField()
          self._forceFields['OBC_RL'].set_strength(params['OBC_RL'])
          if (params['OBC_RL'] > 0):
            self.top_RL.universe.setForceField(self._forceFields['OBC_RL'])

      eval = ForceField.EnergyEvaluator(\
        self.top.universe, self.top.universe._forcefield, None, None, None, None)
      eval.key = evaluator_key
      self.top.universe._evaluator[(None, None, None)] = eval
      self._evaluators[evaluator_key] = eval
    #TODO: need to create EnergyEvaluator using openmm
    else:
      # e.g. fflist=[self._forceFields['OBC']]
      # self._forceFields['OBC'] = OBCForceField()
      compoundFF = fflist[0]
      for ff in fflist[1:]:
        compoundFF += ff

  def energyTerms(self, confs, E=None, process='CD'):
    """Calculates energy terms for a series of configurations

    Units are kJ/mol.

    Parameters
    ----------
    confs : list of np.array
      Configurations
    E : dict of np.array
      Dictionary to add to
    process : str
      Process, either 'BC' or 'CD'

    Returns
    -------
    E : dict of np.array
      Dictionary of energy terms
    """
    if E is None:
      E = {}

    params_full = self.paramsFromAlpha(alpha=1.0,
                                               process=process,
                                               site=(process == 'CD'))
    if process == 'CD':
      for scalable in scalables:
        params_full[scalable] = 1
    self.setParams(params_full)

    # Molecular mechanics and grid interaction energies
    E['MM'] = np.zeros(len(confs), dtype=float)
    if process == 'BC':
      if 'OBC' in params_full.keys():
        E['OBC'] = np.zeros(len(confs), dtype=float)
    if process == 'CD':
      for term in (scalables):
        E[term] = np.zeros(len(confs), dtype=float)
      if self.isForce('site'):
        E['site'] = np.zeros(len(confs), dtype=float)
      if self.isForce('InternalRestraint'):
        E['k_angular_int'] = np.zeros(len(confs), dtype=float)
      if self.isForce('ExternalRestraint'):
        E['k_angular_ext'] = np.zeros(len(confs), dtype=float)
        E['k_spatial_ext'] = np.zeros(len(confs), dtype=float)
    if MMTK:
      for c in range(len(confs)):
        self.top.universe.setConfiguration(
          Configuration(self.top.universe, confs[c]))
        eT = self.top.universe.energyTerms()
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
      return E
    else:
      for c in range(len(confs)):
        self.top.OMM_simulation.context.setPositions(confs[c])
        state = self.top.OMM_simulation.context.getState(getEnergy=True)
        force_groups = self.top.OMM_system.getNumForces()
        for i in range(force_groups):
          self.top.OMM_system.getForce(i).setForceGroup(i)
          group_state = self.top.OMM_simulation.context.getState(getEnergy=True, groups={i})
          group_energy = group_state.getPotentialEnergy()
          force_name = self.top.OMM_system.getForce(i).__class__.__name__
          try:
            E[term_map[force_name]][c] += group_energy
          except KeyError:
            print(force_name)
            print('Keys in term map', term_map.keys())
            print('Keys in E', E.keys())
            raise Exception('key not found in term map or E')
      return E

  def paramsFromAlpha(self,
                       alpha,
                       process='CD',
                       params_o=None,
                       site=True,
                       crossed=False):
    """Creates a parameter dictionary for a given alpha value

    Parameters
    ----------
    alpha : float
      Progress variable for the thermodynamic state
    process : str
      Process, either 'BC' or 'CD'
    params_o : dict of float
      Parameter dictionary to modify
    site : bool
      If True, the ligand is restricted to the binding site
    crossed : bool
      If True, the initial protocol is complete

    Returns
    -------
    params : dict of float
      Parameter dictionary that corresponds to the given alpha
    """
    if (params_o is not None):
      params = copy.deepcopy(params_o)
      if 'steps_per_trial' not in params_o.keys():
        params[
          'steps_per_trial'] = 1 * self.args.params[process]['steps_per_sweep']
    else:
      params = {}
      params[
        'steps_per_trial'] = 1 * self.args.params[process]['steps_per_sweep']

    params['MM'] = True
    if crossed is not None:
      params['crossed'] = crossed

    if process == 'CD':
      alpha_g = 4. * (alpha - 0.5)**2 / (1 + np.exp(-100 * (alpha - 0.5)))
      if alpha_g < 1E-10:
        alpha_g = 0
      if self.args.params['CD']['solvation']=='Desolvated' or \
         self.args.params['CD']['solvation']=='Reduced':
        params['OBC'] = 0
      elif self.args.params['CD']['solvation'] == 'Fractional':
        params['OBC'] = alpha_g  # Scales the solvent with the grid
      elif self.args.params['CD']['solvation'] == 'Full':
        params['OBC'] = 1.0
      if self.args.params['CD']['pose'] > -1:
        # Pose BPMF
        alpha_r = np.tanh(16 * alpha * alpha)
        params['alpha'] = alpha
        params['k_angular_int'] = self.args.params['CD']['k_pose'] * alpha_r
        params['k_angular_ext'] = self.args.params['CD']['k_pose']
        params['k_spatial_ext'] = self.args.params['CD']['k_pose']
        params['sLJr'] = alpha_g
        params['sLJa'] = alpha_g
        params['ELE'] = alpha_g
        params['T'] = alpha_r * (self.args.params['BC']['T_TARGET'] - self.args.params['BC']['T_HIGH']) + self.args.params['BC']['T_HIGH']
      else:
        # BPMF
        alpha_sg = 1. - 4. * (alpha - 0.5)**2
        params['alpha'] = alpha
        params['sLJr'] = alpha_sg
        if self.args.params['CD']['solvation'] != 'Reduced':
          params['sELE'] = alpha_sg
        params['LJr'] = alpha_g
        params['LJa'] = alpha_g
        params['ELE'] = alpha_g \
          if self.args.params['CD']['solvation']!='Reduced' else 0.2*alpha_g
        if site is not None:
          params['site'] = site
        if self.args.params['CD']['temperature_scaling'] == 'Linear':
          params['T'] = alpha * (self.args.params['BC']['T_TARGET'] - self.args.params['BC']['T_HIGH']) + self.args.params['BC']['T_HIGH']
        elif self.args.params['CD']['temperature_scaling'] == 'Quadratic':
          params['T'] = alpha_g * (self.args.params['BC']['T_TARGET'] - self.args.params['BC']['T_HIGH']) + self.args.params['BC']['T_HIGH']
    elif process == 'BC':
      # If alpha = 0.0, T = T_HIGH. If alpha = 1.0, T = T_TARGET.
      params['alpha'] = alpha
      params['T'] = self.args.params['BC']['T_HIGH'] - alpha * (self.args.params['BC']['T_HIGH'] - self.args.params['BC']['T_TARGET'])
      if self.args.params['BC']['solvation'] == 'Desolvated':
        params['OBC'] = alpha
      elif self.args.params['BC']['solvation'] == 'Reduced':
        params['OBC'] = alpha
      elif self.args.params['BC']['solvation'] == 'Fractional':
        params['OBC'] = alpha
      elif self.args.params['BC']['solvation'] == 'Full':
        params['OBC'] = 1.0
    else:
      raise Exception("Unknown process!")

    return params

  def clear_evaluators(self):
    """Deletes the stored evaluators and grids to save memory
    """
    self._evaluators = {}
    for scalable in scalables:
      if (scalable in self._forceFields.keys()):
        del self._forceFields[scalable]

  def isForce(self, val):
    """Determines whether a force named 'val' is defined
    """
    return (val in self._forceFields.keys())

  def getGridParams(self):
    """Returns the counts, center, and spacing used for the electrostatic grid
    """
    self.setParams({'MM': True, 'ELE': 1})
    gd = self._forceFields['ELE'].grid_data
    dims = gd['counts']
    #TODO: what is this factor? Temporarily set it to 1.
    factor = 1
    center = factor * (gd['counts'] * gd['spacing'] / 2. + gd['origin'])
    spacing = factor * gd['spacing'][0]
    return (dims, center, spacing)
