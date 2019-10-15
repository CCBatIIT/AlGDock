import os, sys
import numpy as np

import copy

import MMTK
from MMTK.ForceFields import ForceField

from AlGDock.BindingPMF import scalables
from AlGDock.BindingPMF import R
from AlGDock.BindingPMF import HMStime


class System:
  """Forces and masses

  ...

  Attributes
  ----------
  _evaluators : dict
  _forceFields : dict
  T : float
    Current temperature
  T_HIGH : float
    High temperature
  T_TARGET : float
    Target temperature
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
    self.top = top
    self.top_RL = top_RL
    self.starting_pose = starting_pose

    self.T_HIGH = self.args.params['BC']['T_HIGH']
    self.T_TARGET = self.args.params['BC']['T_TARGET']

    self._evaluators = {}
    self._forceFields = {}

    # Molecular mechanics force fields
    from MMTK.ForceFields import Amber12SBForceField
    self._forceFields['gaff'] = Amber12SBForceField(
      parameter_file=self.args.FNs['forcefield'],
      mod_files=self.args.FNs['frcmodList'])

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
    if evaluator_key in self._evaluators.keys():
      self.top.universe._evaluator[(None,None,None)] = \
        self._evaluators[evaluator_key]
      return

    # Otherwise create a new evaluator
    fflist = []
    if ('MM' in params.keys()) and params['MM']:
      fflist.append(self._forceFields['gaff'])
    if ('site' in params.keys()) and params['site']:
      # Set up the binding site in the force field
      append_site = True  # Whether to append the binding site to the force field
      if not 'site' in self._forceFields.keys():
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
          print 'Binding site not defined!'
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
              # TODO: For conversion to OpenMM, get ParticleProperties from the Force object
              scaling_factors_ELE = np.array([ \
                self.top.molecule.getAtomProperty(a, 'scaling_factor_electrostatic') \
                  for a in self.top.molecule.atomList()],dtype=float)
              scaling_factors_LJr = np.array([ \
                self.top.molecule.getAtomProperty(a, 'scaling_factor_LJr') \
                  for a in self.top.molecule.atomList()],dtype=float)
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
        params['T'] = alpha_r * (self.T_TARGET - self.T_HIGH) + self.T_HIGH
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
          params['T'] = alpha * (self.T_TARGET - self.T_HIGH) + self.T_HIGH
        elif self.args.params['CD']['temperature_scaling'] == 'Quadratic':
          params['T'] = alpha_g * (self.T_TARGET - self.T_HIGH) + self.T_HIGH
    elif process == 'BC':
      # If alpha = 0.0, T = T_HIGH. If alpha = 1.0, T = T_TARGET.
      params['alpha'] = alpha
      params['T'] = self.T_HIGH - alpha * (self.T_HIGH - self.T_TARGET)
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
