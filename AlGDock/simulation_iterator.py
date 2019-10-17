import MMTK
from MMTK.ParticleProperties import Configuration

import numpy as np

import time


class SimulationIterator:
  """SimulationIterators take a molecular configuration and generate a new one.

  ...

  Attributes
  ----------
  args : AlGDock.simulation_arguments.SimulationArguments
    Simulation arguments
  top : AlGDock.topology.Topology
    Topology with ligand
  system : AlGDock.system.System
    System
  _sampler : dict
  """
  def __init__(self, args, top, system):
    """Initializes the class

    Parameters
    ----------
    args : AlGDock.simulation_arguments.SimulationArguments
      Simulation arguments
    top : AlGDock.topology.Topology
      Topology with ligand
    system : AlGDock.system.System
      System
    """
    self.args = args
    self.top = top
    self.system = system

    self._samplers = {}
    # Uses cython class
    # from SmartDarting import SmartDartingIntegrator # @UnresolvedImport
    # Uses python class
    from AlGDock.Integrators.SmartDarting.SmartDarting \
      import SmartDartingIntegrator # @UnresolvedImport
    self._samplers['BC_SmartDarting'] = SmartDartingIntegrator(\
      self.top.universe, self.top.molecule, False)
    self._samplers['CD_SmartDarting'] = SmartDartingIntegrator(\
      self.top.universe, self.top.molecule, True)

    from AlGDock.Integrators.ExternalMC.ExternalMC import ExternalMCIntegrator
    self._samplers['ExternalMC'] = ExternalMCIntegrator(\
      self.top.universe, self.top.molecule, step_size=0.25*MMTK.Units.Ang)

    for p in ['BC', 'CD']:
      if self.args.params[p]['sampler'] == 'HMC':
        from AlGDock.Integrators.HamiltonianMonteCarlo.HamiltonianMonteCarlo \
          import HamiltonianMonteCarloIntegrator
        self._samplers[p] = HamiltonianMonteCarloIntegrator(self.top.universe)
      elif self.args.params[p]['sampler'] == 'NUTS':
        from NUTS import NUTSIntegrator  # @UnresolvedImport
        self._samplers[p] = NUTSIntegrator(self.top.universe)
      elif self.args.params[p]['sampler'] == 'VV':
        from AlGDock.Integrators.VelocityVerlet.VelocityVerlet \
          import VelocityVerletIntegrator
        self._samplers[p] = VelocityVerletIntegrator(self.top.universe)
      else:
        raise Exception('Unrecognized sampler!')

  def iteration(self, seed, process, params_k, \
      initialize=False, reference=0):
    """Performs an iteration for a single thermodynamic state

    Parameters
    ----------
    seed : np.array
      Starting configuration
    process : str
      Process, either 'BC' or 'CD'
    params_k : dict of float
      Parameters describing a thermodynamic state
    """

    self.top.universe.setConfiguration(Configuration(self.top.universe, seed))

    self.system.setParams(params_k)
    if 'delta_t' in params_k.keys():
      delta_t = params_k['delta_t']
    else:
      raise Exception('No time step specified')
    if 'steps_per_trial' in params_k.keys():
      steps_per_trial = params_k['steps_per_trial']
    else:
      steps_per_trial = self.args.params[process]['steps_per_sweep']

    if initialize:
      steps = self.args.params[process]['steps_per_seed']
      ndarts = self.args.params[process]['darts_per_seed']
    else:
      steps = self.args.params[process]['steps_per_sweep']
      ndarts = self.args.params[process]['darts_per_sweep']

    random_seed = reference * reference + int(abs(seed[0][0] * 10000))
    if self.args.random_seed > 0:
      random_seed += self.args.random_seed
    else:
      random_seed += int(time.time() * 1000)

    results = {}

    # Execute external MCMC moves
    if (process == 'CD') and (self.args.params['CD']['MCMC_moves']>0) \
        and (params_k['alpha'] < 0.1) and (self.args.params['CD']['pose']==-1):
      time_start_ExternalMC = time.time()
      dat = self._samplers['ExternalMC'](ntrials=5, T=params_k['T'])
      results['acc_ExternalMC'] = dat[2]
      results['att_ExternalMC'] = dat[3]
      results['time_ExternalMC'] = (time.time() - time_start_ExternalMC)

    # Execute dynamics sampler
    time_start_sampler = time.time()
    dat = self._samplers[process](\
      steps=steps, steps_per_trial=steps_per_trial, \
      T=params_k['T'], delta_t=delta_t, \
      normalize=(process=='BC'), adapt=initialize, random_seed=random_seed)
    results['acc_Sampler'] = dat[2]
    results['att_Sampler'] = dat[3]
    results['delta_t'] = dat[4]
    results['time_Sampler'] = (time.time() - time_start_sampler)

    # Execute smart darting
    if (ndarts > 0) and not ((process == 'CD') and (params_k['alpha'] < 0.1)):
      time_start_SmartDarting = time.time()
      dat = self._samplers[process+'_SmartDarting'](\
        ntrials=ndarts, T=params_k['T'], random_seed=random_seed+5)
      results['acc_SmartDarting'] = dat[2]
      results['att_SmartDarting'] = dat[3]
      results['time_SmartDarting'] = (time.time() - time_start_SmartDarting)

    # Store and return results
    results['confs'] = np.copy(dat[0][-1])
    results['Etot'] = dat[1][-1]
    results['reference'] = reference

    return results

  def iteration_worker(self, input, output):
    """Executes an iteration from a multiprocessing queue

    Parameters
    ----------
    input : multiprocessing.Queue
      Tasks to complete
    output : multiprocessing.Queue
      Completed tasks
    """
    for args in iter(input.get, 'STOP'):
      result = self.iteration(*args)
      output.put(result)

  def initializeSmartDartingConfigurations(self, seeds, process, log, data):
    """Initializes the configurations for Smart Darting

    Parameters
    ----------
    seeds : list of np.array
      Starting configurations
    process : str
      Process, either 'BC' or 'CD'
    log : AlGDock.logger.Logger
      Logger that includes tee function
    data : AlGDock.simulation_data.SimulationData
      Location for minimized configurations
    """
    if self.args.params[process]['darts_per_seed'] > 0:
      outstr = self._samplers[process + '_SmartDarting'].set_confs(seeds)
      data[process].confs['SmartDarting'] = \
        self._samplers[process+'_SmartDarting'].confs
      log.tee(outstr)

  def addSmartDartingConfigurations(self, new_confs, process, log, data):
    """Adds new configurations for Smart Darting

    Parameters
    ----------
    new_confs : list of np.array
      New configurations
    process : str
      Process, either 'BC' or 'CD'
    data : AlGDock.simulation_data.SimulationData object
      Location for minimized configurations
    """
    if self.args.params[process]['darts_per_seed'] > 0:
      confs_SmartDarting = [np.copy(conf) \
        for conf in data[process].confs['samples'][k][-1]]
      outstr = self._samplers[process+'_SmartDarting'].set_confs(\
        new_confs + data[process].confs['SmartDarting'])
      data[process].confs['SmartDarting'] = \
        self._samplers[process+'_SmartDarting'].confs
      log.tee(outstr)

  def clearSmartDartingConfigurations(self, process):
    """Clears the list of configurations for Smart Darting

    Parameters
    ----------
    seeds : list of np.array
      Starting configurations
    process : str
      Process, either 'BC' or 'CD'
    """
    if self.args.params[process]['darts_per_seed'] > 0:
      self._samplers[process + '_SmartDarting'].confs = []
