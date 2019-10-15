from AlGDock.BindingPMF import HMStime

import time
import numpy as np

import copy

import MMTK
import MMTK.Units
from MMTK.ParticleProperties import Configuration

from AlGDock.BindingPMF import R

class Initialization():
  """Establishes the initial protocol of thermodynamic states

  Attributes
  ----------
  args : AlGDock.simulation_arguments.SimulationArguments
    Simulation arguments
  log : AlGDock.logger.Logger
    Simulation log
  top : AlGDock.topology.Topology
    Topology of the ligand
  system : AlGDock.system.System
    Simulation system
  _get_confs_to_rescore : function
    Returns the configurations to rescore
  iterator : AlGDock.simulation_iterator.SimulationIterator
    Performs an iteration on one thermodynamic state
  data : AlGDock.simulation_data.SimulationData
    Stores results from the simulation
  """
  def __init__(self, args, log, top, system, iterator, data, save, _u_kln):
    """Initializes the class

    Parameters
    ----------
    args : AlGDock.simulation_arguments.SimulationArguments
      Simulation arguments
    log : AlGDock.logger.Logger
      Simulation log
    top : AlGDock.topology.Topology
      Topology of the ligand
    system : AlGDock.system.System
      Simulation system
    iterator : AlGDock.simulation_iterator.SimulationIterator
      Performs an iteration on one thermodynamic state
    data : AlGDock.simulation_data.SimulationData
      Stores results from the simulation
    save : AlGDock.BindingPMF.save
      Saves the data
    _u_kln : AlGDock.BindingPMF._u_kln
      Evaluates energies in different thermodynamic states
    """
    pass
    self.args = args
    self.log = log
    self.top = top
    self.system = system
    self.iterator = iterator
    self.data = data
    self.save = save
    self._u_kln = _u_kln

  def run(self, process, seeds):
    """Performs the initialization

    Parameters
    ----------
    process : str
      Process, either 'BC' or 'CD'
    seeds : list of np.array
      The starting configurations
    """
    self.log.set_lock(process)
    self.log.recordStart(process+' initialization')
    self.log.tee("\n>>> Initialization for %s, starting at "%process + \
      time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "\n")
    if process == 'BC':
      seeds = self._initialize_BC(seeds)
    elif process == 'CD':
      seeds = self._initialize_CD(seeds)
    self.log.tee("Elapsed time for initialization for %s:"%process + \
      HMStime(self.log.timeSince(process+' initialization')))
    self.log.clear_lock(process)

  def _initial_sim_state(self, seeds, process, params_k):
    """
    Initializes a state, returning the configurations and potential energy.
    Attempts simulation up to 12 times, adjusting the time step.
    """

    if not 'delta_t' in params_k.keys():
      params_k[
        'delta_t'] = 1. * self.args.params[process]['delta_t'] * MMTK.Units.fs
    params_k['steps_per_trial'] = self.args.params[process]['steps_per_sweep']

    attempts_left = 12
    while (attempts_left > 0):
      # Get initial potential energy
      Es_o = []
      for seed in seeds:
        self.top.universe.setConfiguration(
          Configuration(self.top.universe, seed))
        Es_o.append(self.top.universe.energy())
      Es_o = np.array(Es_o)

      # Perform simulation
      results = []
      if self.args.cores > 1:
        # Multiprocessing code
        import multiprocessing
        m = multiprocessing.Manager()
        task_queue = m.Queue()
        done_queue = m.Queue()
        for k in range(len(seeds)):
          task_queue.put((seeds[k], process, params_k, True, k))
        processes = [multiprocessing.Process(target=self.iterator.iteration_worker, \
            args=(task_queue, done_queue)) for p in range(self.args.cores)]
        for p in range(self.args.cores):
          task_queue.put('STOP')
        for p in processes:
          p.start()
        for p in processes:
          p.join()
        results = [done_queue.get() for seed in seeds]
        for p in processes:
          p.terminate()
      else:
        # Single process code
        results = [self.iterator.iteration(\
          seeds[k], process, params_k, True, k) for k in range(len(seeds))]

      seeds = [result['confs'] for result in results]
      Es_n = np.array([result['Etot'] for result in results])
      deltaEs = Es_n - Es_o
      attempts_left -= 1

      # Get the time step
      delta_t = np.array([result['delta_t'] for result in results])
      if np.std(delta_t) > 1E-3:
        # If the integrator adapts the time step, take an average
        delta_t = min(max(np.mean(delta_t), \
          self.args.params[process]['delta_t']/5.0*MMTK.Units.fs), \
          self.args.params[process]['delta_t']*0.1*MMTK.Units.fs)
      else:
        delta_t = delta_t[0]

      # Adjust the time step
      if self.args.params[process]['sampler'] == 'HMC':
        # Adjust the time step for Hamiltonian Monte Carlo
        acc_rate = float(np.sum([r['acc_Sampler'] for r in results]))/\
          np.sum([r['att_Sampler'] for r in results])
        if acc_rate > 0.8:
          delta_t += 0.125 * MMTK.Units.fs
        elif acc_rate < 0.4:
          if delta_t < 2.0 * MMTK.Units.fs:
            params_k['steps_per_trial'] = max(
              int(params_k['steps_per_trial'] / 2.), 1)
          delta_t -= 0.25 * MMTK.Units.fs
          if acc_rate < 0.1:
            delta_t -= 0.25 * MMTK.Units.fs
        else:
          attempts_left = 0
      else:
        # For other integrators, make sure the time step
        # is small enough to see changes in the energy
        if (np.std(deltaEs) < 1E-3):
          delta_t -= 0.25 * MMTK.Units.fs
        else:
          attempts_left = 0

      if delta_t < 0.1 * MMTK.Units.fs:
        delta_t = 0.1 * MMTK.Units.fs

      params_k['delta_t'] = delta_t

    sampler_metrics = ''
    for s in ['ExternalMC', 'SmartDarting', 'Sampler']:
      if np.array(['acc_' + s in r.keys() for r in results]).any():
        acc = np.sum([r['acc_' + s] for r in results])
        att = np.sum([r['att_' + s] for r in results])
        time = np.sum([r['time_' + s] for r in results])
        if att > 0:
          sampler_metrics += '%s %d/%d=%.2f (%.1f s); '%(\
            s,acc,att,float(acc)/att,time)
    return (seeds, Es_n - Es_o, delta_t, sampler_metrics)


  def _tL_tensor(self, E, params_c, process='CD'):
    # Metric tensor for the thermodynamic length
    T = params_c['T']
    deltaT = np.abs(self.args.params['BC']['T_HIGH'] - self.args.params['BC']['T_TARGET'])
    if process == 'CD':
      alpha = params_c['alpha']
      alpha_g = 4. * (alpha - 0.5)**2 / (1 + np.exp(-100 * (alpha - 0.5)))
      if alpha_g < 1E-10:
        alpha_g = 0
      da_g_da = (400.*(alpha-0.5)**2*np.exp(-100.*(alpha-0.5)))/(\
        1+np.exp(-100.*(alpha-0.5)))**2 + \
        (8.*(alpha-0.5))/(1 + np.exp(-100.*(alpha-0.5)))

      # Psi_g are terms that are scaled in with alpha_g
      # OBC is the strength of the OBC scaling in the current state
      s_ELE = 0.2 if self.args.params['CD']['solvation'] == 'Reduced' else 1.0
      if self.args.params['CD']['solvation'] == 'Desolvated':
        Psi_g = self._u_kln([E], [{'LJr': 1, 'LJa': 1, 'ELE': 1}], noBeta=True)
        OBC = 0.0
      elif self.args.params['CD']['solvation'] == 'Reduced':
        Psi_g = self._u_kln([E], [{
          'LJr': 1,
          'LJa': 1,
          'ELE': 0.2
        }],
                            noBeta=True)
        OBC = 0.0
      elif self.args.params['CD']['solvation'] == 'Fractional':
        Psi_g = self._u_kln([E], [{
          'LJr': 1,
          'LJa': 1,
          'ELE': 1,
          'OBC': 1
        }],
                            noBeta=True)
        OBC = alpha_g
      elif self.args.params['CD']['solvation'] == 'Full':
        Psi_g = self._u_kln([E], [{'LJr': 1, 'LJa': 1, 'ELE': 1}], noBeta=True)
        OBC = 1.0

      if self.args.params['CD']['pose'] > -1:  # Pose BPMF
        alpha_r = np.tanh(16 * alpha * alpha)
        da_r_da = 32. * alpha / np.cosh(16. * alpha * alpha)**2
        # Scaled in with alpha_r
        U_r = self._u_kln([E], \
          [{'k_angular_int':self.args.params['CD']['k_pose']}], noBeta=True)
        # Total potential energy
        U_RL_g = self._u_kln([E],
          [{'MM':True, 'OBC':OBC, 'T':T, \
            'k_angular_ext':params_c['k_angular_ext'], \
            'k_spatial_ext':params_c['k_spatial_ext'], \
            'k_angular_int':params_c['k_angular_int'], \
            'sLJr':alpha_g, 'sLJa':alpha_g, 'ELE':alpha_g}], noBeta=True)
        return np.abs(da_r_da)*U_r.std()/(R*T) + \
               np.abs(da_g_da)*Psi_g.std()/(R*T) + \
               deltaT*U_RL_g.std()/(R*T*T)
      else:  # BPMF
        alpha_sg = 1. - 4. * (alpha - 0.5)**2
        da_sg_da = -8 * (alpha - 0.5)
        if self.args.params['CD']['solvation'] == 'Reduced':
          Psi_sg = self._u_kln([E], [{'sLJr': 1}], noBeta=True)
          U_RL_g = self._u_kln([E],
            [{'MM':True, 'OBC':OBC, 'site':True, 'T':T,\
            'sLJr':alpha_sg, 'LJr':alpha_g, 'LJa':alpha_g, 'ELE':s_ELE*alpha_g}], noBeta=True)
        else:
          # Scaled in with soft grid
          Psi_sg = self._u_kln([E], [{'sLJr': 1, 'sELE': 1}], noBeta=True)
          # Total potential energy
          U_RL_g = self._u_kln([E],
            [{'MM':True, 'OBC':OBC, 'site':True, 'T':T,\
            'sLJr':alpha_sg, 'sELE':alpha_sg, \
            'LJr':alpha_g, 'LJa':alpha_g, 'ELE':alpha_g}], noBeta=True)
        return np.abs(da_sg_da)*Psi_sg.std()/(R*T) + \
               np.abs(da_g_da)*Psi_g.std()/(R*T) + \
               deltaT*U_RL_g.std()/(R*T*T)
    elif process == 'BC':
      if self.args.params['BC']['solvation'] == 'Full':
        # OBC is always on
        return deltaT * self._u_kln([E], [{
          'MM': True,
          'OBC': 1.0
        }],
                                    noBeta=True).std() / (R * T * T)
      else:
        # OBC is scaled with the progress variable
        alpha = params_c['alpha']
        return self._u_kln([E],[{'OBC':1.0}], noBeta=True).std()/(R*T) + \
          deltaT*self._u_kln([E],[{'MM':True, 'OBC':alpha}], noBeta=True).std()/(R*T*T)
    else:
      raise Exception("Unknown process!")

  def _initialize_BC(self, seeds):
    """Initializes the BC protocol

    Returns
    -------
    TODO
    """
    self.log.recordStart('BCsave')

    if seeds is not None:
      self.log.tee("\n>>> Initial warming, starting at " + \
        time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "\n")

      params_o = self.system.paramsFromAlpha(1.0, 'BC', site=False)
      self.data['BC'].protocol = [params_o]
      self.system.setParams(params_o)

      # Run at starting temperature
      self.log.recordStart('BC_state')
      (confs, DeltaEs, params_o['delta_t'], sampler_metrics) = \
        self._initial_sim_state(seeds, 'BC', params_o)
      E = self.system.energyTerms(confs, process='BC')
      self.data['BC'].confs['replicas'] = [
        confs[np.random.randint(len(confs))]
      ]
      self.data['BC'].confs['samples'] = [[confs]]
      self.data['BC'].Es = [[E]]

      self.log.tee("  at %d K in %s: %s" %
                   (params_o['T'], HMStime(
                     self.log.timeSince('BC_state')), sampler_metrics))
      self.log.tee("    dt=%.2f fs; tL_tensor=%.3e"%(\
        params_o['delta_t']*1000., self._tL_tensor(E,params_o,process='BC')))
    else:
      self.log.tee("\n>>> Initial warming, continuing at " + \
        time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()))
      params_o = self.data['BC'].protocol[-1]
      confs = self.data['BC'].confs['samples'][-1][0]
      E = self.data['BC'].Es[-1][0]

    if self.args.params['BC']['darts_per_seed'] > 0:
      self.data['BC'].confs['SmartDarting'] += confs

    # Main loop for initial BC:
    # choose new temperature, randomly select seeds, simulate
    rejectStage = 0
    while (not self.data['BC'].protocol[-1]['crossed']):
      # Choose new temperature
      params_n = self._next_BC_state(E = E, params_o = params_o, \
        pow = rejectStage)
      self.data['BC'].protocol.append(params_n)

      # Randomly select seeds for new trajectory
      u_o = self._u_kln([E], [params_o])
      u_n = self._u_kln([E], [params_n])
      du = u_n - u_o
      weights = np.exp(-du + min(du))
      seedIndicies = np.random.choice(len(u_o), \
        size = self.args.params['BC']['seeds_per_state'], \
        p = weights/sum(weights))
      seeds = [np.copy(confs[s]) for s in seedIndicies]

      # Store old data
      confs_o = confs
      E_o = E

      # Simulate
      self.log.recordStart('BC_state')
      self.system.setParams(params_n)
      self.log.tee(
        self.iterator.initializeSmartDartingConfigurations(
          seeds, 'BC', self.data))
      (confs, DeltaEs, params_n['delta_t'], sampler_metrics) = \
        self._initial_sim_state(seeds, 'BC', params_n)

      if self.args.params['BC']['darts_per_seed'] > 0:
        self.data['BC'].confs['SmartDarting'] += confs

      # Get state energies
      E = self.system.energyTerms(confs, process='BC')

      # Estimate the mean replica exchange acceptance rate
      # between the previous and new state
      (u_kln, N_k) = self._u_kln([[E_o], [E]], self.data['BC'].protocol[-2:])
      N = min(N_k)
      acc = np.exp(-u_kln[0, 1, :N] - u_kln[1, 0, :N] + u_kln[0, 0, :N] +
                   u_kln[1, 1, :N])
      mean_acc = np.mean(np.minimum(acc, np.ones(acc.shape)))

      self.log.tee("  at %d K in %s: %s"%(\
        params_n['T'], HMStime(self.log.timeSince('BC_state')), sampler_metrics))
      self.log.tee("    dt=%.2f fs; tL_tensor=%.2e; <acc>=%.2f"%(\
        params_n['delta_t']*1000., \
        self._tL_tensor(E,params_o,process='BC'), mean_acc))

      if (mean_acc<self.args.params['BC']['min_repX_acc']) and \
          (self.args.params['BC']['protocol']=='Adaptive'):
        # If the acceptance probability is too low,
        # reject the state and restart
        self.data['BC'].protocol.pop()
        confs = confs_o
        E = E_o
        rejectStage += 1
        self.log.tee("  rejected new state with low acceptance rate")
      elif (len(self.data['BC'].protocol)>2) and \
          (mean_acc>0.99) and (not params_n['crossed']) and \
          (self.args.params['BC']['protocol'] == 'Adaptive'):
        # If the acceptance probability is too high,
        # reject the previous state and restart
        self.data['BC'].confs['replicas'][-1] = confs[np.random.randint(
          len(confs))]
        self.data['BC'].protocol.pop()
        self.data['BC'].protocol[-1] = copy.deepcopy(params_n)
        rejectStage -= 1
        params_o = params_n
        self.log.tee("  rejected previous state with high acceptance rate")
      else:
        # Store data and continue with initialization
        self.data['BC'].confs['replicas'].append(confs[np.random.randint(
          len(confs))])
        self.data['BC'].confs['samples'].append([confs])
        self.data['BC'].Es.append([E])
        self.data['BC'].protocol[-1] = copy.deepcopy(params_n)
        rejectStage = 0
        params_o = params_n

      # Special tasks after the last stage
      if self.data['BC'].protocol[-1]['crossed']:
        # For warming, reverse protocol and energies
        self.log.tee("  reversing replicas, samples, and protocol")
        self.data['BC'].confs['replicas'].reverse()
        self.data['BC'].confs['samples'].reverse()
        self.data['BC'].Es.reverse()
        self.data['BC'].protocol.reverse()
        self.data['BC'].protocol[0]['crossed'] = False
        self.data['BC'].protocol[-1]['crossed'] = True

        if (not self.args.params['BC']['keep_intermediate']):
          for k in range(1, len(self.data['BC'].protocol) - 1):
            self.data['BC'].confs['samples'][k] = []

        self.data['BC'].cycle += 1

      # Save progress every 5 minutes
      if (self.log.timeSince('BCsave') > 5 * 60):
        self.save('BC')
        self.log.recordStart('BCsave')
        saved = True
      else:
        saved = False

      if not self.log.isTimeRemaining():
        if not saved:
          self.save('BC')
          self.log.tee("")
        self.log.tee("  no time remaining for initial BC")
        self.log.clear_lock('BC')
        return False

    # Save data
    if not saved:
      self.save('BC')
      self.log.tee("")

    self.log.tee("Elapsed time for initial warming of " + \
      "%d states: "%len(self.data['BC'].protocol) + \
      HMStime(self.log.timeSince('BC')))
    self.log.clear_lock('BC')
    self.iterator.clearSmartDartingConfigurations('BC')
    return True

  def _next_BC_state(self, E=None, params_o=None, pow=None, warm=True):
    if E is None:
      E = self.data['BC'].Es[-1]

    if params_o is None:
      params_o = self.data['BC'].protocol[-1]

    if self.args.params['BC']['protocol'] == 'Adaptive':
      tL_tensor = self._tL_tensor(E, params_o, process='BC')
      crossed = params_o['crossed']
      if pow is not None:
        tL_tensor = tL_tensor * (1.25**pow)
      if tL_tensor > 1E-7:
        dL = self.args.params['BC']['therm_speed'] / tL_tensor
        if warm:
          T = params_o['T'] + dL
          if T > self.args.params['BC']['T_HIGH']:
            T = self.args.params['BC']['T_HIGH']
            crossed = True
        else:
          T = params_o['T'] - dL
          if T < self.args.params['BC']['T_TARGET']:
            T = self.args.params['BC']['T_TARGET']
            crossed = True
      else:
        raise Exception('No variance in configuration energies')
    elif self.args.params['BC']['protocol'] == 'Geometric':
      T_GEOMETRIC = np.exp(
        np.linspace(np.log(self.args.params['BC']['T_TARGET']), np.log(self.args.params['BC']['T_HIGH']),
                    int(1 / self.args.params['BC']['therm_speed'])))
      if warm:
        T_START, T_END = self.args.params['BC']['T_TARGET'], self.args.params['BC']['T_HIGH']
      else:
        T_START, T_END = self.args.params['BC']['T_HIGH'], self.args.params['BC']['T_TARGET']
        T_GEOMETRIC = T_GEOMETRIC[::-1]
      T = T_GEOMETRIC[len(self.data['BC'].protocol)]
      crossed = (len(self.data['BC'].protocol) == (len(T_GEOMETRIC) - 1))
    alpha = (self.args.params['BC']['T_HIGH'] - T) / (self.args.params['BC']['T_HIGH'] - self.args.params['BC']['T_TARGET'])
    return self.system.paramsFromAlpha(alpha,
                                        process='BC',
                                        params_o=params_o,
                                        crossed=crossed)

  def _initialize_CD(self, seeds):
    self.log.recordStart('initial_CD')
    self.log.recordStart('CDsave')

    if self.data['CD'].protocol == []:
      if seeds is not None:
        decoupling = True

        params_o = self.system.paramsFromAlpha(1.0, 'CD')
        self.data['CD'].protocol = [params_o]
        self.system.setParams(params_o)

        attempts = 0
        DeltaEs = np.array([0.])
        # Iteratively simulate and modify the time step
        # until the simulations lead to different energies
        while np.std(DeltaEs) < 1E-7:
          attempts += 1
          if attempts == 5:
            self._store_infinite_f_RL()
            raise Exception('Unable to ramp temperature')

          (confs, DeltaEs, params_o['delta_t'], sampler_metrics) = \
            self._initial_sim_state(seeds, 'CD', params_o)

        # Get state energies
        E = self.system.energyTerms(confs)
        self.data['CD'].confs['replicas'] = [
          confs[np.random.randint(len(confs))]
        ]
        self.data['CD'].confs['samples'] = [[confs]]
        self.data['CD'].Es = [[E]]

        self.log.tee("\n  at a=%.3e in %s: %s"%(\
          params_o['alpha'], HMStime(self.log.timeSince('initial_CD')), sampler_metrics))
        self.log.tee("    dt=%.2f fs, tL_tensor=%.3e"%(\
          params_o['delta_t']*1000., \
          self._tL_tensor(E,params_o)))
      else:
        decoupling = False

        self.log.tee("\n>>> Initial CD, starting at " + \
          time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "\n")

        # Select samples from the high T unbound state and ensure there are enough
        confs_HT = []
        for k in range(1, len(self.data['BC'].Es[0])):
          confs_HT += list(self.data['BC'].confs['samples'][0][k])
        while len(confs_HT) < self.args.params['CD']['seeds_per_state']:
          self.log.tee(
            "More samples from high temperature ligand simulation needed")
          self.log.clear_lock('CD')
          self._replica_exchange('BC')
          self.log.set_lock('CD')
          confs_HT = []
          for k in range(1, len(self.data['BC'].Es[0])):
            confs_HT += list(self.data['BC'].confs['samples'][0][k])
        confs_HT = confs_HT[:self.args.params['CD']['seeds_per_state']]

        if (self.args.params['CD']['pose'] == -1):
          (confs, E) = self._random_CD()
          self.log.tee("  random CD complete in " +
                       HMStime(self.log.timeSince('initial_CD')))
          if randomOnly:
            self.log.clear_lock('CD')
            return
        else:
          params_o = self.system.paramsFromAlpha(0, 'CD')
          self.system.setParams(params_o)
          confs = confs_HT
          E = self.system.energyTerms(confs)
          self.data['CD'].confs['replicas'] = [
            confs[np.random.randint(len(confs))]
          ]
          self.data['CD'].confs['samples'] = [confs]
          self.data['CD'].Es = [[E]]
          self.data['CD'].protocol = [params_o]
    else:
      # Continuing from a previous CD instance
      decoupling = self.data['CD'].protocol[0]['alpha'] > self.data['CD'].protocol[1]['alpha']
      self.log.tee("\n>>> Initial %sCD, "%({True:'un',False:''}[decoupling]) + \
        "continuing at " + \
        time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()))
      confs = self.data['CD'].confs['samples'][-1][0]
      E = self.data['CD'].Es[-1][0]

    if self.args.params['CD']['darts_per_seed'] > 0:
      self.data['CD'].confs['SmartDarting'] += confs

    params_o = self.data['CD'].protocol[-1]

    # Main loop for initial CD:
    # choose new thermodynamic variables,
    # randomly select seeds,
    # simulate
    rejectStage = 0
    while (not self.data['CD'].protocol[-1]['crossed']):
      # Determine next value of the protocol
      params_n = self._next_CD_state(E = E, params_o = params_o, \
          pow = rejectStage, decoupling = decoupling)
      self.data['CD'].protocol.append(params_n)
      if len(self.data['CD'].protocol) > 1000:
        self._clear('CD')
        self.save('CD')
        self._store_infinite_f_RL()
        raise Exception('Too many replicas!')
      if abs(rejectStage) > 20:
        self._clear('CD')
        self.save('CD')
        self._store_infinite_f_RL()
        raise Exception('Too many consecutive rejected stages!')

      # Randomly select seeds for new trajectory
      u_o = self._u_kln([E], [params_o])
      u_n = self._u_kln([E], [params_n])
      du = u_n - u_o
      weights = np.exp(-du + min(du))
      seedIndicies = np.random.choice(len(u_o), \
        size = self.args.params['CD']['seeds_per_state'], \
        p=weights/sum(weights))

      if (not decoupling) and (self.args.params['CD']['pose'] == -1) \
          and (len(self.data['CD'].protocol)==2):
        # BC state 0 configurations, randomly oriented
        # Use the lowest energy configuration
        # in the first CD state for replica exchange
        ind = np.argmin(u_n)
        (c,i_rot,i_trans) = np.unravel_index(ind, \
          (self.args.params['CD']['seeds_per_state'], self._n_rot, self.data['CD']._n_trans))
        repX_conf = np.add(np.dot(confs[c], self._random_rotT[i_rot,:,:]),\
                           self.data['CD']._random_trans[i_trans].array)
        self.data['CD'].confs['replicas'] = [repX_conf]
        self.data['CD'].confs['samples'] = [[repX_conf]]
        self.data['CD'].Es = [[dict([(key,np.array([val[ind]])) \
          for (key,val) in E.iteritems()])]]
        seeds = []
        for ind in seedIndicies:
          (c,i_rot,i_trans) = np.unravel_index(ind, \
            (self.args.params['CD']['seeds_per_state'], self._n_rot, self._n_trans))
          seeds.append(np.add(np.dot(confs[c], self._random_rotT[i_rot,:,:]), \
            self.data['CD']._random_trans[i_trans].array))
        confs = None
        E = {}
      else:  # Seeds from last state
        seeds = [np.copy(confs[ind]) for ind in seedIndicies]
      self.data['CD'].confs['seeds'] = seeds

      # Store old data
      confs_o = confs
      E_o = E

      # Simulate
      self.log.recordStart('CD_state')
      self.system.setParams(params_n)
      if self.args.params['CD']['darts_per_seed'] > 0 and params_n['alpha'] > 0.1:
        self.log.tee(
          self.iterator.initializeSmartDartingConfigurations(
            self.data['CD'].confs['SmartDarting'], 'CD', self.data))
      (confs, DeltaEs, params_n['delta_t'], sampler_metrics) = \
        self._initial_sim_state(seeds, 'CD', params_n)
      if np.std(DeltaEs) < 1E-7:
        self._store_infinite_f_RL()
        raise Exception('Unable to initialize simulation')

      if self.args.params['CD']['darts_per_seed'] > 0:
        self.data['CD'].confs['SmartDarting'] += confs

      # Get state energies
      E = self.system.energyTerms(confs)

      # Estimate the mean replica exchange acceptance rate
      # between the previous and new state
      self.log.tee("  at a=%.3e in %s: %s"%(\
        params_n['alpha'], \
        HMStime(self.log.timeSince('CD_state')), sampler_metrics))

      if E_o != {}:
        (u_kln, N_k) = self._u_kln([[E_o], [E]], self.data['CD'].protocol[-2:])
        N = min(N_k)
        acc = np.exp(-u_kln[0, 1, :N] - u_kln[1, 0, :N] + u_kln[0, 0, :N] +
                     u_kln[1, 1, :N])
        mean_acc = np.mean(np.minimum(acc, np.ones(acc.shape)))

        self.log.tee("    dt=%.2f fs; tL_tensor=%.3e; <acc>=%.2f"%(\
          params_n['delta_t']*1000., self._tL_tensor(E,params_o), mean_acc))
      else:
        self.log.tee("    dt=%.2f fs; tL_tensor=%.3e"%(\
          params_n['delta_t']*1000., self._tL_tensor(E,params_o)))

      # Decide whether to keep the state
      if len(self.data['CD'].protocol) > (1 + (not decoupling)):
        if (mean_acc<self.args.params['CD']['min_repX_acc']) and \
            (self.args.params['CD']['protocol']=='Adaptive'):
          # If the acceptance probability is too low,
          # reject the state and restart
          self.data['CD'].protocol.pop()
          confs = confs_o
          E = E_o
          rejectStage += 1
          self.log.tee(
            "  rejected new state with low estimated acceptance rate")
        elif len(self.data['CD'].protocol)>(2+(not decoupling)) and \
            (mean_acc>0.99) and (not params_n['crossed']) and \
            (self.args.params['CD']['protocol']=='Adaptive'):
          # If the acceptance probability is too high,
          # reject the previous state and restart
          self.data['CD'].confs['replicas'][-1] = confs[np.random.randint(
            len(confs))]
          self.data['CD'].protocol.pop()
          self.data['CD'].protocol[-1] = copy.deepcopy(params_n)
          rejectStage -= 1
          params_o = params_n
          self.log.tee(
            "  rejected past state with high estimated acceptance rate")
        else:
          # Store data and continue with initialization
          self.data['CD'].confs['replicas'].append(confs[np.random.randint(
            len(confs))])
          self.data['CD'].confs['samples'].append([confs])
          self.data['CD'].Es.append([E])
          self.data['CD'].protocol[-1] = copy.deepcopy(params_n)
          rejectStage = 0
          params_o = params_n
      else:
        # Store data and continue with initialization (first time)
        self.data['CD'].confs['replicas'].append(confs[np.random.randint(
          len(confs))])
        self.data['CD'].confs['samples'].append([confs])
        self.data['CD'].Es.append([E])
        self.data['CD'].protocol[-1] = copy.deepcopy(params_n)
        rejectStage = 0
        params_o = params_n

      # Special tasks after the last stage
      if (self.data['CD'].protocol[-1]['crossed']):
        # For decoupling, reverse protocol and energies
        if decoupling:
          self.log.tee("  reversing replicas, samples, and protocol")
          self.data['CD'].confs['replicas'].reverse()
          self.data['CD'].confs['samples'].reverse()
          self.data['CD'].confs['seeds'] = None
          self.data['CD'].Es.reverse()
          self.data['CD'].protocol.reverse()
          self.data['CD'].protocol[0]['crossed'] = False
          self.data['CD'].protocol[-1]['crossed'] = True

        if (not self.args.params['CD']['keep_intermediate']):
          for k in range(1, len(self.data['CD'].protocol) - 1):
            self.data['CD'].confs['samples'][k] = []

        self.data['CD'].cycle += 1

      # Save progress every 10 minutes
      if (self.log.timeSince('CDsave') > (10 * 60)):
        self.save('CD')
        self.log.recordStart('CDsave')
        saved = True
      else:
        saved = False

      if not self.log.isTimeRemaining():
        if not saved:
          self.save('CD')
        self.log.tee("  no time remaining for initial CD")
        self.log.clear_lock('CD')
        return False

    if not saved:
      self.save('CD')

    self.log.tee("\nElapsed time for initial CD of " + \
      "%d states: "%len(self.data['CD'].protocol) + \
      HMStime(self.log.timeSince('initial_CD')))
    self.log.clear_lock('CD')
    self.iterator.clearSmartDartingConfigurations('CD')
    return True

  def _random_CD(self):
    """
      Randomly places the ligand into the receptor and evaluates energies

      The first state of CD is sampled by randomly placing configurations
      from the high temperature ligand simulation into the binding site.
    """
    # Select samples from the high T unbound state
    E_MM = []
    E_OBC = []
    confs = []
    for k in range(1, len(self.data['BC'].Es[0])):
      E_MM += list(self.data['BC'].Es[0][k]['MM'])
      if ('OBC' in self.data['BC'].Es[0][k].keys()):
        E_OBC += list(self.data['BC'].Es[0][k]['OBC'])
      confs += list(self.data['BC'].confs['samples'][0][k])

    random_CD_inds = np.array(np.linspace(0,len(E_MM), \
      self.args.params['CD']['seeds_per_state'],endpoint=False),dtype=int)
    BC0_Es_MM = [E_MM[ind] for ind in random_CD_inds]
    BC0_Es_OBC = []
    if E_OBC != []:
      BC0_Es_OBC = [E_OBC[ind] for ind in random_CD_inds]
    BC0_confs = [confs[ind] for ind in random_CD_inds]

    # Do the random CD
    params_o = self.system.paramsFromAlpha(0.0, 'CD')
    params_o['delta_t'] = 1. * self.data['BC'].protocol[0]['delta_t']
    self.data['CD'].protocol = [params_o]

    # Set up the force field with full interaction grids
    self.system.setParams(self.system.paramsFromAlpha(1.0, 'CD'))

    # Either loads or generates the random translations and rotations for the first state of CD
    if not (hasattr(self, '_random_trans') and hasattr(self, '_random_rotT')):
      self.data['CD']._max_n_trans = 10000
      # Default density of points is 50 per nm**3
      self.data['CD']._n_trans = max(
        min(
          np.int(
            np.ceil(self._forceFields['site'].volume *
                    self.args.params['CD']['site_density'])),
          self.data['CD']._max_n_trans), 5)
      self.data['CD']._random_trans = np.ndarray(
        (self.data['CD']._max_n_trans), dtype=Vector)
      for ind in range(self.data['CD']._max_n_trans):
        self.data['CD']._random_trans[ind] = Vector(
          self._forceFields['site'].randomPoint())
      self.data['CD']._max_n_rot = 100
      self._n_rot = 100
      self._random_rotT = np.ndarray((self.data['CD']._max_n_rot, 3, 3))
      from AlGDock.Integrators.ExternalMC.ExternalMC import random_rotate
      for ind in range(self.data['CD']._max_n_rot):
        self._random_rotT[ind, :, :] = np.transpose(random_rotate())
    else:
      self.data['CD']._max_n_trans = self.data['CD']._random_trans.shape[0]
      self._n_rot = self._random_rotT.shape[0]

    # Get interaction energies.
    # Loop over configurations, random rotations, and random translations
    E = {}
    for term in (['MM', 'site'] + scalables):
      # Large array creation may cause MemoryError
      E[term] = np.zeros((self.args.params['CD']['seeds_per_state'], \
        self.data['CD']._max_n_rot,self.data['CD']._n_trans))
    self.log.tee("  allocated memory for interaction energies")

    from AlGDock.system import term_map

    converged = False
    n_trans_o = 0
    n_trans_n = self.data['CD']._n_trans
    while not converged:
      for c in range(self.args.params['CD']['seeds_per_state']):
        E['MM'][c, :, :] = BC0_Es_MM[c]
        if BC0_Es_OBC != []:
          E['OBC'][c, :, :] = BC0_Es_OBC[c]
        for i_rot in range(self._n_rot):
          conf_rot = Configuration(self.top.universe,\
            np.dot(BC0_confs[c], self._random_rotT[i_rot,:,:]))
          for i_trans in range(n_trans_o, n_trans_n):
            self.top.universe.setConfiguration(conf_rot)
            self.top.universe.translateTo(
              self.data['CD']._random_trans[i_trans])
            eT = self.top.universe.energyTerms()
            for (key, value) in eT.iteritems():
              if key != 'electrostatic':  # For some reason, MMTK double-counts electrostatic energies
                E[term_map[key]][c, i_rot, i_trans] += value
      E_c = {}
      for term in E.keys():
        # Large array creation may cause MemoryError
        E_c[term] = np.ravel(E[term][:, :self._n_rot, :n_trans_n])
      self.log.tee("  allocated memory for %d translations" % n_trans_n)
      (u_kln,N_k) = self._u_kln([E_c],\
        [params_o,self._next_CD_state(E=E_c, params_o=params_o, decoupling=False)])
      du = u_kln[0, 1, :] - u_kln[0, 0, :]
      bootstrap_reps = 50
      f_grid0 = np.zeros(bootstrap_reps)
      for b in range(bootstrap_reps):
        du_b = du[np.random.randint(0, len(du), len(du))]
        f_grid0[b] = -np.log(np.exp(-du_b + min(du_b)).mean()) + min(du_b)
      f_grid0_std = f_grid0.std()
      converged = f_grid0_std < 0.1
      if not converged:
        self.log.tee("  with %s translations "%n_trans_n + \
                 "the predicted free energy difference is %.5g (%.5g)"%(\
                 f_grid0.mean(),f_grid0_std))
        if n_trans_n == self.data['CD']._max_n_trans:
          break
        n_trans_o = n_trans_n
        n_trans_n = min(n_trans_n + 25, self.data['CD']._max_n_trans)
        for term in (['MM', 'site'] + scalables):
          # Large array creation may cause MemoryError
          E[term] = np.dstack((E[term], \
            np.zeros((self.args.params['CD']['seeds_per_state'],\
              self.data['CD']._max_n_rot,25))))

    if self.data['CD']._n_trans != n_trans_n:
      self.data['CD']._n_trans = n_trans_n

    self.log.tee("  %d ligand configurations "%len(BC0_Es_MM) + \
             "were randomly docked into the binding site using "+ \
             "%d translations and %d rotations "%(n_trans_n,self._n_rot))
    self.log.tee("  the predicted free energy difference between the" + \
             " first and second CD states is " + \
             "%.5g (%.5g)"%(f_grid0.mean(),f_grid0_std))

    self.log.recordStart('ravel')
    for term in E.keys():
      E[term] = np.ravel(E[term][:, :self._n_rot, :self.data['CD']._n_trans])
    self.log.tee("  raveled energy terms in " + \
      HMStime(self.log.timeSince('ravel')))

    return (BC0_confs, E)

  def _next_CD_state(self, E=None, params_o=None, pow=None, decoupling=False):
    """
    Determines the parameters for the next CD state
    """

    if E is None:
      E = self.data['CD'].Es[-1]

    if params_o is None:
      params_o = self.data['CD'].protocol[-1]

    if self.args.params['CD']['protocol'] == 'Adaptive':
      # Change grid scaling and temperature simultaneously
      tL_tensor = self._tL_tensor(E, params_o)
      crossed = params_o['crossed']
      # Calculate the change in the progress variable, capping at 0.05
      if pow is None:
        dL = min(self.args.params['CD']['therm_speed'] / tL_tensor, 0.05)
      else:
        # If there have been rejected stages, reduce dL
        dL = min(
          self.args.params['CD']['therm_speed'] / tL_tensor / (1.25**pow),
          0.05)
      if decoupling:
        alpha = params_o['alpha'] - dL
        if (self.args.params['CD']['pose'] > -1) and \
           (params_o['alpha'] > 0.5) and (alpha < 0.5):
          # Stop at 0.5 to facilitate entropy-energy decomposition
          alpha = 0.5
        elif alpha < 0.0:
          if pow > 0:
            alpha = params_o['alpha'] * (1 - 0.8**pow)
          else:
            alpha = 0.0
            crossed = True
      else:
        alpha = params_o['alpha'] + dL
        if (self.args.params['CD']['pose'] > -1) and \
           (params_o['alpha'] < 0.5) and (alpha > 0.5):
          # Stop at 0.5 to facilitate entropy-energy decomposition
          alpha = 0.5
        elif alpha > 1.0:
          if pow > 0:
            alpha = params_o['alpha'] + (1 - params_o['alpha']) * 0.8**pow
          else:
            alpha = 1.0
            crossed = True
      return self.system.paramsFromAlpha(alpha,
                                          process='CD',
                                          params_o=params_o,
                                          crossed=crossed)
    elif self.args.params['CD']['protocol'] == 'Geometric':
      A_GEOMETRIC = [0.] + list(
        np.exp(
          np.linspace(np.log(1E-10), np.log(1.0),
                      int(1 / self.args.params['CD']['therm_speed']))))
      if decoupling:
        A_GEOMETRIC.reverse()
      alpha = A_GEOMETRIC[len(self.data['CD'].protocol)]
      crossed = len(self.data['CD'].protocol) == (len(alpha_GEOMETRIC) - 1)
      return self.system.paramsFromAlpha(alpha,
                                          process='CD',
                                          params_o=params_o,
                                          crossed=crossed)
