# The example is 1of6

import os, shutil, glob
os.system('rm -rf dock')

import AlGDock.BindingPMF_plots
from AlGDock.BindingPMF_plots import *

self = AlGDock.BindingPMF_plots.BPMF_plots(\
  dir_dock='dock', dir_cool='cool',\
  ligand_tarball='prmtopcrd/ligand.tar.gz', \
  ligand_database='ligand.db', \
  forcefield='prmtopcrd/gaff.dat', \
  ligand_prmtop='ligand.prmtop', \
  ligand_inpcrd='ligand.trans.inpcrd', \
  receptor_tarball='prmtopcrd/receptor.tar.gz', \
  receptor_prmtop='receptor.prmtop', \
  receptor_inpcrd='receptor.trans.inpcrd', \
  receptor_fixed_atoms='receptor.pdb', \
  complex_tarball='prmtopcrd/complex.tar.gz', \
  complex_prmtop='complex.prmtop', \
  complex_inpcrd='complex.trans.inpcrd', \
  complex_fixed_atoms='complex.pdb', \
  score = 'prmtopcrd/anchor_and_grow_scored.mol2', \
  dir_grid='grids', \
  protocol='Adaptive', cool_therm_speed=1.5, dock_therm_speed=1.5,\
  sampler='NUTS', \
  MCMC_moves=1, \
  seeds_per_state=10, steps_per_seed=200, darts_per_seed=5, \
  sweeps_per_cycle=25, attempts_per_sweep=100, \
  steps_per_sweep=50, darts_per_sweep=1, \
  cool_repX_cycles=3, dock_repX_cycles=4, \
  site='Sphere', site_center=[1.74395, 1.74395, 1.74395], site_max_R=0.6, \
  site_density=10., \
  phases=['NAMD_Gas','NAMD_GBSA'], \
  cores=-1, \
  rmsd=True,
  random_seed=50)

randomOnly=False
undock=True

dock_start_time = time.time()

if self.dock_protocol==[]:
  self.tee("\n>>> Initial docking, starting at " + \
    time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime()))
  if undock:
    lambda_o = self._lambda(1.0, 'dock', MM=True, site=True, crossed=False)
    self.dock_protocol = [lambda_o]
    self._set_universe_evaluator(lambda_o)
    seeds = self._get_confs_to_rescore(site=True, minimize=True)[0]

    if seeds==[]:
      undock = False
    else:
      self.confs['dock']['starting_poses'] = seeds
      # initializes smart darting for docking and sets the universe
      # to the lowest energy configuration
      if self.params['dock']['darts_per_seed']>0:
        self.tee(self.sampler['dock_SmartDarting'].set_confs(seeds))
      elif len(seeds)>0:
        self.universe.setConfiguration(Configuration(self.universe,seeds[-1]))
      
      # Ramp up the temperature
      T_LOW = 20.
      T_SERIES = T_LOW*(self.T_TARGET/T_LOW)**(np.arange(30)/29.)
      for T in T_SERIES:
        random_seed = int(abs(seeds[0][0][0]*10000)) + int(T*10000)
        if self._random_seed==0:
          random_seed += int(time.time())
        self.sampler['dock'](steps = 500, steps_per_trial = 100, T=T,\
                             delta_t=self.delta_t, random_seed=random_seed)
      seeds = [self.universe.configuration().array]

      # Simulate
      sim_start_time = time.time()
      (confs, Es_tot, lambda_o['delta_t'], sampler_metrics) = \
        self._initial_sim_state(\
          seeds*self.params['dock']['seeds_per_state'], 'dock', lambda_o)

      # Get state energies
      E = self._calc_E(confs)
      self.confs['dock']['replicas'] = [confs[np.random.randint(len(confs))]]
      self.confs['dock']['samples'] = [[confs]]
      self.dock_Es = [[E]]

      self.tee("  generated %d configurations "%len(confs) + \
               "with progress %e "%lambda_o['a'] + \
               "in " + HMStime(time.time()-sim_start_time))
      self.tee(sampler_metrics)
      self.tee("  dt=%.3f ps, tL_tensor=%.3e"%(\
        lambda_o['delta_t']*1000., \
        self._tL_tensor(E,lambda_o)))

  if not undock:
    (cool0_confs, E) = self.random_dock()
    self.tee("  random docking complete in " + \
             HMStime(time.time()-dock_start_time))
    if randomOnly:
      self._clear_lock('dock')
else:
  # Continuing from a previous docking instance
  self.tee("\n>>> Initial docking, continuing at " + \
    time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime()))
  confs = self.confs['dock']['samples'][-1][0]
  E = self.dock_Es[-1][0]

lambda_o = self.dock_protocol[-1]

# Main loop for initial docking:
# choose new thermodynamic variables,
# randomly select seeds,
# simulate
rejectStage = 0
# while (not self.dock_protocol[-1]['crossed']):
for protocol_index in range(3):
  # Determine next value of the protocol
  lambda_n = self._next_dock_state(E = E, lambda_o = lambda_o, \
      pow = rejectStage, undock = undock)
  self.dock_protocol.append(lambda_n)
  if len(self.dock_protocol)>1000:
    self._clear('dock')
    self._save('dock')
    self._store_infinite_f_RL()
    raise Exception('Too many replicas!')
  if abs(rejectStage)>20:
    self._clear('dock')
    self._save('dock')
    self._store_infinite_f_RL()
    raise Exception('Too many consecutive rejected stages!')

  # Randomly select seeds for new trajectory
  u_o = self._u_kln([E],[lambda_o])
  u_n = self._u_kln([E],[lambda_n])
  du = u_n-u_o
  weights = np.exp(-du+min(du))
  seedIndicies = np.random.choice(len(u_o), \
    size = self.params['dock']['seeds_per_state'], \
    p=weights/sum(weights))

  if (not undock) and (len(self.dock_protocol)==2):
    # Cooling state 0 configurations, randomly oriented
    # Use the lowest energy configuration in the first docking state for replica exchange
    ind = np.argmin(u_n)
    (c,i_rot,i_trans) = np.unravel_index(ind, (self.params['dock']['seeds_per_state'], self._n_rot, self._n_trans))
    repX_conf = np.add(np.dot(cool0_confs[c], self._random_rotT[i_rot,:,:]),\
                       self._random_trans[i_trans].array)
    self.confs['dock']['replicas'] = [repX_conf]
    self.confs['dock']['samples'] = [[repX_conf]]
    self.dock_Es = [[dict([(key,np.array([val[ind]])) for (key,val) in E.iteritems()])]]
    seeds = []
    for ind in seedIndicies:
      (c,i_rot,i_trans) = np.unravel_index(ind, (self.params['dock']['seeds_per_state'], self._n_rot, self._n_trans))
      seeds.append(np.add(np.dot(cool0_confs[c], self._random_rotT[i_rot,:,:]), self._random_trans[i_trans].array))
    confs = None
    E = {}
  else: # Seeds from last state
    seeds = [np.copy(confs[ind]) for ind in seedIndicies]
  self.confs['dock']['seeds'] = seeds

  # Store old data
  confs_o = confs
  E_o = E

  # Simulate
  sim_start_time = time.time()
  self._set_universe_evaluator(lambda_n)
  if self.params['dock']['darts_per_seed']>0:
    self.tee(self.sampler['dock_SmartDarting'].set_confs(confs, append=True))
  (confs, Es_tot, lambda_n['delta_t'], sampler_metrics) = \
    self._initial_sim_state(seeds, 'dock', lambda_n)

  # Get state energies
  E = self._calc_E(confs)

  self.tee("  generated %d configurations "%len(confs) + \
           "with progress %f "%lambda_n['a'] + \
           "in " + HMStime(time.time()-sim_start_time))
  self.tee(sampler_metrics)
  self.tee("  dt=%.3f ps, tL_tensor=%.3e"%(\
    lambda_n['delta_t']*1000.,
    self._tL_tensor(E,lambda_n)))

  # Decide whether to keep the state
  if len(self.dock_protocol)>(1+(not undock)):
    # Estimate the mean replica exchange acceptance rate
    # between the previous and new state
    (u_kln,N_k) = self._u_kln([[E_o],[E]], self.dock_protocol[-2:])
    N = min(N_k)
    acc = np.exp(-u_kln[0,1,:N]-u_kln[1,0,:N]+u_kln[0,0,:N]+u_kln[1,1,:N])
    mean_acc = np.mean(np.minimum(acc,np.ones(acc.shape)))
    
    if (mean_acc<self.params['dock']['min_repX_acc']):
      # If the acceptance probability is too low,
      # reject the state and restart
      self.dock_protocol.pop()
      confs = confs_o
      E = E_o
      rejectStage += 1
      self.tee("  rejected new state, as estimated replica exchange acceptance rate of %f is too low"%mean_acc)
    elif (mean_acc>0.99) and (not lambda_n['crossed']):
      # If the acceptance probability is too high,
      # reject the previous state and restart
      self.confs['dock']['replicas'][-1] = confs[np.random.randint(len(confs))]
      self.dock_protocol.pop()
      self.dock_protocol[-1] = copy.deepcopy(lambda_n)
      rejectStage -= 1
      lambda_o = lambda_n
      self.tee("  rejected previous state, as estimated replica exchange acceptance rate of %f is too high"%mean_acc)
    else:
      # Store data and continue with initialization
      self.confs['dock']['replicas'].append(confs[np.random.randint(len(confs))])
      self.confs['dock']['samples'].append([confs])
      self.dock_Es.append([E])
      self.dock_protocol[-1] = copy.deepcopy(lambda_n)
      rejectStage = 0
      lambda_o = lambda_n
      self.tee("  the estimated replica exchange acceptance rate is %f\n"%mean_acc)

      if (not self.params['dock']['keep_intermediate']):
        if len(self.dock_protocol)>(2+(not undock)):
          self.confs['dock']['samples'][-2] = []
  else:
    # Store data and continue with initialization (first time)
    self.confs['dock']['replicas'].append(confs[np.random.randint(len(confs))])
    self.confs['dock']['samples'].append([confs])
    self.dock_Es.append([E])
    self.dock_protocol[-1] = copy.deepcopy(lambda_n)
    rejectStage = 0
    lambda_o = lambda_n

#
from MMTK import Units
RT = 8.3144621*Units.J/Units.mol/Units.K*self.T

BPMF = self
self = BPMF.sampler['dock_SmartDarting']

self.universe.setConfiguration(Configuration(self.universe, self.confs[0]))

#
ntrials = 10

acc = 0.
energies = []
closest_poses = []

xo_Cartesian = np.copy(self.universe.configuration().array)
xo_BAT = self._BAT_util.BAT(xo_Cartesian, extended=self.extended)
eo = self.universe.energy()
if self.extended:
  closest_pose_o = self._closest_pose_Cartesian(\
    xo_Cartesian[self.molecule.heavy_atoms,:])
else:
  closest_pose_o = self._closest_pose_BAT(xo_BAT[self._BAT_to_perturb])

report = ''
for t in range(ntrials):
  # Choose a pose to dart towards
  dart_towards = closest_pose_o
  while dart_towards==closest_pose_o:
    dart_towards = np.random.choice(len(self.weights), p=self.weights)
  # Generate a trial move
  xn_BAT = np.copy(xo_BAT)
  xn_BAT[self._BAT_to_perturb] = xo_BAT[self._BAT_to_perturb] + self.darts[closest_pose_o][dart_towards]
  
  # Check that the trial move is closest to dart_towards
  if self.extended:
    xn_Cartesian = self._BAT_util.Cartesian(xn_BAT)
    closest_pose_n = self._closest_pose_Cartesian(\
      xn_Cartesian[self.molecule.heavy_atoms,:])
    if (closest_pose_n!=dart_towards):
      continue
  else:
    closest_pose_n = self._closest_pose_BAT(xn_BAT[self._BAT_to_perturb])
    if (closest_pose_n!=dart_towards):
      continue
    xn_Cartesian = self._BAT_util.Cartesian(xn_BAT)

  # Determine energy of new state
  self.universe.setConfiguration(Configuration(self.universe, xn_Cartesian))
  en = self.universe.energy()
  report += 'Attempting move from near pose %d with energy %f to pose %d with energy %f. '%(closest_pose_o,eo,closest_pose_n,en)

  # Accept or reject the trial move
  if (abs(en-eo)<1000) and \
      ((en<eo) or (np.random.random()<np.exp(-(en-eo)/RT))):
    xo_Cartesian = xn_Cartesian
    xo_BAT = xn_BAT
    eo = 1.*en
    closest_pose_o = closest_pose_n
    acc += 1
    report += 'Accepted.\n'
  else:
    report += 'Rejected.\n'

print report
self.universe.setConfiguration(Configuration(self.universe, xo_Cartesian))

#  confs = self.confs + [xo_Cartesian, xn_Cartesian]
#  import AlGDock.IO
#  IO_dcd = AlGDock.IO.dcd(self.molecule)
#  IO_dcd.write('confs.dcd', confs)
#  self._BAT_util.showMolecule(dcdFN='confs.dcd')
#  os.remove('confs.dcd')