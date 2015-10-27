# The example is 1of6

run_types = ['cool']

import AlGDock.BindingPMF_plots
import os, shutil, glob

self = None
for run_type in run_types:
  del self
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
    seeds_per_state=10, steps_per_seed=200, darts_per_seed=0, \
    sweeps_per_cycle=25, attempts_per_sweep=100, \
    steps_per_sweep=50, darts_per_sweep=0, \
    cool_repX_cycles=3, dock_repX_cycles=4, \
    site='Sphere', site_center=[1.74395, 1.74395, 1.74395], site_max_R=0.6, \
    site_density=10., \
    phases=['NAMD_Gas','NAMD_GBSA'], \
    cores=-1, \
    rmsd=True)
  self._run(run_type)

process = 'cool'
from AlGDock.BindingPMF_plots import *

if not process in ['dock','cool']:
  raise Exception('Process must be dock or cool')
# GMC
def gMC_initial_setup():
  """
  Initialize BAT converter object.
  Decide which internal coord to crossover. Here, only the soft torsions will be crossovered.
  Produce a list of replica (state) index pairs to be swaped. Only Neighbor pairs will be swaped.
  Assume that self.universe, self.molecule and K (number of states) exist
  as global variables when the function is called.
  """
  from AlGDock.RigidBodies import identifier
  import itertools
  BAT_converter = identifier( self.universe, self.molecule )
  BAT = BAT_converter.BAT( extended = True )
  # this assumes that the torsional angles are stored in the tail of BAT
  softTorsionId = [ i + len(BAT) - BAT_converter.ntorsions for i in BAT_converter._softTorsionInd ]
  torsions_to_crossover = []
  for i in range(1, len(softTorsionId) ):
    combinations = itertools.combinations( softTorsionId, i )
    for c in combinations:
      torsions_to_crossover.append( list(c) )
  #
  BAT_converter.BAT_to_crossover = torsions_to_crossover
  if len( BAT_converter.BAT_to_crossover ) == 0:
    self.tee('  GMC No BAT to crossover')
  state_indices = range( K )
  state_indices_to_swap = zip( state_indices[0::2], state_indices[1::2] ) + \
                  zip( state_indices[1::2], state_indices[2::2] )
  #
  return BAT_converter, state_indices_to_swap
#
def do_gMC( nr_attempts, BAT_converter, state_indices_to_swap, torsion_threshold ):
  """
  Assume self.universe, confs, lambdas, state_inds, inv_state_inds exist as global variables
  when the function is called.
  If at least one of the torsions in the combination chosen for an crossover attempt
  changes more than torsion_threshold, the crossover will be attempted.
  The function will update confs.
  It returns the number of attempts and the number of accepted moves.
  """
  if nr_attempts < 0:
    raise Exception('Number of attempts must be nonnegative!')
  if torsion_threshold < 0.:
    raise Exception('Torsion threshold must be nonnegative!')
  #
  if len( BAT_converter.BAT_to_crossover ) == 0:
    return 0., 0.
  #
  from random import randrange
  # get reduced energies and BAT for all configurations in confs
  BATs = []
  energies = np.zeros( K, dtype = float )
  for c_ind in range(K):
    s_ind = state_inds[ c_ind ]
    self.universe.setConfiguration( Configuration( self.universe, confs[c_ind] ) )
    BATs.append( np.array( BAT_converter.BAT( extended = True ) , dtype = float ) )
    self._set_universe_evaluator( lambdas[ s_ind ] )
    reduced_e = self.universe.energy() / ( R*lambdas[ s_ind ]['T'] )
    energies[ c_ind ] = reduced_e
  #
  nr_sets_of_torsions = len( BAT_converter.BAT_to_crossover )
  #
  attempt_count , acc_count = 0 , 0
  sweep_count = 0
  while True:
    sweep_count += 1
    if (sweep_count * K) > (1000 * nr_attempts):
      self.tee('  GMC Sweep too many times, but few attempted. Consider reducing torsion_threshold.')
      return attempt_count, acc_count
    #
    for state_pair in state_indices_to_swap:
      conf_ind_k0 = inv_state_inds[ state_pair[0] ]
      conf_ind_k1 = inv_state_inds[ state_pair[1] ]
      # check if it should attempt for this pair of states
      ran_set_torsions = BAT_converter.BAT_to_crossover[ randrange( nr_sets_of_torsions ) ]
      do_crossover = np.any(np.abs(BATs[conf_ind_k0][ran_set_torsions] - BATs[conf_ind_k1][ran_set_torsions]) >= torsion_threshold)
      if do_crossover:
        attempt_count += 1
        # BAT and reduced energies before crossover
        BAT_k0_be = copy.deepcopy( BATs[conf_ind_k0] )
        BAT_k1_be = copy.deepcopy( BATs[conf_ind_k1] )
        e_k0_be = energies[conf_ind_k0]
        e_k1_be = energies[conf_ind_k1]
        # BAT after crossover
        BAT_k0_af = copy.deepcopy( BAT_k0_be )
        BAT_k1_af = copy.deepcopy( BAT_k1_be )
        for index in ran_set_torsions:
          tmp = BAT_k0_af[ index ]
          BAT_k0_af[ index ] = BAT_k1_af[ index ]
          BAT_k1_af[ index ] = tmp
        # Cartesian coord and reduced energies after crossover.
        BAT_converter.Cartesian( BAT_k0_af )
        self._set_universe_evaluator( lambdas[ state_pair[0] ] )
        e_k0_af = self.universe.energy() / ( R*lambdas[ state_pair[0] ]['T'] )
        conf_k0_af = copy.deepcopy( self.universe.configuration().array )
        #
        BAT_converter.Cartesian( BAT_k1_af )
        self._set_universe_evaluator( lambdas[ state_pair[1] ] )
        e_k1_af = self.universe.energy() / ( R*lambdas[ state_pair[1] ]['T'] )
        conf_k1_af = copy.deepcopy( self.universe.configuration().array )
        #
        de = ( e_k0_be - e_k0_af ) + ( e_k1_be - e_k1_af )
        # update confs, energies, BATS
        if (de > 0) or ( np.random.uniform() < np.exp(de) ):
          acc_count += 1
          confs[conf_ind_k0] = conf_k0_af
          confs[conf_ind_k1] = conf_k1_af
          #
          energies[conf_ind_k0] = e_k0_af
          energies[conf_ind_k1] = e_k1_af
          #
          BATs[conf_ind_k0] = BAT_k0_af
          BATs[conf_ind_k1] = BAT_k1_af
        #
        if attempt_count == nr_attempts:
          return attempt_count, acc_count

if process=='cool':
  terms = ['MM']
else:
  terms = ['MM','site','misc'] + self._scalables

cycle = getattr(self,'_%s_cycle'%process)
confs = self.confs[process]['replicas']
lambdas = getattr(self,process+'_protocol')

# A list of pairs of replica indicies
K = len(lambdas)
pairs_to_swap = []
for interval in range(1,min(5,K)):
  lower_inds = []
  for lowest_index in range(interval):
    lower_inds += range(lowest_index,K-interval,interval)
  upper_inds = np.array(lower_inds) + interval
  pairs_to_swap += zip(lower_inds,upper_inds)

from repX import attempt_swaps

# Setting the force field will load grids
# before multiple processes are spawned
for k in range(K):
  self._set_universe_evaluator(lambdas[k])

storage = {}
for var in ['confs','state_inds','energies']:
  storage[var] = []

cycle_start_time = time.time()

if self._cores>1:
  # Multiprocessing setup
  m = multiprocessing.Manager()
  task_queue = m.Queue()
  done_queue = m.Queue()

# GMC
do_gMC = self.params[process]['GMC_attempts'] > 0
if do_gMC:
  self.tee('  Using GMC for %s' %process)
  nr_gMC_attempts = K * self.params[process]['GMC_attempts']
  torsion_threshold = self.params[process]['GMC_tors_threshold']
  gMC_attempt_count = 0
  gMC_acc_count     = 0
  time_gMC = 0.0
  BAT_converter, state_indices_to_swap = gMC_initial_setup()

# Do replica exchange
time_ExternalMC = 0.0
time_SmartDarting = 0.0
time_repX = 0.0
state_inds = range(K)
inv_state_inds = range(K)
for sweep in range(self.params[process]['sweeps_per_cycle']):
  E = {}
  for term in terms:
    E[term] = np.zeros(K, dtype=float)
  E['acc_ExternalMC'] = np.zeros(K, dtype=float)
  E['acc_SmartDarting'] = np.zeros(K, dtype=float)
  acc_Sampler = np.zeros(K, dtype=float)
  # Sample within each state
  if self._cores>1:
    for k in range(K):
      task_queue.put((confs[k], process, lambdas[state_inds[k]], False, k))
    for p in range(self._cores):
      task_queue.put('STOP')
    processes = [multiprocessing.Process(target=self._sim_one_state_worker, \
        args=(task_queue, done_queue)) for p in range(self._cores)]
    for p in processes:
      p.start()
    for p in processes:
      p.join()
    unordered_results = [done_queue.get() for k in range(K)]
    results = sorted(unordered_results, key=lambda d: d['reference'])
    for p in processes:
      p.terminate()
  else:
    # Single process code
    results = [self._sim_one_state(confs[k], process, \
        lambdas[state_inds[k]], False, k) for k in range(K)]
  # GMC
  if do_gMC:
    time_start_gMC = time.time()
    att_count, acc_count = do_gMC( nr_gMC_attempts, BAT_converter, state_indices_to_swap, torsion_threshold )
    gMC_attempt_count += att_count
    gMC_acc_count     += acc_count
    time_gMC =+ ( time.time() - time_start_gMC )
  # Store results
  for k in range(K):
    if 'acc_ExternalMC' in results[k].keys():
      E['acc_ExternalMC'][k] = results[k]['acc_ExternalMC']
      time_ExternalMC += results[k]['time_ExternalMC']
    if 'acc_SmartDarting' in results[k].keys():
      E['acc_SmartDarting'][k] = results[k]['acc_SmartDarting']
      time_SmartDarting += results[k]['time_SmartDarting']
    confs[k] = results[k]['confs']
    if process=='cool':
        E['MM'][k] = results[k]['Etot']
    acc_Sampler[k] += results[k]['acc_Sampler']
  if process=='dock':
    E = self._calc_E(confs, E) # Get energies for scalables
    # Get rmsd values
    if self.params['dock']['rmsd'] is not False:
      E['rmsd'] = np.array([np.sqrt(((confs[k][self.molecule.heavy_atoms,:] - \
        self.confs['rmsd'])**2).sum()/self.molecule.nhatoms) for k in range(K)])
  # Calculate u_ij (i is the replica, and j is the configuration),
  #    a list of arrays
  (u_ij,N_k) = self._u_kln(E, [lambdas[state_inds[c]] for c in range(K)])
  # Do the replica exchange
  repX_start_time = time.time()
  (state_inds, inv_state_inds) = \
    attempt_swaps(state_inds, inv_state_inds, u_ij, pairs_to_swap, \
      self.params[process]['attempts_per_sweep'])
  time_repX += (time.time()-repX_start_time)
  # Store data in local variables
  storage['confs'].append(list(confs))
  storage['state_inds'].append(list(state_inds))
  storage['energies'].append(copy.deepcopy(E))

# GMC
if do_gMC:
  self.tee('  {0}/{1} crossover attempts ({2:.3g}) accepted in {3}'.format(\
    gMC_acc_count, gMC_attempt_count, \
    float(gMC_acc_count)/float(gMC_attempt_count) \
      if gMC_attempt_count > 0 else 0, \
    HMStime(time_gMC)))

# Estimate relaxation time from autocorrelation
state_inds = np.array(storage['state_inds'])
tau_ac = pymbar.timeseries.integratedAutocorrelationTimeMultiple(state_inds.T)
# There will be at least per_independent and up to sweeps_per_cycle saved samples
# max(int(np.ceil((1+2*tau_ac)/per_independent)),1) is the minimum stride,
# which is based on per_independent samples per autocorrelation time.
# max(self.params['dock']['sweeps_per_cycle']/per_independent)
# is the maximum stride, which gives per_independent samples if possible.
per_independent = self.params[process]['snaps_per_independent']
stride = min(max(int(np.ceil((1+2*tau_ac)/per_independent)),1), \
             max(int(np.ceil(self.params[process]['sweeps_per_cycle']/per_independent)),1))

store_indicies = np.array(\
  range(min(stride-1,self.params[process]['sweeps_per_cycle']-1), \
  self.params[process]['sweeps_per_cycle'], stride), dtype=int)
nsaved = len(store_indicies)

self.tee("  storing %d configurations for %d replicas"%(nsaved, len(confs)) + \
  " in cycle %d"%cycle + \
  " (tau_ac=%f)"%(tau_ac))
self.tee("  with %s for external MC"%(HMStime(time_ExternalMC)) + \
  " and %s for smart darting"%(HMStime(time_SmartDarting)) + \
  " and %s for replica exchange"%(HMStime(time_repX)) + \
  " in " + HMStime(time.time()-cycle_start_time))


# Get indicies for storing global variables
inv_state_inds = np.zeros((nsaved,K),dtype=int)
for snap in range(nsaved):
  state_inds = storage['state_inds'][store_indicies[snap]]
  for state in range(K):
    inv_state_inds[snap][state_inds[state]] = state

# Reorder energies and replicas for storage
if process=='dock':
  terms.append('acc_ExternalMC') # Make sure to save the acceptance probability
  if self.params['dock']['rmsd'] is not False:
    terms.append('rmsd') # Make sure to save the rmsd
terms.append('acc_SmartDarting')
Es = []
for state in range(K):
  E_state = {}
  if state==0:
    E_state['repXpath'] = storage['state_inds']
    E_state['acc_Sampler'] = acc_Sampler
  for term in terms:
    E_state[term] = np.array([storage['energies'][store_indicies[snap]][term][inv_state_inds[snap][state]] for snap in range(nsaved)])
  Es.append([E_state])

self.confs[process]['replicas'] = \
  [storage['confs'][store_indicies[-1]][inv_state_inds[-1][state]] \
   for state in range(K)]

for state in range(K):
  getattr(self,process+'_Es')[state].append(Es[state][0])

for state in range(K):
  if self.params[process]['keep_intermediate'] or \
      ((process=='cool') and (state==0)) or \
      (state==(K-1)):
    confs = [storage['confs'][store_indicies[snap]][inv_state_inds[snap][state]] for snap in range(nsaved)]
    self.confs[process]['samples'][state].append(confs)
  else:
    self.confs[process]['samples'][state].append([])

if self.params[process]['darts_per_sweep']>0:
  self._set_universe_evaluator(getattr(self,process+'_protocol')[-1])
  confs_SmartDarting = [np.copy(conf) \
    for conf in self.confs[process]['samples'][state][-1]]
  self.tee(self.sampler[process+'_SmartDarting'].set_confs(\
    confs_SmartDarting, append=True))

setattr(self,'_%s_cycle'%process,cycle + 1)
