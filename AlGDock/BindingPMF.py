#!/usr/bin/env python

# TODO: Free energy of external confinement for poseBPMFs

import os
import cPickle as pickle
import gzip
import copy

from AlGDock.IO import load_pkl_gz
from AlGDock.IO import write_pkl_gz
from AlGDock.logger import NullDevice

import sys
import time
import numpy as np

from collections import OrderedDict

from AlGDock import dictionary_tools
from AlGDock import path_tools

import MMTK
import MMTK.Units
from MMTK.ParticleProperties import Configuration
from MMTK.ForceFields import ForceField

import Scientific
try:
  from Scientific._vector import Vector
except:
  from Scientific.Geometry.VectorModule import Vector

import pymbar.timeseries

import multiprocessing
from multiprocessing import Process

try:
  import requests  # for downloading additional files
except:
  print '  no requests module for downloading additional files'

# For profiling. Unnecessary for normal execution.
# from memory_profiler import profile

#############
# Constants #
#############

R = 8.3144621 * MMTK.Units.J / MMTK.Units.mol / MMTK.Units.K

scalables = ['OBC', 'sLJr', 'sELE', 'LJr', 'LJa', 'ELE']

# In APBS, minimum ratio of PB grid length to maximum dimension of solute
LFILLRATIO = 4.0  # For the ligand
RFILLRATIO = 2.0  # For the receptor/complex

DEBUG = False


def HMStime(s):
  """
  Given the time in seconds, an appropriately formatted string.
  """
  if s < 60.:
    return '%.2f s' % s
  elif s < 3600.:
    return '%d:%.2f' % (int(s / 60 % 60), s % 60)
  else:
    return '%d:%d:%.2f' % (int(s / 3600), int(s / 60 % 60), s % 60)


##############
# Main Class #
##############


class BPMF:
  def __init__(self, **kwargs):
    """Parses the input arguments and runs the requested CD calculation"""

    #         mod_path = os.path.join(os.path.dirname(a.__file__), 'BindingPMF.py')
    #         print """###########
    # # AlGDock #
    # ###########
    # Molecular CD with adaptively scaled alchemical interaction grids
    #
    # in {0}
    # last modified {1}
    #     """.format(mod_path, time.ctime(os.path.getmtime(mod_path)))

    from AlGDock.argument_parser import SimulationArguments
    self.args = SimulationArguments(**kwargs)

    from AlGDock.simulation_data import SimulationData
    self.data = {}
    self.data['BC'] = SimulationData(self.args.dir['BC'], 'BC', \
      self.args.params['CD']['pose'])
    self.data['CD'] = SimulationData(self.args.dir['CD'], 'CD', \
      self.args.params['CD']['pose'])

    if not 'run_type' in kwargs.keys():
      kwargs['run_type'] = None

    from AlGDock.logger import Logger
    self.log = Logger(self.args, \
      max_time=kwargs['max_time'], run_type=kwargs['run_type'])

    self.T_HIGH = self.args.params['BC']['T_HIGH']
    self.T_TARGET = self.args.params['BC']['T_TARGET']

    self._setup()

    print '\n*** Simulation parameters and constants ***'
    for p in ['BC', 'CD']:
      print '\nfor %s:' % p
      print dictionary_tools.dict_view(self.args.params[p])[:-1]

    self._run(kwargs['run_type'])

  def _setup(self):
    """Creates an MMTK InfiniteUniverse and adds the ligand"""

    from AlGDock.topology import Topology
    self.top = Topology(self.args)
    self.top_RL = Topology(self.args, includeReceptor=True)

    # Initialize rmsd calculation function
    from AlGDock.RMSD import hRMSD
    self.get_rmsds = hRMSD(self.args.FNs['prmtop']['L'], \
      self.top.inv_prmtop_atom_order_L)

    # Obtain reference pose
    if self.data['CD'].pose > -1:
      if ('starting_poses' in self.data['CD'].confs.keys()) and \
         (self.data['CD'].confs['starting_poses'] is not None):
        starting_pose = np.copy(self.data['CD'].confs['starting_poses'][0])
      else:
        (confs, Es) = self._get_confs_to_rescore(site=False, \
          minimize=False, sort=False)
        if self.args.params['CD']['pose'] < len(confs):
          starting_pose = np.copy(confs[self.args.params['CD']['pose']])
          self.data['CD'].confs['starting_poses'] = [np.copy(starting_pose)]
        else:
          self._clear('CD')
          self._store_infinite_f_RL()
          raise Exception('Pose index greater than number of poses')
    else:
      starting_pose = None

    from AlGDock.system import System
    self.system = System(self.args,
                         self.log,
                         self.top,
                         self.top_RL,
                         starting_pose=starting_pose)

    # Measure the binding site
    if (self.args.params['CD']['site'] == 'Measure'):
      self.args.params['CD']['site'] = 'Sphere'
      if self.args.params['CD']['site_measured'] is not None:
        (self.args.params['CD']['site_max_R'],self.args.params['CD']['site_center']) = \
          self.args.params['CD']['site_measured']
      else:
        print '\n*** Measuring the binding site ***'
        self.system.setParams(
          self.system.paramsFromAlpha(1.0, 'CD', site=False))
        (confs, Es) = self._get_confs_to_rescore(site=False, minimize=True)
        if len(confs) > 0:
          # Use the center of mass for configurations
          # within 20 RT of the lowest energy
          cutoffE = Es['total'][-1] + 20 * (R * self.T)
          coms = []
          for (conf, E) in reversed(zip(confs, Es['total'])):
            if E <= cutoffE:
              self.top.universe.setConfiguration(
                Configuration(self.top.universe, conf))
              coms.append(np.array(self.top.universe.centerOfMass()))
            else:
              break
          print '  %d configurations fit in the binding site' % len(coms)
          coms = np.array(coms)
          center = (np.min(coms, 0) + np.max(coms, 0)) / 2
          max_R = max(
            np.ceil(np.max(np.sqrt(np.sum(
              (coms - center)**2, 1))) * 10.) / 10., 0.6)
          self.args.params['CD']['site_max_R'] = max_R
          self.args.params['CD']['site_center'] = center
          self.top.universe.setConfiguration(
            Configuration(self.top.universe, confs[-1]))
        if ((self.args.params['CD']['site_max_R'] is None) or \
            (self.args.params['CD']['site_center'] is None)):
          raise Exception('No binding site parameters!')
        else:
          self.args.params['CD']['site_measured'] = \
            (self.args.params['CD']['site_max_R'], \
             self.args.params['CD']['site_center'])

    # Read the reference ligand and receptor coordinates
    import AlGDock.IO
    IO_crd = AlGDock.IO.crd()
    if self.args.FNs['inpcrd']['R'] is not None:
      if os.path.isfile(self.args.FNs['inpcrd']['L']):
        lig_crd = IO_crd.read(self.args.FNs['inpcrd']['L'], multiplier=0.1)
      self.data['CD'].confs['receptor'] = IO_crd.read(\
        self.args.FNs['inpcrd']['R'], multiplier=0.1)
    elif self.args.FNs['inpcrd']['RL'] is not None:
      complex_crd = IO_crd.read(self.args.FNs['inpcrd']['RL'], multiplier=0.1)
      lig_crd = complex_crd[self.top_RL.L_first_atom:self.top_RL.L_first_atom + \
        self.top.universe.numberOfAtoms(),:]
      self.data['CD'].confs['receptor'] = np.vstack(\
        (complex_crd[:self.top_RL.L_first_atom,:],\
         complex_crd[self.top_RL.L_first_atom + self.top.universe.numberOfAtoms():,:]))
    elif self.args.FNs['inpcrd']['L'] is not None:
      self.data['CD'].confs['receptor'] = None
      if os.path.isfile(self.args.FNs['inpcrd']['L']):
        lig_crd = IO_crd.read(self.args.FNs['inpcrd']['L'], multiplier=0.1)
    else:
      lig_crd = None

    if lig_crd is not None:
      self.data['CD'].confs['ligand'] = lig_crd[self.top.
                                                inv_prmtop_atom_order_L, :]
      self.top.universe.setConfiguration(\
        Configuration(self.top.universe,self.data['CD'].confs['ligand']))
      if self.top_RL.universe is not None:
        self.top_RL.universe.setConfiguration(\
          Configuration(self.top_RL.universe, \
          np.vstack((self.data['CD'].confs['receptor'],self.data['CD'].confs['ligand']))))

    if self.args.params['CD']['rmsd'] is not False:
      if self.args.params['CD']['rmsd'] is True:
        if lig_crd is not None:
          rmsd_crd = lig_crd[self.top.inv_prmtop_atom_order_L, :]
        else:
          raise Exception('Reference structure for rmsd calculations unknown')
      else:
        rmsd_crd = IO_crd.read(self.args.params['CD']['rmsd'], \
          natoms=self.top.universe.numberOfAtoms(), multiplier=0.1)
        rmsd_crd = rmsd_crd[self.top.inv_prmtop_atom_order_L, :]
      self.data['CD'].confs['rmsd'] = rmsd_crd
      self.get_rmsds.set_ref_configuration(self.data['CD'].confs['rmsd'])

    # TODO: Remove after postprocessin
    # Locate programs for postprocessing
    all_phases = self.args.params['CD']['phases'] + self.args.params['BC'][
      'phases']
    self._load_programs(all_phases)

    # Determine APBS grid spacing
    if 'APBS_PBSA' in self.args.params['CD']['phases'] or \
       'APBS_PBSA' in self.args.params['BC']['phases']:
      self._get_APBS_grid_spacing()

    # Determines receptor electrostatic size
    if np.array([p.find('ALPB') > -1 for p in all_phases]).any():
      self.elsize = self._get_elsize()

    # If configurations are being rescored, start with a docked structure
    (confs, Es) = self._get_confs_to_rescore(site=False, minimize=False)
    if len(confs) > 0:
      self.top.universe.setConfiguration(
        Configuration(self.top.universe, confs[-1]))

    # Samplers may accept the following options:
    # steps - number of MD steps
    # T - temperature in K
    # delta_t - MD time step
    # normalize - normalizes configurations
    # adapt - uses an adaptive time step

    from AlGDock.simulation_iterator import SimulationIterator
    self.iterator = SimulationIterator(self.args, self.top, self.system)

    # Load progress
    self._postprocess(readOnly=True)
    self.calc_f_L(readOnly=True)
    self.calc_f_RL(readOnly=True)

    if self.args.random_seed > 0:
      np.random.seed(self.args.random_seed)

  def _run(self, run_type):
    self.log.recordStart('run')
    self.log.run_type = run_type
    if run_type=='configuration_energies' or \
       run_type=='minimized_configuration_energies':
      self.configuration_energies(\
        minimize = (run_type=='minimized_configuration_energies'), \
        max_confs = 50)
    elif run_type == 'store_params':
      self.save('BC', keys=['progress'])
      self.save('CD', keys=['progress'])
    elif run_type == 'initial_BC':
      self.initial_BC()
    elif run_type == 'BC':  # Sample the BC process
      self.sim_process('BC')
      self._postprocess([('BC', -1, -1, 'L')])
      self.calc_f_L()
    elif run_type == 'initial_CD':
      self.initial_CD()
    elif run_type == 'CD':  # Sample the CD process
      self.sim_process('CD')
      self._postprocess()
      self.calc_f_RL()
      # self.targeted_FEP()
    elif run_type == 'timed':  # Timed replica exchange sampling
      BC_complete = self.sim_process('BC')
      if BC_complete:
        pp_complete = self._postprocess([('BC', -1, -1, 'L')])
        if pp_complete:
          self.calc_f_L()
          CD_complete = self.sim_process('CD')
          if CD_complete:
            pp_complete = self._postprocess()
            if pp_complete:
              self.calc_f_RL()
              # self.targeted_FEP()
    elif run_type == 'timed_BC':  # Timed BC only
      BC_complete = self.sim_process('BC')
      if BC_complete:
        pp_complete = self._postprocess([('BC', -1, -1, 'L')])
        if pp_complete:
          self.calc_f_L()
    elif run_type == 'timed_CD':  # Timed CD only
      CD_complete = self.sim_process('CD')
      if CD_complete:
        pp_complete = self._postprocess()
        if pp_complete:
          self.calc_f_RL()
          # self.targeted_FEP()
    elif run_type == 'postprocess':  # Postprocessing
      self._postprocess()
    elif run_type == 'redo_postprocess':
      self._postprocess(redo_CD=True)
    elif run_type == 'redo_pose_prediction':
      self.calc_f_RL(readOnly=True)
      # Predict native pose
      if self.args.params['CD']['pose'] == -1:
        (self.stats_RL['pose_inds'], self.stats_RL['scores']) = \
          self._get_pose_prediction()
        f_RL_FN = os.path.join(self.args.dir['CD'], 'f_RL.pkl.gz')
        self.log.tee(
          write_pkl_gz(f_RL_FN, (self.f_L, self.stats_RL, self.f_RL, self.B)))
      # self.targeted_FEP()
    elif (run_type == 'free_energies') or (run_type == 'redo_free_energies'):
      self.calc_f_L(redo=(run_type == 'redo_free_energies'))
      self.calc_f_RL(redo=(run_type == 'redo_free_energies'))
      # self.targeted_FEP()
    elif run_type == 'all':
      self.sim_process('BC')
      self._postprocess([('BC', -1, -1, 'L')])
      self.calc_f_L()
      self.sim_process('CD')
      self._postprocess()
      self.calc_f_RL()
      # self.targeted_FEP()
    elif run_type == 'render_docked':
      # For 4 figures
      # 1002*4/600. = 6.68 in at 600 dpi
      #  996*4/600. = 6.64 in at 600 dpi
      view_args = {'axes_off':True, 'size':[996,996], 'scale_by':0.80, \
                   'render':'TachyonInternal'}
      if hasattr(self, '_view_args_rotate_matrix'):
        view_args['rotate_matrix'] = getattr(self, '_view_args_rotate_matrix')
      self.show_samples(prefix='docked', \
        show_ref_ligand=True, show_starting_pose=True, \
        show_receptor=True, save_image=True, execute=True, quit=True, \
        view_args=view_args)
      if self.args.params['CD']['pose'] == -1:
        (self.stats_RL['pose_inds'], self.stats_RL['scores']) = \
          self._get_pose_prediction()
        self.show_pose_prediction(score='grid_fe_u',
          show_ref_ligand=True, show_starting_pose=False, \
          show_receptor=True, save_image=True, execute=True, quit=True, \
          view_args=view_args)
        self.show_pose_prediction(score='OpenMM_OBC2_fe_u',
          show_ref_ligand=True, show_starting_pose=False, \
          show_receptor=True, save_image=True, execute=True, quit=True, \
          view_args=view_args)
    elif run_type == 'render_intermediates':
      view_args = {'axes_off':True, 'size':[996,996], 'scale_by':0.80, \
                   'render':'TachyonInternal'}
      if hasattr(self, '_view_args_rotate_matrix'):
        view_args['rotate_matrix'] = getattr(self, '_view_args_rotate_matrix')
#      self.render_intermediates(\
#        movie_name=os.path.join(self.args.dir['CD'],'CD-intermediates.gif'), \
#        view_args=view_args)
      self.render_intermediates(nframes=8, view_args=view_args)
    elif run_type == 'clear_intermediates':
      for process in ['BC', 'CD']:
        print 'Clearing intermediates for ' + process
        for state_ind in range(1,
                               len(self.data[process].confs['samples']) - 1):
          for cycle_ind in range(
              len(self.data[process].confs['samples'][state_ind])):
            self.data[process].confs['samples'][state_ind][cycle_ind] = []
        self.save(process)
    if run_type is not None:
      print "\nElapsed time for execution of %s: %s" % (
        run_type, HMStime(self.log.timeSince('run')))

  ###########
  # BC #
  ###########
  def initial_BC(self):
    """
    Warms the ligand from self.T_TARGET to self.T_HIGH

    Intermediate thermodynamic states are chosen such that
    thermodynamic length intervals are approximately constant.
    Configurations from each state are subsampled to seed the next simulation.
    """

    if (len(self.data['BC'].protocol) >
        0) and (self.data['BC'].protocol[-1]['crossed']):
      return  # Initial BC is already complete

    self.log.recordStart('BC')

    from AlGDock.ligand_preparation import LigandPreparation
    seeds = LigandPreparation(self.args, self.log, self.top, self.system,
                              self._get_confs_to_rescore, self.iterator,
                              self.data).run('BC')

    from AlGDock.initialization import Initialization
    Initialization(self.args, self.log, self.top, self.system,
                  self.iterator, self.data, self.save, self._u_kln).run('BC', seeds)

    return True

  def calc_f_L(self, readOnly=False, do_solvation=True, redo=False):
    """
    Calculates ligand-specific free energies:
    1. reduced free energy of BC the ligand
       from self.T_HIGH to self.T_TARGET
    2. solvation free energy of the ligand using single-step
       free energy perturbation
    redo does not do anything now; it is an option for debugging
    """
    # Initialize variables as empty lists or by loading data
    f_L_FN = os.path.join(self.args.dir['BC'], 'f_L.pkl.gz')
    dat = load_pkl_gz(f_L_FN)
    if dat is not None:
      (self.stats_L, self.f_L) = dat
    else:
      self.stats_L = dict(\
        [(item,[]) for item in ['equilibrated_cycle','mean_acc']])
      self.stats_L['protocol'] = self.data['BC'].protocol
      self.f_L = dict([(key,[]) for key in ['BC_MBAR'] + \
        [phase+'_solv' for phase in self.args.params['BC']['phases']]])
    if readOnly or self.data['BC'].protocol == []:
      return

    K = len(self.data['BC'].protocol)

    # Make sure all the energies are available
    for c in range(self.data['BC'].cycle):
      if len(self.data['BC'].Es[-1][c].keys()) == 0:
        self.log.tee("  skipping the BC free energy calculation")
        return

    start_string = "\n>>> Ligand free energy calculations, starting at " + \
      time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "\n"
    self.log.recordStart('free energy')

    # Store stats_L internal energies
    self.stats_L['u_K_sampled'] = \
      [self._u_kln([self.data['BC'].Es[-1][c]],[self.data['BC'].protocol[-1]]) \
        for c in range(self.data['BC'].cycle)]
    self.stats_L['u_KK'] = \
      [np.sum([self._u_kln([self.data['BC'].Es[k][c]],[self.data['BC'].protocol[k]]) \
        for k in range(len(self.data['BC'].protocol))],0) \
          for c in range(self.data['BC'].cycle)]

    self.stats_L['equilibrated_cycle'] = self._get_equilibrated_cycle('BC')

    # Calculate BC free energies that have not already been calculated,
    # in units of RT
    updated = False
    for c in range(len(self.f_L['BC_MBAR']), self.data['BC'].cycle):
      if not updated:
        self.log.set_lock('BC')
        if do_solvation:
          self.log.tee(start_string)
        updated = True

      fromCycle = self.stats_L['equilibrated_cycle'][c]
      toCycle = c + 1

      # BC free energy
      BC_Es = []
      for BC_Es_state in self.data['BC'].Es:
        BC_Es.append(BC_Es_state[fromCycle:toCycle])
      (u_kln, N_k) = self._u_kln(BC_Es, self.data['BC'].protocol)
      MBAR = self._run_MBAR(u_kln, N_k)[0]
      self.f_L['BC_MBAR'].append(MBAR)

      # Average acceptance probabilities
      BC_mean_acc = np.zeros(K - 1)
      for k in range(0, K - 1):
        (u_kln, N_k) = self._u_kln(BC_Es[k:k + 2],
                                   self.data['BC'].protocol[k:k + 2])
        N = min(N_k)
        acc = np.exp(-u_kln[0, 1, :N] - u_kln[1, 0, :N] + u_kln[0, 0, :N] +
                     u_kln[1, 1, :N])
        BC_mean_acc[k] = np.mean(np.minimum(acc, np.ones(acc.shape)))
      self.stats_L['mean_acc'].append(BC_mean_acc)

      self.log.tee("  calculated BC free energy of %.2f RT "%(\
                  self.f_L['BC_MBAR'][-1][-1])+\
               "using cycles %d to %d"%(fromCycle, c))

    if not do_solvation:
      if updated:
        if not self.log.run_type.startswith('timed'):
          write_pkl_gz(f_L_FN, (self.stats_L, self.f_L))
        self.log.clear_lock('BC')
      return True

    # Make sure postprocessing is complete
    pp_complete = self._postprocess([('BC', -1, -1, 'L')])
    if not pp_complete:
      return False

    # Store stats_L internal energies
    for phase in self.args.params['BC']['phases']:
      self.stats_L['u_K_'+phase] = \
        [self.data['BC'].Es[-1][c]['L'+phase][:,-1]/(R*self.T_TARGET) \
          for c in range(self.data['BC'].cycle)]

    # Calculate solvation free energies that have not already been calculated,
    # in units of RT
    for phase in self.args.params['BC']['phases']:
      if not phase + '_solv' in self.f_L:
        self.f_L[phase + '_solv'] = []
      if not 'mean_' + phase in self.f_L:
        self.f_L['mean_' + phase] = []

      for c in range(len(self.f_L[phase + '_solv']), self.data['BC'].cycle):
        if not updated:
          self.log.set_lock('BC')
          self.log.tee(start_string)
          updated = True

        fromCycle = self.stats_L['equilibrated_cycle'][c]
        toCycle = c + 1

        if not ('L' + phase) in self.data['BC'].Es[-1][c].keys():
          raise Exception('L%s energies not found in cycle %d' % (phase, c))

        # Arbitrarily, solvation is the
        # 'forward' direction and desolvation the 'reverse'
        u_L = np.concatenate([self.data['BC'].Es[-1][n]['L'+phase] \
          for n in range(fromCycle,toCycle)])/(R*self.T_TARGET)
        u_sampled = np.concatenate(\
          [self._u_kln([self.data['BC'].Es[-1][c]],[self.data['BC'].protocol[-1]]) \
            for c in range(fromCycle,toCycle)])
        du_F = (u_L[:, -1] - u_sampled)
        min_du_F = min(du_F)
        w_L = np.exp(-du_F + min_du_F)
        f_L_solv = -np.log(np.mean(w_L)) + min_du_F
        mean_u_phase = np.sum(u_L[:, -1] * w_L) / np.sum(w_L)

        self.f_L[phase + '_solv'].append(f_L_solv)
        self.f_L['mean_' + phase].append(mean_u_phase)
        self.log.tee("  calculated " + phase + " solvation free energy of " + \
                 "%.5g RT "%(f_L_solv) + \
                 "using cycles %d to %d"%(fromCycle, toCycle-1))

    if updated:
      self.log.tee(write_pkl_gz(f_L_FN, (self.stats_L, self.f_L)))
      self.log.tee("\nElapsed time for free energy calculation: " + \
        HMStime(self.log.timeSince('free energy')))
      self.log.clear_lock('BC')
    return True

  ###########
  # Docking #
  ###########
  def initial_CD(self, randomOnly=False):
    """
      Docks the ligand into the receptor

      Intermediate thermodynamic states are chosen such that
      thermodynamic length intervals are approximately constant.
      Configurations from each state are subsampled to seed the next simulation.
    """

    if (len(self.data['CD'].protocol) >
        0) and (self.data['CD'].protocol[-1]['crossed']):
      return  # Initial CD already complete

    from AlGDock.ligand_preparation import LigandPreparation
    seeds = LigandPreparation(self.args, self.log, self.top, self.system,
                              self._get_confs_to_rescore, self.iterator,
                              self.data).run('CD')

    from AlGDock.initialization import Initialization
    Initialization(self.args, self.log, self.top, self.system,
                  self.iterator, self.data, self.save, self._u_kln).run('CD', seeds)

    return True

  def calc_f_RL(self, readOnly=False, do_solvation=True, redo=False):
    """
    Calculates the binding potential of mean force
    redo recalculates f_RL and B except grid_MBAR
    """
    if self.data['CD'].protocol == []:
      return  # Initial CD is incomplete

    # Initialize variables as empty lists or by loading data
    if self.args.params['CD']['pose'] == -1:
      f_RL_FN = os.path.join(self.args.dir['CD'], 'f_RL.pkl.gz')
    else:
      f_RL_FN = os.path.join(self.args.dir['CD'], \
        'f_RL_pose%03d.pkl.gz'%self.args.params['CD']['pose'])

    dat = load_pkl_gz(f_RL_FN)
    if (dat is not None):
      (self.f_L, self.stats_RL, self.f_RL, self.B) = dat
    else:
      self._clear_f_RL()
    if readOnly:
      return True

    if redo:
      for key in self.f_RL.keys():
        if key != 'grid_MBAR':
          self.f_RL[key] = []
      self.B = {'MMTK_MBAR': []}
      for phase in self.args.params['CD']['phases']:
        for method in ['min_Psi', 'mean_Psi', 'EXP', 'MBAR']:
          self.B[phase + '_' + method] = []

    # Make sure all the energies are available
    for c in range(self.data['CD'].cycle):
      if len(self.data['CD'].Es[-1][c].keys()) == 0:
        self.log.tee("  skipping the binding PMF calculation")
        return
    if not hasattr(self, 'f_L'):
      self.log.tee("  skipping the binding PMF calculation")
      return

    start_string = "\n>>> Complex free energy calculations, starting at " + \
      time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "\n"
    self.log.recordStart('BPMF')

    updated = False

    def set_updated_to_True(updated, start_string, quiet=False):
      if (updated is False):
        self.log.set_lock('CD')
        if not quiet:
          self.log.tee(start_string)
      return True

    K = len(self.data['CD'].protocol)

    # Store stats_RL
    # Internal energies
    self.stats_RL['u_K_sampled'] = \
      [self._u_kln([self.data['CD'].Es[-1][c]],[self.data['CD'].protocol[-1]]) \
        for c in range(self.data['CD'].cycle)]
    self.stats_RL['u_KK'] = \
      [np.sum([self._u_kln([self.data['CD'].Es[k][c]],[self.data['CD'].protocol[k]]) \
        for k in range(len(self.data['CD'].protocol))],0) \
          for c in range(self.data['CD'].cycle)]

    # Interaction energies
    for c in range(len(self.stats_RL['Psi_grid']), self.data['CD'].cycle):
      self.stats_RL['Psi_grid'].append(
          (self.data['CD'].Es[-1][c]['LJr'] + \
           self.data['CD'].Es[-1][c]['LJa'] + \
           self.data['CD'].Es[-1][c]['ELE'])/(R*self.T_TARGET))
      updated = set_updated_to_True(updated,
                                    start_string,
                                    quiet=not do_solvation)

    # Estimate cycle at which simulation has equilibrated
    eqc_o = self.stats_RL['equilibrated_cycle']
    self.stats_RL['equilibrated_cycle'] = self._get_equilibrated_cycle('CD')
    if self.stats_RL['equilibrated_cycle'] != eqc_o:
      updated = set_updated_to_True(updated,
                                    start_string,
                                    quiet=not do_solvation)

    # Store rmsd values
    if (self.args.params['CD']['rmsd'] is not False):
      k = len(self.data['CD'].protocol) - 1
      for c in range(self.data['CD'].cycle):
        if not 'rmsd' in self.data['CD'].Es[k][c].keys():
          confs = [conf for conf in self.data['CD'].confs['samples'][k][c]]
          self.data['CD'].Es[k][c]['rmsd'] = self.get_rmsds(confs)
    self.stats_RL['rmsd'] = [(np.hstack([self.data['CD'].Es[k][c]['rmsd']
      if 'rmsd' in self.data['CD'].Es[k][c].keys() else [] \
        for c in range(self.stats_RL['equilibrated_cycle'][-1], \
                       self.data['CD'].cycle)])) \
          for k in range(len(self.data['CD'].protocol))]

    # Calculate CD free energies that have not already been calculated
    while len(self.f_RL['grid_MBAR']) < self.data['CD'].cycle:
      self.f_RL['grid_MBAR'].append([])
    while len(self.stats_RL['mean_acc']) < self.data['CD'].cycle:
      self.stats_RL['mean_acc'].append([])

    for c in range(self.data['CD'].cycle):
      # If solvation free energies are not being calculated,
      # only calculate the grid free energy for the current cycle
      if (not do_solvation) and c < (self.data['CD'].cycle - 1):
        continue
      if self.f_RL['grid_MBAR'][c] != []:
        continue

      fromCycle = self.stats_RL['equilibrated_cycle'][c]
      extractCycles = range(fromCycle, c + 1)

      # Extract relevant energies
      CD_Es = [Es[fromCycle:c+1] \
        for Es in self.data['CD'].Es]

      # Use MBAR for the grid scaling free energy estimate
      (u_kln, N_k) = self._u_kln(CD_Es, self.data['CD'].protocol)
      MBAR = self._run_MBAR(u_kln, N_k)[0]
      self.f_RL['grid_MBAR'][c] = MBAR
      updated = set_updated_to_True(updated,
                                    start_string,
                                    quiet=not do_solvation)

      self.log.tee("  calculated grid scaling free energy of %.2f RT "%(\
                  self.f_RL['grid_MBAR'][c][-1])+\
               "using cycles %d to %d"%(fromCycle, c))

      # Average acceptance probabilities
      mean_acc = np.zeros(K - 1)
      for k in range(0, K - 1):
        (u_kln, N_k) = self._u_kln(CD_Es[k:k + 2],
                                   self.data['CD'].protocol[k:k + 2])
        N = min(N_k)
        acc = np.exp(-u_kln[0, 1, :N] - u_kln[1, 0, :N] + u_kln[0, 0, :N] +
                     u_kln[1, 1, :N])
        mean_acc[k] = np.mean(np.minimum(acc, np.ones(acc.shape)))
      self.stats_RL['mean_acc'][c] = mean_acc

    if not do_solvation:
      if updated:
        if not self.log.run_type.startswith('timed'):
          self.log.tee(write_pkl_gz(f_RL_FN, \
            (self.f_L, self.stats_RL, self.f_RL, self.B)))
        self.log.clear_lock('CD')
      return True

    # Make sure postprocessing is complete
    pp_complete = self._postprocess()
    if not pp_complete:
      return False
    self.calc_f_L()

    # Make sure all the phase energies are available
    for c in range(self.data['CD'].cycle):
      for phase in self.args.params['CD']['phases']:
        for prefix in ['L', 'RL']:
          if not prefix + phase in self.data['CD'].Es[-1][c].keys():
            self.log.tee("  postprocessed energies for %s unavailable" % phase)
            return

    # Store stats_RL internal energies for phases
    for phase in self.args.params['CD']['phases']:
      self.stats_RL['u_K_'+phase] = \
        [self.data['CD'].Es[-1][c]['RL'+phase][:,-1]/(R*self.T_TARGET) \
          for c in range(self.data['CD'].cycle)]

    # Interaction energies
    for phase in self.args.params['CD']['phases']:
      if (not 'Psi_' + phase in self.stats_RL):
        self.stats_RL['Psi_' + phase] = []
      for c in range(len(self.stats_RL['Psi_' + phase]),
                     self.data['CD'].cycle):
        self.stats_RL['Psi_'+phase].append(
          (self.data['CD'].Es[-1][c]['RL'+phase][:,-1] - \
           self.data['CD'].Es[-1][c]['L'+phase][:,-1] - \
           self.args.original_Es[0][0]['R'+phase][:,-1])/(R*self.T_TARGET))

    # Predict native pose
    if self.args.params['CD']['pose'] == -1:
      (self.stats_RL['pose_inds'], self.stats_RL['scores']) = \
        self._get_pose_prediction()

    # BPMF assuming receptor and complex solvation cancel
    self.B['MMTK_MBAR'] = [-self.f_L['BC_MBAR'][-1][-1] + \
      self.f_RL['grid_MBAR'][c][-1] for c in range(len(self.f_RL['grid_MBAR']))]

    # BPMFs
    for phase in self.args.params['CD']['phases']:
      for key in [phase + '_solv']:
        if not key in self.f_RL:
          self.f_RL[key] = []
      for method in ['min_Psi', 'mean_Psi', 'EXP', 'MBAR']:
        if not phase + '_' + method in self.B:
          self.B[phase + '_' + method] = []

      # Receptor solvation
      f_R_solv = self.args.original_Es[0][0]['R' + phase][:, -1] / (
        R * self.T_TARGET)

      for c in range(len(self.B[phase + '_MBAR']), self.data['CD'].cycle):
        updated = set_updated_to_True(updated, start_string)
        extractCycles = range(self.stats_RL['equilibrated_cycle'][c], c + 1)

        # From the full grid to the fully bound complex in phase
        u_RL = np.concatenate([\
          self.data['CD'].Es[-1][c]['RL'+phase][:,-1]/(R*self.T_TARGET) \
          for c in extractCycles])
        u_sampled = np.concatenate([\
          self.stats_RL['u_K_sampled'][c] for c in extractCycles])

        du = u_RL - u_sampled
        min_du = min(du)
        weights = np.exp(-du + min_du)

        # Filter outliers
        if self.args.params['CD']['pose'] > -1:
          toKeep = du > (np.mean(du) - 3 * np.std(du))
          du = du[toKeep]
          weights[~toKeep] = 0.

        weights = weights / sum(weights)

        # Exponential average
        f_RL_solv = -np.log(np.exp(-du + min_du).mean()) + min_du - f_R_solv

        # Interaction energies
        Psi = np.concatenate([self.stats_RL['Psi_'+phase][c] \
          for c in extractCycles])
        min_Psi = min(Psi)
        max_Psi = max(Psi)

        # Complex solvation
        self.f_RL[phase + '_solv'].append(f_RL_solv)

        # Various BPMF estimates
        self.B[phase + '_min_Psi'].append(min_Psi)
        self.B[phase + '_mean_Psi'].append(np.sum(weights * Psi))
        self.B[phase+'_EXP'].append(\
          np.log(sum(weights*np.exp(Psi-max_Psi))) + max_Psi)

        self.B[phase+'_MBAR'].append(\
          - self.f_L[phase+'_solv'][-1] - self.f_L['BC_MBAR'][-1][-1] \
          + self.f_RL['grid_MBAR'][-1][-1] + f_RL_solv)

        self.log.tee("  calculated %s binding PMF of %.5g RT with cycles %d to %d"%(\
          phase, self.B[phase+'_MBAR'][-1], \
          self.stats_RL['equilibrated_cycle'][c], c))

    if updated:
      self.log.tee(
        write_pkl_gz(f_RL_FN, (self.f_L, self.stats_RL, self.f_RL, self.B)))
      self.log.tee("\nElapsed time for binding PMF estimation: " + \
        HMStime(self.log.timeSince('BPMF')))
    self.log.clear_lock('CD')

  def _store_infinite_f_RL(self):
    if self.args.params['CD']['pose'] == -1:
      f_RL_FN = os.path.join(self.args.dir['CD'], 'f_RL.pkl.gz')
    else:
      f_RL_FN = os.path.join(self.args.dir['CD'],\
        'f_RL_pose%03d.pkl.gz'%self.args.params['CD']['pose'])
    self.log.tee(write_pkl_gz(f_RL_FN, (self.f_L, [], np.inf, np.inf)))

  def _get_equilibrated_cycle(self, process):
    # Get previous results, if any
    if process == 'BC':
      if hasattr(self,'stats_L') and \
          ('equilibrated_cycle' in self.stats_L.keys()) and \
          self.stats_L['equilibrated_cycle']!=[]:
        equilibrated_cycle = self.stats_L['equilibrated_cycle']
      else:
        equilibrated_cycle = [0]
    elif process == 'CD':
      if hasattr(self,'stats_RL') and \
          ('equilibrated_cycle' in self.stats_RL.keys()) and \
          self.stats_RL['equilibrated_cycle']!=[]:
        equilibrated_cycle = self.stats_RL['equilibrated_cycle']
      else:
        equilibrated_cycle = [0]

    # Estimate equilibrated cycle
    for last_c in range(len(equilibrated_cycle), \
        self.data[process].cycle):
      correlation_times = [np.inf] + [\
        pymbar.timeseries.integratedAutocorrelationTime(\
          np.concatenate([self.data[process].Es[0][c]['mean_energies'] \
            for c in range(start_c,len(self.data[process].Es[0])) \
            if 'mean_energies' in self.data[process].Es[0][c].keys()])) \
               for start_c in range(1,last_c)]
      g = 2 * np.array(correlation_times) + 1
      nsamples_tot = [n for n in reversed(np.cumsum([len(self.data[process].Es[0][c]['MM']) \
        for c in reversed(range(last_c))]))]
      nsamples_ind = nsamples_tot / g
      equilibrated_cycle_last_c = max(np.argmax(nsamples_ind), 1)
      equilibrated_cycle.append(equilibrated_cycle_last_c)

    return equilibrated_cycle

  def _get_rmsd_matrix(self):
    process = 'CD'
    equilibrated_cycle = self.stats_RL['equilibrated_cycle'][-1]

    # Gather snapshots
    for k in range(equilibrated_cycle, self.data[process].cycle):
      if not isinstance(self.data[process].confs['samples'][-1][k], list):
        self.data[process].confs['samples'][-1][k] = [
          self.data[process].confs['samples'][-1][k]
        ]
    import itertools
    confs = np.array([conf for conf in itertools.chain.from_iterable(\
      [self.data[process].confs['samples'][-1][c] \
        for c in range(equilibrated_cycle,self.data[process].cycle)])])

    cum_Nk = np.cumsum([0] + [len(self.data['CD'].confs['samples'][-1][c]) \
      for c in range(self.data['CD'].cycle)])
    nsamples = cum_Nk[-1]

    # Obtain a full rmsd matrix
    # TODO: Check this
    if ('rmsd_matrix' in self.stats_RL.keys()) and \
        (len(self.stats_RL['rmsd_matrix'])==(nsamples*(nsamples-1)/2)):
      rmsd_matrix = stats_RL['rmsd_matrix']
    else:
      # Create a new matrix
      rmsd_matrix = []
      for c in range(len(confs)):
        rmsd_matrix.extend(self.get_rmsds(confs[c + 1:], confs[c]))
      rmsd_matrix = np.clip(rmsd_matrix, 0., None)
      self.stats_RL['rmsd_matrix'] = rmsd_matrix

    # TODO: Write code to extend previous matrix
    # Extend a previous matrix
    # rmsd_matrix = self.stats_RL['rmsd_matrix']
    # from scipy.spatial.distance import squareform
    # rmsd_matrix_sq = squareform(rmsd_matrix)
    #
    # for c in range(len(confs)):
    #   rmsd_matrix.extend(self.get_rmsds(confs[c+1:], confs[c]))
    # rmsd_matrix = np.clip(rmsd_matrix, 0., None)
    # self.stats_RL['rmsd_matrix'] = rmsd_matrix

    return rmsd_matrix

  def _cluster_samples(self, rmsd_matrix):
    # Clustering
    import scipy.cluster
    Z = scipy.cluster.hierarchy.linkage(rmsd_matrix, method='complete')
    assignments = np.array(\
      scipy.cluster.hierarchy.fcluster(Z, 0.1, criterion='distance'))

    # Reindexes the assignments in order of appearance
    new_index = 0
    mapping_to_new_index = {}
    for assignment in assignments:
      if not assignment in mapping_to_new_index.keys():
        mapping_to_new_index[assignment] = new_index
        new_index += 1
    assignments = [mapping_to_new_index[a] for a in assignments]
    return assignments

  def _get_pose_prediction(self, representative='medoid'):
    process = 'CD'
    equilibrated_cycle = self.stats_RL['equilibrated_cycle'][-1]
    stats = self.stats_RL

    rmsd_matrix = self._get_rmsd_matrix()
    assignments = self._cluster_samples(rmsd_matrix)

    cum_Nk = np.cumsum([0] + [len(self.data[process].confs['samples'][-1][c]) \
      for c in range(equilibrated_cycle,self.data[process].cycle)])

    def linear_index_to_pair(ind):
      cycle = list(ind < cum_Nk).index(True) - 1
      n = ind - cum_Nk[cycle]
      return (cycle + equilibrated_cycle, n)

    # Select a representative of each cluster
    pose_inds = []
    scores = {}

    if representative == 'medoid':
      # based on the medoid
      from scipy.spatial.distance import squareform
      rmsd_matrix_sq = squareform(rmsd_matrix)
      for n in range(max(assignments) + 1):
        inds = [i for i in range(len(assignments)) if assignments[i] == n]
        rmsd_matrix_n = rmsd_matrix_sq[inds][:, inds]
        (cycle,
         n) = linear_index_to_pair(inds[np.argmin(np.mean(rmsd_matrix_n, 0))])
        pose_inds.append((cycle, n))
    else:
      if 'Psi_' + representative in stats.keys():
        # based on the lowest interaction energy in specified phase
        phase = representative
        Psi_n = np.concatenate([stats['Psi_'+phase][c] \
                  for c in range(equilibrated_cycle,self.data[process].cycle)])
        for n in range(max(assignments) + 1):
          inds = [i for i in range(len(assignments)) if assignments[i] == n]
          (cycle, n) = linear_index_to_pair(inds[np.argmin(Psi_n[inds])])
          pose_inds.append((cycle, n))

    # If relevant, store the rmsd of the representatives
    if self.args.params['CD']['rmsd']:
      scores['rmsd'] = []
      for (cycle, n) in pose_inds:
        scores['rmsd'].append(self.data['CD'].Es[-1][cycle]['rmsd'][n])

    # Score clusters based on total energy
    uo = np.concatenate([stats['u_K_sampled'][c] \
      for c in range(equilibrated_cycle,self.data[process].cycle)])
    for phase in (['grid'] + self.args.params[process]['phases']):
      if phase != 'grid':
        un = np.concatenate([stats['u_K_'+phase][c] \
          for c in range(equilibrated_cycle,self.data[process].cycle)])
        du = un - uo
        min_du = min(du)
        weights = np.exp(-du + min_du)
      else:
        un = uo
        weights = np.ones(len(assignments))
      cluster_counts = np.histogram(assignments, \
        bins=np.arange(len(set(assignments))+1)-0.5,
        weights=weights)[0]
      # by free energy
      cluster_fe = -np.log(cluster_counts)
      cluster_fe -= np.min(cluster_fe)
      scores[phase + '_fe_u'] = cluster_fe
      # by minimum and mean energy
      scores[phase + '_min_u'] = []
      scores[phase + '_mean_u'] = []
      for n in range(max(assignments) + 1):
        un_n = [un[i] for i in range(len(assignments)) if assignments[i] == n]
        scores[phase + '_min_u'].append(np.min(un_n))
        scores[phase + '_mean_u'].append(np.mean(un_n))

    if process == 'CD':
      # Score clusters based on interaction energy
      Psi_o = np.concatenate([stats['Psi_grid'][c] \
        for c in range(equilibrated_cycle,self.data[process].cycle)])
      for phase in (['grid'] + self.args.params[process]['phases']):
        if phase != 'grid':
          Psi_n = np.concatenate([stats['Psi_'+phase][c] \
            for c in range(equilibrated_cycle,self.data[process].cycle)])
          dPsi = Psi_n - Psi_o
          min_dPsi = min(dPsi)
          weights = np.exp(-dPsi + min_dPsi)
        else:
          Psi_n = Psi_o
          weights = np.ones(len(assignments))
        cluster_counts = np.histogram(assignments, \
          bins=np.arange(len(set(assignments))+1)-0.5,
          weights=weights)[0]
        # by free energy
        cluster_fe = -np.log(cluster_counts)
        cluster_fe -= np.min(cluster_fe)
        scores[phase + '_fe_Psi'] = cluster_fe
        # by minimum and mean energy
        scores[phase + '_min_Psi'] = []
        scores[phase + '_mean_Psi'] = []
        for n in range(max(assignments) + 1):
          Psi_n_n = [
            Psi_n[i] for i in range(len(assignments)) if assignments[i] == n
          ]
          scores[phase + '_min_Psi'].append(np.min(Psi_n_n))
          scores[phase + '_mean_Psi'].append(np.mean(Psi_n_n))

    for key in scores.keys():
      scores[key] = np.array(scores[key])

    return (pose_inds, scores)

  def configuration_energies(self, minimize=False, max_confs=None):
    """
    Calculates the energy for configurations from self.args.FNs['score']
    """
    # Determine the name of the file
    prefix = 'xtal' if self.args.FNs['score']=='default' else \
      os.path.basename(self.args.FNs['score']).split('.')[0]
    if minimize:
      prefix = 'min_' + prefix
    energyFN = os.path.join(self.args.dir['CD'], prefix + '.pkl.gz')

    # Set the force field to fully interacting
    params_full = self.system.paramsFromAlpha(1.0, 'CD')
    self.system.setParams(params_full)

    # Load the configurations
    if os.path.isfile(energyFN):
      (confs, Es) = load_pkl_gz(energyFN)
    else:
      (confs, Es) = self._get_confs_to_rescore(site=False, \
        minimize=minimize, sort=False)

    self.log.set_lock('CD')
    self.log.tee("\n>>> Calculating energies for %d configurations, "%len(confs) + \
      "starting at " + \
      time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "\n")
    self.log.recordStart('configuration_energies')

    updated = False
    # Calculate MM and OBC energies
    if not 'MM' in Es.keys():
      Es = self.system.energyTerms(confs, Es)
      solvation_o = self.args.params['CD']['solvation']
      self.args.params['CD']['solvation'] = 'Full'
      if self.system.isForce('OBC'):
        del self._forceFields['OBC']
      self.system.clear_evaluators()
      self.system.setParams(params_full)
      Es = self.system.energyTerms(confs, Es)
      self.args.params['CD']['solvation'] = solvation_o
      updated = True

    # Direct electrostatic energy
    FN = os.path.join(os.path.dirname(self.args.FNs['grids']['ELE']),
                      'direct_ele.nc')
    if not 'direct_ELE' in Es.keys() and os.path.isfile(FN):
      key = 'direct_ELE'
      Es[key] = np.zeros(len(confs))
      from AlGDock.ForceFields.Grid.Interpolation import InterpolationForceField
      FF = InterpolationForceField(FN, \
        scaling_property='scaling_factor_electrostatic')
      self.top.universe.setForceField(FF)
      for c in range(len(confs)):
        self.top.universe.setConfiguration(
          Configuration(self.top.universe, confs[c]))
        Es[key][c] = self.top.universe.energy()
      updated = True

    # Calculate symmetry-corrected RMSD
    if not 'rmsd' in Es.keys() and (self.args.params['CD']['rmsd'] is
                                    not False):
      Es['rmsd'] = self.get_rmsds(confs)
      updated = True

    if updated:
      self.log.tee("\nElapsed time for ligand MM, OBC, and grid energies: " + \
        HMStime(self.log.timeSince('configuration_energies')), \
        process='CD')
    self.log.clear_lock('CD')

    # Reduce the number of conformations
    if max_confs is not None:
      confs = confs[:max_confs]

    # Implicit solvent energies
    self._load_programs(self.args.params['CD']['phases'])

    self.data['CD'].confs['starting_poses'] = None
    self._postprocess([('original', 0, 0, 'R')])

    for phase in self.args.params['CD']['phases']:
      if not 'R' + phase in Es.keys():
        Es['R' + phase] = self.args.params['CD']['receptor_' + phase]

    toClear = []
    for phase in self.args.params['CD']['phases']:
      for moiety in ['L', 'RL']:
        if not moiety + phase in Es.keys():
          outputname = os.path.join(self.args.dir['CD'],
                                    '%s.%s%s' % (prefix, moiety, phase))
          if phase.startswith('NAMD'):
            traj_FN = os.path.join(self.args.dir['CD'],
                                   '%s.%s.dcd' % (prefix, moiety))
            self._write_traj(traj_FN, confs, moiety)
          elif phase.startswith('sander'):
            traj_FN = os.path.join(self.args.dir['CD'],
                                   '%s.%s.mdcrd' % (prefix, moiety))
            self._write_traj(traj_FN, confs, moiety)
          elif phase.startswith('gbnsr6'):
            traj_FN = os.path.join(self.args.dir['CD'], \
              '%s.%s%s'%(prefix,moiety,phase),'in.crd')
          elif phase.startswith('OpenMM'):
            traj_FN = None
          elif phase in ['APBS_PBSA']:
            traj_FN = os.path.join(self.args.dir['CD'],
                                   '%s.%s.pqr' % (prefix, moiety))
          else:
            raise Exception('Unknown phase!')
          if not traj_FN in toClear:
            toClear.append(traj_FN)
          for program in ['NAMD', 'sander', 'gbnsr6', 'OpenMM', 'APBS']:
            if phase.startswith(program):
              # TODO: Mechanism to do partial calculation
              Es[moiety+phase] = getattr(self,'_%s_Energy'%program)(confs, \
                moiety, phase, traj_FN, outputname, debug=DEBUG)
              updated = True
              # Get any data added since the calculation started
              if os.path.isfile(energyFN):
                (confs_o, Es_o) = load_pkl_gz(energyFN)
                for key in Es_o.keys():
                  if key not in Es.keys():
                    Es[key] = Es_o[key]
              # Store the data
              self.log.tee(write_pkl_gz(energyFN, (confs, Es)))
              break
    for FN in toClear:
      if (FN is not None) and os.path.isfile(FN):
        os.remove(FN)

    for key in Es.keys():
      Es[key] = np.array(Es[key])
    self._combine_MM_and_solvent(Es)

    if updated:
      self.log.set_lock('CD')
      self.log.tee("\nElapsed time for energies: " + \
        HMStime(self.log.timeSince('configuration_energies')), \
        process='CD')
      self.log.clear_lock('CD')

      # Get any data added since the calculation started
      if os.path.isfile(energyFN):
        (confs_o, Es_o) = load_pkl_gz(energyFN)
        for key in Es_o.keys():
          if key not in Es.keys():
            Es[key] = Es_o[key]

      # Store the data
      self.log.tee(write_pkl_gz(energyFN, (confs, Es)))
    return (confs, Es)

  ######################
  # Internal Functions #
  ######################

  def _replica_exchange(self, process):
    """
    Performs a cycle of replica exchange
    """
    if not process in ['CD', 'BC']:
      raise Exception('Process must be CD or BC')


# GMC

    def gMC_initial_setup():
      """
      Initialize BAT converter object.
      Decide which internal coord to crossover. Here, only the soft torsions will be crossovered.
      Produce a list of replica (state) index pairs to be swaped. Only Neighbor pairs will be swaped.
      Assume that self.top.universe, self.top.molecule and K (number of states) exist
      as global variables when the function is called.
      """
      from AlGDock.rigid_bodies import identifier
      import itertools
      BAT_converter = identifier(self.top.universe, self.top.molecule)
      BAT = BAT_converter.BAT(extended=True)
      # this assumes that the torsional angles are stored in the tail of BAT
      softTorsionId = [
        i + len(BAT) - BAT_converter.ntorsions
        for i in BAT_converter._softTorsionInd
      ]
      torsions_to_crossover = []
      for i in range(1, len(softTorsionId)):
        combinations = itertools.combinations(softTorsionId, i)
        for c in combinations:
          torsions_to_crossover.append(list(c))
      #
      BAT_converter.BAT_to_crossover = torsions_to_crossover
      if len(BAT_converter.BAT_to_crossover) == 0:
        self.log.tee('  GMC No BAT to crossover')
      state_indices = range(K)
      state_indices_to_swap = zip( state_indices[0::2], state_indices[1::2] ) + \
                      zip( state_indices[1::2], state_indices[2::2] )
      #
      return BAT_converter, state_indices_to_swap

    #
    def do_gMC(nr_attempts, BAT_converter, state_indices_to_swap,
               torsion_threshold):
      """
      Assume self.top.universe, confs, paramss, state_inds, inv_state_inds exist as global variables
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
      if len(BAT_converter.BAT_to_crossover) == 0:
        return 0., 0.
      #
      from random import randrange
      # get reduced energies and BAT for all configurations in confs
      BATs = []
      energies = np.zeros(K, dtype=float)
      for c_ind in range(K):
        s_ind = state_inds[c_ind]
        self.top.universe.setConfiguration(
          Configuration(self.top.universe, confs[c_ind]))
        BATs.append(np.array(BAT_converter.BAT(extended=True), dtype=float))
        self.system.setParams(paramss[s_ind])
        reduced_e = self.top.universe.energy() / (R * paramss[s_ind]['T'])
        energies[c_ind] = reduced_e
      #
      nr_sets_of_torsions = len(BAT_converter.BAT_to_crossover)
      #
      attempt_count, acc_count = 0, 0
      sweep_count = 0
      while True:
        sweep_count += 1
        if (sweep_count * K) > (1000 * nr_attempts):
          self.log.tee(
            '  GMC Sweep too many times, but few attempted. Consider reducing torsion_threshold.'
          )
          return attempt_count, acc_count
        #
        for state_pair in state_indices_to_swap:
          conf_ind_k0 = inv_state_inds[state_pair[0]]
          conf_ind_k1 = inv_state_inds[state_pair[1]]
          # check if it should attempt for this pair of states
          ran_set_torsions = BAT_converter.BAT_to_crossover[randrange(
            nr_sets_of_torsions)]
          do_crossover = np.any(
            np.abs(BATs[conf_ind_k0][ran_set_torsions] -
                   BATs[conf_ind_k1][ran_set_torsions]) >= torsion_threshold)
          if do_crossover:
            attempt_count += 1
            # BAT and reduced energies before crossover
            BAT_k0_be = copy.deepcopy(BATs[conf_ind_k0])
            BAT_k1_be = copy.deepcopy(BATs[conf_ind_k1])
            e_k0_be = energies[conf_ind_k0]
            e_k1_be = energies[conf_ind_k1]
            # BAT after crossover
            BAT_k0_af = copy.deepcopy(BAT_k0_be)
            BAT_k1_af = copy.deepcopy(BAT_k1_be)
            for index in ran_set_torsions:
              tmp = BAT_k0_af[index]
              BAT_k0_af[index] = BAT_k1_af[index]
              BAT_k1_af[index] = tmp
            # Cartesian coord and reduced energies after crossover.
            BAT_converter.Cartesian(BAT_k0_af)
            self.system.setParams(paramss[state_pair[0]])
            e_k0_af = self.top.universe.energy() / (
              R * paramss[state_pair[0]]['T'])
            conf_k0_af = copy.deepcopy(self.top.universe.configuration().array)
            #
            BAT_converter.Cartesian(BAT_k1_af)
            self.system.setParams(paramss[state_pair[1]])
            e_k1_af = self.top.universe.energy() / (
              R * paramss[state_pair[1]]['T'])
            conf_k1_af = copy.deepcopy(self.top.universe.configuration().array)
            #
            de = (e_k0_be - e_k0_af) + (e_k1_be - e_k1_af)
            # update confs, energies, BATS
            if (de > 0) or (np.random.uniform() < np.exp(de)):
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

    #
    self.log.set_lock(process)

    confs = self.data[process].confs['replicas']
    paramss = self.data[process].protocol

    terms = ['MM']
    if process == 'BC':
      terms += ['OBC']
    elif process == 'CD':
      if self.args.params['CD']['pose'] > -1:
        # Pose BPMF
        terms += ['k_angular_ext', 'k_spatial_ext', 'k_angular_int']
      else:
        terms += ['site']
      terms += scalables

    # A list of pairs of replica indicies
    K = len(paramss)
    pairs_to_swap = []
    for interval in range(1, min(5, K)):
      lower_inds = []
      for lowest_index in range(interval):
        lower_inds += range(lowest_index, K - interval, interval)
      upper_inds = np.array(lower_inds) + interval
      pairs_to_swap += zip(lower_inds, upper_inds)

    from repX import attempt_swaps

    # Setting the force field will load grids
    # before multiple processes are spawned
    for k in range(K):
      self.system.setParams(paramss[k])

    # If it has not been set up, set up Smart Darting
    self.iterator.initializeSmartDartingConfigurations(
      self.data[process].confs['SmartDarting'], process, self.data)

    # storage[key][sweep_index][state_index] will contain data
    # from the replica exchange sweeps
    storage = {}
    for var in ['confs', 'state_inds', 'energies']:
      storage[var] = []

    self.log.recordStart('repX cycle')

    if self.args.cores > 1:
      # Multiprocessing setup
      m = multiprocessing.Manager()
      task_queue = m.Queue()
      done_queue = m.Queue()

    # GMC
    do_gMC = self.args.params[process]['GMC_attempts'] > 0
    if do_gMC:
      self.log.tee('  Using GMC for %s' % process)
      nr_gMC_attempts = K * self.args.params[process]['GMC_attempts']
      torsion_threshold = self.args.params[process]['GMC_tors_threshold']
      gMC_attempt_count = 0
      gMC_acc_count = 0
      time_gMC = 0.0
      BAT_converter, state_indices_to_swap = gMC_initial_setup()

    # MC move statistics
    acc = {}
    att = {}
    for move_type in ['ExternalMC', 'SmartDarting', 'Sampler']:
      acc[move_type] = np.zeros(K, dtype=int)
      att[move_type] = np.zeros(K, dtype=int)
      self.log.timings[move_type] = 0.
    self.log.timings['repX'] = 0.

    mean_energies = []

    # Do replica exchange
    state_inds = range(K)
    inv_state_inds = range(K)
    nsweeps = self.args.params[process]['sweeps_per_cycle']
    nsnaps = nsweeps / self.args.params[process]['snaps_per_cycle']
    for sweep in range(nsweeps):
      E = {}
      for term in terms:
        E[term] = np.zeros(K, dtype=float)
      # Sample within each state
      if self.args.cores > 1:
        for k in range(K):
          task_queue.put((confs[k], process, paramss[state_inds[k]], False, k))
        for p in range(self.args.cores):
          task_queue.put('STOP')
        processes = [multiprocessing.Process(target=self.iterator.iteration_worker, \
            args=(task_queue, done_queue)) for p in range(self.args.cores)]
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
        results = [self.iterator.iteration(confs[k], process, \
            paramss[state_inds[k]], False, k) for k in range(K)]

      # GMC
      if do_gMC:
        time_start_gMC = time.time()
        att_count, acc_count = do_gMC(nr_gMC_attempts, BAT_converter,
                                      state_indices_to_swap, torsion_threshold)
        gMC_attempt_count += att_count
        gMC_acc_count += acc_count
        time_gMC = +(time.time() - time_start_gMC)

      # Store energies
      for k in range(K):
        confs[k] = results[k]['confs']
      mean_energies.append(np.mean([results[k]['Etot'] for k in range(K)]))
      E = self.system.energyTerms(confs, E, process=process)

      # Store MC move statistics
      for k in range(K):
        for move_type in ['ExternalMC', 'SmartDarting', 'Sampler']:
          key = 'acc_' + move_type
          if key in results[k].keys():
            acc[move_type][state_inds[k]] += results[k][key]
            att[move_type][state_inds[k]] += results[k]['att_' + move_type]
            self.log.timings[move_type] += results[k]['time_' + move_type]

      # Calculate u_ij (i is the replica, and j is the configuration),
      #    a list of arrays
      (u_ij, N_k) = self._u_kln(E, [paramss[state_inds[c]] for c in range(K)])
      # Do the replica exchange
      repX_start_time = time.time()
      (state_inds, inv_state_inds) = \
        attempt_swaps(state_inds, inv_state_inds, u_ij, pairs_to_swap, \
          self.args.params[process]['attempts_per_sweep'])
      self.log.timings['repX'] += (time.time() - repX_start_time)

      # Store data in local variables
      if (sweep + 1) % self.args.params[process]['snaps_per_cycle'] == 0:
        if (process == 'CD') and (self.args.params['CD']['rmsd'] is not False):
          E['rmsd'] = self.get_rmsds(confs)
        storage['confs'].append(list(confs))
        storage['state_inds'].append(list(state_inds))
        storage['energies'].append(copy.deepcopy(E))

    # GMC
    if do_gMC:
      self.log.tee('  {0}/{1} crossover attempts ({2:.3g}) accepted in {3}'.format(\
        gMC_acc_count, gMC_attempt_count, \
        float(gMC_acc_count)/float(gMC_attempt_count) \
          if gMC_attempt_count > 0 else 0, \
        HMStime(time_gMC)))

    # Report
    self.log.tee("  completed cycle %d in %s"%(self.data[process].cycle, \
      HMStime(self.log.timeSince('repX cycle'))))
    MC_report = " "
    for move_type in ['ExternalMC', 'SmartDarting', 'Sampler']:
      total_acc = np.sum(acc[move_type])
      total_att = np.sum(att[move_type])
      if total_att > 0:
        MC_report += " %s %d/%d=%.2f (%.1f s);"%(move_type, \
          total_acc, total_att, float(total_acc)/total_att, \
          self.log.timings[move_type])
    MC_report += " repX t %.1f s" % self.log.timings['repX']
    self.log.tee(MC_report)

    # Adapt HamiltonianMonteCarlo parameters
    if self.args.params[process]['sampler'] == 'HMC':
      acc_rates = np.array(acc['Sampler'], dtype=np.float) / att['Sampler']
      for k in range(K):
        acc_rate = acc_rates[k]
        if acc_rate > 0.8:
          paramss[k]['delta_t'] += 0.125 * MMTK.Units.fs
          paramss[k]['steps_per_trial'] = min(paramss[k]['steps_per_trial']*2,\
            self.args.params[process]['steps_per_sweep'])
        elif acc_rate < 0.4:
          if paramss[k]['delta_t'] < 2.0 * MMTK.Units.fs:
            paramss[k]['steps_per_trial'] = max(
              int(paramss[k]['steps_per_trial'] / 2.), 1)
          paramss[k]['delta_t'] -= 0.25 * MMTK.Units.fs
          if acc_rate < 0.1:
            paramss[k]['delta_t'] -= 0.25 * MMTK.Units.fs
        if paramss[k]['delta_t'] < 0.1 * MMTK.Units.fs:
          paramss[k]['delta_t'] = 0.1 * MMTK.Units.fs

    # Get indicies for sorting by thermodynamic state, not replica
    inv_state_inds = np.zeros((nsnaps, K), dtype=int)
    for snap in range(nsnaps):
      state_inds = storage['state_inds'][snap]
      for k in range(K):
        inv_state_inds[snap][state_inds[k]] = k

    # Sort energies and conformations by thermodynamic state
    # and store in global variables
    #   self.data[process].Es and self.data[process].confs['samples']
    # and also local variables
    #   Es_repX and confs_repX
    if (process == 'CD') and (self.args.params['CD']['rmsd'] is not False):
      terms.append('rmsd')  # Make sure to save the rmsd
    Es_repX = []
    for k in range(K):
      E_k = {}
      E_k_repX = {}
      if k == 0:
        E_k['acc'] = acc
        E_k['att'] = att
        E_k['mean_energies'] = mean_energies
      for term in terms:
        E_term = np.array([storage['energies'][snap][term][\
          inv_state_inds[snap][k]] for snap in range(nsnaps)])
        E_k[term] = E_term
        E_k_repX[term] = E_term
      self.data[process].Es[k].append(E_k)
      Es_repX.append([E_k_repX])

    confs_repX = []
    for k in range(K):
      confs_k = [storage['confs'][snap][inv_state_inds[snap][k]] \
        for snap in range(nsnaps)]
      if self.args.params[process]['keep_intermediate'] or \
          ((process=='BC') and (k==0)) or (k==(K-1)):
        self.data[process].confs['samples'][k].append(confs_k)
      confs_repX.append(confs_k)

    # Store final conformation of each replica
    self.data[process].confs['replicas'] = \
      [np.copy(storage['confs'][-1][inv_state_inds[-1][k]]) \
       for k in range(K)]

    if self.args.params[process]['darts_per_sweep'] > 0:
      self.system.setParams(self.data[process].protocol[-1])
      new_confs = [np.copy(conf) \
        for conf in data[process].confs['samples'][k][-1]]
      self.iterator.addSmartDartingConfigurations(new_confs, process,
                                                  self.data)

    self.data[process].cycle += 1
    self.save(process)
    self.log.tee("")
    self.log.clear_lock(process)

    # The code below is only for sampling importance resampling
    if not self.args.params[process]['sampling_importance_resampling']:
      return

    # Calculate appropriate free energy
    if process == 'BC':
      self.calc_f_L(do_solvation=False)
      f_k = self.f_L['BC_MBAR'][-1]
    elif process == 'CD':
      self.calc_f_RL(do_solvation=False)
      f_k = self.f_RL['grid_MBAR'][-1]

    # Get weights for sampling importance resampling
    # MBAR weights for replica exchange configurations
    (u_kln, N_k) = self._u_kln(Es_repX, paramss)

    # This is a more direct way to get the weights
    from pymbar.utils import kln_to_kn
    u_kn = kln_to_kn(u_kln, N_k=N_k)

    from pymbar.utils import logsumexp
    log_denominator_n = logsumexp(f_k - u_kn.T, b=N_k, axis=1)
    logW = f_k - u_kn.T - log_denominator_n[:, np.newaxis]
    W_nl = np.exp(logW)
    for k in range(K):
      W_nl[:, k] = W_nl[:, k] / np.sum(W_nl[:, k])

    # This is for conversion to 2 indicies: state and snapshot
    cum_N_state = np.cumsum([0] + list(N_k))

    def linear_index_to_snapshot_index(ind):
      state_index = list(ind < cum_N_state).index(True) - 1
      nis_index = ind - cum_N_state[state_index]
      return (state_index, nis_index)

    # Selects new replica exchange snapshots
    self.data[process].confs['replicas'] = []
    for k in range(K):
      (s,n) = linear_index_to_snapshot_index(\
        np.random.choice(range(W_nl.shape[0]), size = 1, p = W_nl[:,k])[0])
      self.data[process].confs['replicas'].append(np.copy(confs_repX[s][n]))

  def sim_process(self, process):
    """
    Simulate and analyze a BC or CD process.

    As necessary, first conduct an initial BC or CD
    and then run a desired number of replica exchange cycles.
    """
    if (self.data[process].protocol==[]) or \
       (not self.data[process].protocol[-1]['crossed']):
      time_left = getattr(self, 'initial_' + process)()
      if not time_left:
        return False

    # Main loop for replica exchange
    if (self.args.params[process]['repX_cycles'] is not None) and \
       ((self.data[process].cycle < \
         self.args.params[process]['repX_cycles'])):

      # Load configurations to score from another program
      if (process=='CD') and (self.data['CD'].cycle==1) and \
         (self.args.params['CD']['pose'] == -1) and \
         (self.args.FNs['score'] is not None) and \
         (self.args.FNs['score']!='default'):
        self.log.set_lock('CD')
        self.log.tee("\n>>> Reinitializing replica exchange configurations")
        self.system.setParams(self.system.paramsFromAlpha(1.0, 'CD'))
        confs = self._get_confs_to_rescore(\
          nconfs=len(self.data['CD'].protocol), site=True, minimize=True)[0]
        self.log.clear_lock('CD')
        if len(confs) > 0:
          self.data['CD'].confs['replicas'] = confs

      self.log.tee("\n>>> Replica exchange for {0}, starting at {1}\n".format(\
        process, time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())), \
        process=process)
      self.log.recordStart(process + '_repX_start')
      start_cycle = self.data[process].cycle
      cycle_times = []
      while (self.data[process].cycle <
             self.args.params[process]['repX_cycles']):
        self._replica_exchange(process)
        cycle_times.append(self.log.timeSince('repX cycle'))
        if process == 'CD':
          self._insert_CD_state_between_low_acc()
        if not self.log.isTimeForTask(cycle_times):
          return False
      self.log.tee("Elapsed time for %d cycles of replica exchange: %s"%(\
         (self.data[process].cycle - start_cycle), \
          HMStime(self.log.timeSince(process+'_repX_start'))), \
          process=process)

    # If there are insufficient configurations,
    #   do additional replica exchange on the BC process
    if (process == 'BC'):
      E_MM = []
      for k in range(len(self.data['BC'].Es[0])):
        E_MM += list(self.data['BC'].Es[0][k]['MM'])
      while len(E_MM) < self.args.params['CD']['seeds_per_state']:
        self.log.tee(
          "More samples from high temperature ligand simulation needed",
          process='BC')
        self._replica_exchange('BC')
        cycle_times.append(self.log.timeSince('repX cycle'))
        if not self.log.isTimeForTask(cycle_times):
          return False
        E_MM = []
        for k in range(len(self.data['BC'].Es[0])):
          E_MM += list(self.data['BC'].Es[0][k]['MM'])

    # Clear evaluators to save memory
    self.system.clear_evaluators()

    return True  # The process has completed

  def _insert_CD_state(self, alpha, clear=True):
    """
    Inserts a new thermodynamic state into the CD protocol.
    Samples for previous cycles are added by sampling importance resampling.
    Clears grid_MBAR.
    """
    # Defines a new thermodynamic state based on the neighboring state
    neighbor_ind = [alpha < p['alpha']
                    for p in self.data['CD'].protocol].index(True) - 1
    params_n = self.system.paramsFromAlpha(
      alpha, params_o=self.data['CD'].protocol[neighbor_ind])

    # For sampling importance resampling,
    # prepare an augmented matrix for pymbar calculations
    # with a new thermodynamic state
    (u_kln_s, N_k) = self._u_kln(self.data['CD'].Es, self.data['CD'].protocol)
    (K, L, N) = u_kln_s.shape

    u_kln_n = self._u_kln(self.data['CD'].Es, [params_n])[0]
    L += 1
    N_k = np.append(N_k, [0])

    u_kln = np.zeros([K, L, N])
    u_kln[:, :-1, :] = u_kln_s
    for k in range(K):
      u_kln[k, -1, :] = u_kln_n[k, 0, :]

    # Determine SIR weights
    weights = self._run_MBAR(u_kln, N_k, augmented=True)[1][:, -1]
    weights = weights / sum(weights)

    # Resampling
    # Convert linear indices to 3 indicies: state, cycle, and snapshot
    cum_N_state = np.cumsum([0] + list(N_k))
    cum_N_cycle = [np.cumsum([0] + [self.data['CD'].Es[k][c]['MM'].shape[0] \
      for c in range(len(self.data['CD'].Es[k]))]) for k in range(len(self.data['CD'].Es))]

    def linear_index_to_snapshot_index(ind):
      state_index = list(ind < cum_N_state).index(True) - 1
      nis_index = ind - cum_N_state[state_index]
      cycle_index = list(nis_index < cum_N_cycle[state_index]).index(True) - 1
      nic_index = nis_index - cum_N_cycle[state_index][cycle_index]
      return (state_index, cycle_index, nic_index)

    def snapshot_index_to_linear_index(state_index, cycle_index, nic_index):
      return cum_N_state[state_index] + cum_N_cycle[state_index][
        cycle_index] + nic_index

    # Terms to copy
    if self.args.params['CD']['pose'] > -1:
      # Pose BPMF
      terms = ['MM',\
        'k_angular_ext','k_spatial_ext','k_angular_int'] + scalables
    else:
      # BPMF
      terms = ['MM', 'site'] + scalables

    CD_Es_s = []
    confs_s = []
    for c in range(len(self.data['CD'].Es[0])):
      CD_Es_c = dict([(term, []) for term in terms])
      confs_c = []
      for n_in_c in range(len(self.data['CD'].Es[-1][c]['MM'])):
        if (cum_N_cycle[-1][c] == 0):
          (snapshot_s,snapshot_c,snapshot_n) = linear_index_to_snapshot_index(\
           np.random.choice(range(len(weights)), size = 1, p = weights)[0])
        else:
          snapshot_c = np.inf
          while (snapshot_c > c):
            (snapshot_s,snapshot_c,snapshot_n) = linear_index_to_snapshot_index(\
             np.random.choice(range(len(weights)), size = 1, p = weights)[0])
        for term in terms:
          CD_Es_c[term].append(\
            np.copy(self.data['CD'].Es[snapshot_s][snapshot_c][term][snapshot_n]))
        if self.args.params['CD']['keep_intermediate']:
          # Has not been tested:
          confs_c.append(\
            np.copy(self.data['CD'].confs['samples'][snapshot_s][snapshot_c]))
      for term in terms:
        CD_Es_c[term] = np.array(CD_Es_c[term])
      CD_Es_s.append(CD_Es_c)
      confs_s.append(confs_c)

    # Insert resampled values
    self.data['CD'].protocol.insert(neighbor_ind + 1, params_n)
    self.data['CD'].Es.insert(neighbor_ind + 1, CD_Es_s)
    self.data['CD'].confs['samples'].insert(neighbor_ind + 1, confs_s)
    self.data['CD'].confs['replicas'].insert(neighbor_ind+1, \
      np.copy(self.data['CD'].confs['replicas'][neighbor_ind]))

    if clear:
      self._clear_f_RL()

  def _insert_CD_state_between_low_acc(self):
    # Insert thermodynamic states between those with low acceptance probabilities
    eq_c = self._get_equilibrated_cycle('CD')[-1]

    def calc_mean_acc(k):
      CD_Es = [Es[eq_c:self.data['CD'].cycle] for Es in self.data['CD'].Es]
      (u_kln,N_k) = self._u_kln(CD_Es[k:k+2],\
                                self.data['CD'].protocol[k:k+2])
      N = min(N_k)
      acc = np.exp(-u_kln[0, 1, :N] - u_kln[1, 0, :N] + u_kln[0, 0, :N] +
                   u_kln[1, 1, :N])
      return np.mean(np.minimum(acc, np.ones(acc.shape)))

    updated = False
    k = 0
    while k < len(self.data['CD'].protocol) - 1:
      mean_acc = calc_mean_acc(k)
      # print k, self.data['CD'].protocol[k]['alpha'], self.data['CD'].protocol[k+1]['alpha'], mean_acc
      while mean_acc < 0.4:
        if not updated:
          updated = True
          self.log.set_lock('CD')
        alpha_k = self.data['CD'].protocol[k]['alpha']
        alpha_kp = self.data['CD'].protocol[k + 1]['alpha']
        alpha_n = (alpha_k + alpha_kp) / 2.
        report = '  inserted state'
        report += ' between %.5g and %.5g at %.5g\n' % (alpha_k, alpha_kp, alpha_n)
        report += '  to improve acceptance rate from %.5g ' % mean_acc
        self._insert_CD_state(alpha_n, clear=False)
        mean_acc = calc_mean_acc(k)
        report += 'to %.5g' % mean_acc
        # print k, self.data['CD'].protocol[k]['alpha'], self.data['CD'].protocol[k+1]['alpha'], mean_acc
        self.log.tee(report)
      k += 1
    if updated:
      self._clear_f_RL()
      self.save('CD')
      self.log.tee("")
      self.log.clear_lock('CD')

  def _get_confs_to_rescore(self,
                            nconfs=None,
                            site=False,
                            minimize=True,
                            sort=True):
    """Returns configurations to rescore and their corresponding energies

    Parameters
    ----------
    nconfs : int or None
      Number of configurations to keep. If it is smaller than the number
      of unique configurations, then the lowest energy configurations will
      be kept. If it is larger, then the lowest energy configuration will be
      duplicated. If it is None, then all unique configurations will be kept.
    site : bool
      If True, configurations that are outside of the binding site
      will be discarded.
    minimize : bool
      If True, the configurations will be minimized
    sort : bool
      If True, configurations and energies will be sorted by DECREASING energy.

    Returns
    -------
    confs : list of np.array
      Configurations
    energies : list of float
      Energies of the configurations
    """
    # Get configurations
    count = {'xtal': 0, 'dock6': 0, 'initial_CD': 0, 'duplicated': 0}

    # based on the score option
    if self.args.FNs['score'] == 'default':
      confs = [np.copy(self.data['CD'].confs['ligand'])]
      count['xtal'] = 1
      Es = {}
      if nconfs is None:
        nconfs = 1
    elif (self.args.FNs['score'] is None) or (not os.path.isfile(
        self.args.FNs['score'])):
      confs = []
      Es = {}
    elif self.args.FNs['score'].endswith('.mol2') or \
         self.args.FNs['score'].endswith('.mol2.gz'):
      import AlGDock.IO
      IO_dock6_mol2 = AlGDock.IO.dock6_mol2()
      (confs, Es) = IO_dock6_mol2.read(self.args.FNs['score'], \
        reorder=self.top.inv_prmtop_atom_order_L,
        multiplier=0.1) # to convert Angstroms to nanometers
      count['dock6'] = len(confs)
    elif self.args.FNs['score'].endswith('.mdcrd'):
      import AlGDock.IO
      IO_crd = AlGDock.IO.crd()
      lig_crds = IO_crd.read(self.args.FNs['score'], \
        multiplier=0.1) # to convert Angstroms to nanometers
      confs = np.array_split(
        lig_crds, lig_crds.shape[0] / self.top.universe.numberOfAtoms())
      confs = [conf[self.top.inv_prmtop_atom_order_L, :] for conf in confs]
      Es = {}
    elif self.args.FNs['score'].endswith('.nc'):
      from netCDF4 import Dataset
      dock6_nc = Dataset(self.args.FNs['score'], 'r')
      confs = [
        dock6_nc.variables['confs'][n][self.top.inv_prmtop_atom_order_L, :]
        for n in range(dock6_nc.variables['confs'].shape[0])
      ]
      Es = dict([(key, dock6_nc.variables[key][:])
                 for key in dock6_nc.variables.keys() if key != 'confs'])
      dock6_nc.close()
      count['dock6'] = len(confs)
    elif self.args.FNs['score'].endswith('.pkl.gz'):
      F = gzip.open(self.args.FNs['score'], 'r')
      confs = pickle.load(F)
      F.close()
      if not isinstance(confs, list):
        confs = [confs]
      Es = {}
    else:
      raise Exception('Input configuration format not recognized')

    # based on the seeds
    # TODO: Use CD seeds for BC
    if (self.data['CD'].confs['seeds'] is not None) and \
       (self.args.params['CD']['pose']==-1):
      confs = confs + self.data['CD'].confs['seeds']
      Es = {}
      count['initial_CD'] = len(self.data['CD'].confs['seeds'])

    if len(confs) == 0:
      return ([], {})

    if site:
      # Filters out configurations not in the binding site
      confs_in_site = []
      Es_in_site = dict([(label, []) for label in Es.keys()])
      old_eval = None
      if (None, None, None) in self.top.universe._evaluator.keys():
        old_eval = self.top.universe._evaluator[(None, None, None)]
      self.system.setParams({'site': True, 'T': self.T_TARGET})
      for n in range(len(confs)):
        self.top.universe.setConfiguration(
          Configuration(self.top.universe, confs[n]))
        if self.top.universe.energy() < 1.:
          confs_in_site.append(confs[n])
          for label in Es.keys():
            Es_in_site[label].append(Es[label][n])
      if old_eval is not None:
        self.top.universe._evaluator[(None, None, None)] = old_eval
      confs = confs_in_site
      Es = Es_in_site

    try:
      self.top.universe.energy()
    except ValueError:
      return (confs, {})

    if minimize:
      Es = {}
      (confs, energies) = self._checkedMinimizer(confs)
    else:
      # Evaluate energies
      energies = []
      for conf in confs:
        self.top.universe.setConfiguration(
          Configuration(self.top.universe, conf))
        energies.append(self.top.universe.energy())

    if sort and len(confs) > 0:
      # Sort configurations by DECREASING energy
      energies, confs = (list(l) for l in zip(*sorted(zip(energies, confs), \
        key=lambda p:p[0], reverse=True)))

    # Shrink or extend configuration and energy array
    if nconfs is not None:
      confs = confs[-nconfs:]
      energies = energies[-nconfs:]
      while len(confs) < nconfs:
        confs.append(confs[-1])
        energies.append(energies[-1])
        count['duplicated'] += 1
      count['nconfs'] = nconfs
    else:
      count['nconfs'] = len(confs)
    count['minimized'] = {True: ' minimized', False: ''}[minimize]
    Es['total'] = np.array(energies)

    self.log.tee(
      "  keeping {nconfs}{minimized} configurations out of\n  {xtal} from xtal, {dock6} from dock6, {initial_CD} from initial CD, and {duplicated} duplicated"
      .format(**count))
    return (confs, Es)

  def _checkedMinimizer(self, confs):
    """Minimizes configurations while checking for crashes and overflows

    Parameters
    ----------
    confs : list of np.array
      Configurations to minimize

    Returns
    -------
    confs : list of np.array
      Minimized configurations
    energies : list of float
      Energies of the minimized configurations
    """
    from MMTK.Minimization import SteepestDescentMinimizer  # @UnresolvedImport
    minimizer = SteepestDescentMinimizer(self.top.universe)

    original_stderr = sys.stderr
    sys.stderr = NullDevice()  # Suppresses warnings for minimization

    minimized_confs = []
    minimized_energies = []
    self.log.recordStart('minimization')
    for conf in confs:
      self.top.universe.setConfiguration(
        Configuration(self.top.universe, conf))
      x_o = np.copy(self.top.universe.configuration().array)
      e_o = self.top.universe.energy()
      for rep in range(50):
        minimizer(steps=25)
        x_n = np.copy(self.top.universe.configuration().array)
        e_n = self.top.universe.energy()
        diff = abs(e_o - e_n)
        if np.isnan(e_n) or diff < 0.05 or diff > 1000.:
          self.top.universe.setConfiguration(
            Configuration(self.top.universe, x_o))
          break
        else:
          x_o = x_n
          e_o = e_n
      if not np.isnan(e_o):
        minimized_confs.append(x_o)
        minimized_energies.append(e_o)

    sys.stderr = original_stderr  # Restores error reporting

    confs = minimized_confs
    energies = minimized_energies
    self.log.tee("  minimized %d configurations in "%len(confs) + \
      HMStime(self.log.timeSince('minimization')) + \
      "\n  the first %d energies are:\n  "%min(len(confs),10) + \
      ', '.join(['%.2f'%e for e in energies[:10]]))
    return confs, energies

  def _run_MBAR(self, u_kln, N_k, augmented=False):
    """
    Estimates the free energy of a transition using BAR and MBAR
    """
    import pymbar
    K = len(N_k) - 1 if augmented else len(N_k)
    f_k_FEPF = np.zeros(K)
    f_k_BAR = np.zeros(K)
    W_nl = None
    for k in range(K - 1):
      w_F = u_kln[k, k + 1, :N_k[k]] - u_kln[k, k, :N_k[k]]
      min_w_F = min(w_F)
      w_R = u_kln[k + 1, k, :N_k[k + 1]] - u_kln[k + 1, k + 1, :N_k[k + 1]]
      min_w_R = min(w_R)
      f_k_FEPF[k + 1] = -np.log(np.mean(np.exp(-w_F + min_w_F))) + min_w_F
      try:
        f_k_BAR[k+1] = pymbar.BAR(w_F, w_R, \
                       relative_tolerance=1.0E-5, \
                       verbose=False, \
                       compute_uncertainty=False)
      except:
        f_k_BAR[k + 1] = f_k_FEPF[k + 1]
        print 'Error with BAR. Using FEP.'
    f_k_FEPF = np.cumsum(f_k_FEPF)
    f_k_BAR = np.cumsum(f_k_BAR)
    try:
      if augmented:
        f_k_BAR = np.append(f_k_BAR, [0])
      f_k_pyMBAR = pymbar.MBAR(u_kln, N_k, \
        relative_tolerance=1.0E-5, \
        verbose = False, \
        initial_f_k = f_k_BAR, \
        maximum_iterations = 20)
      f_k_MBAR = f_k_pyMBAR.f_k
      W_nl = f_k_pyMBAR.getWeights()
    except:
      print N_k, f_k_BAR
      f_k_MBAR = f_k_BAR
      print 'Error with MBAR. Using BAR.'
    if np.isnan(f_k_MBAR).any():
      f_k_MBAR = f_k_BAR
      print 'Error with MBAR. Using BAR.'
    return (f_k_MBAR, W_nl)

  def _u_kln(self, eTs, paramss, noBeta=False):
    """
    Computes a reduced potential energy matrix.  k is the sampled state.  l is the state for which energies are evaluated.

    Input:
    eT is a
      -dictionary (of mapped energy terms) of numpy arrays (over states)
      -list (over states) of dictionaries (of mapped energy terms) of numpy arrays (over configurations), or a
      -list (over states) of lists (over cycles) of dictionaries (of mapped energy terms) of numpy arrays (over configurations)
    paramss is a list of thermodynamic states
    noBeta means that the energy will not be divided by RT

    Output: u_kln or (u_kln, N_k)
    u_kln is the matrix (as a numpy array)
    N_k is an array of sample sizes
    """
    L = len(paramss)

    addMM = ('MM' in paramss[0].keys()) and (paramss[0]['MM'])
    addSite = ('site' in paramss[0].keys()) and (paramss[0]['site'])
    probe_keys = ['MM','k_angular_ext','k_spatial_ext','k_angular_int'] + \
      scalables
    probe_key = [key for key in paramss[0].keys() if key in probe_keys][0]

    if isinstance(eTs, dict):
      # There is one configuration per state
      K = len(eTs[probe_key])
      N_k = np.ones(K, dtype=int)
      u_kln = []
      E_base = np.zeros(K)
      if addMM:
        E_base += eTs['MM']
      if addSite:
        E_base += eTs['site']
      for l in range(L):
        E = 1. * E_base
        for scalable in scalables:
          if scalable in paramss[l].keys():
            E += paramss[l][scalable] * eTs[scalable]
        for key in ['k_angular_ext', 'k_spatial_ext', 'k_angular_int']:
          if key in paramss[l].keys():
            E += paramss[l][key] * eTs[key]
        if noBeta:
          u_kln.append(E)
        else:
          u_kln.append(E / (R * paramss[l]['T']))
    elif isinstance(eTs[0], dict):
      K = len(eTs)
      N_k = np.array([len(eTs[k][probe_key]) for k in range(K)])
      u_kln = np.zeros([K, L, N_k.max()], np.float)

      for k in range(K):
        E_base = 0.0
        if addMM:
          E_base += eTs[k]['MM']
        if addSite:
          E_base += eTs[k]['site']
        for l in range(L):
          E = 1. * E_base
          for scalable in scalables:
            if scalable in paramss[l].keys():
              E += paramss[l][scalable] * eTs[k][scalable]
          for key in ['k_angular_ext', 'k_spatial_ext', 'k_angular_int']:
            if key in paramss[l].keys():
              E += paramss[l][key] * eTs[k][key]
          if noBeta:
            u_kln[k, l, :N_k[k]] = E
          else:
            u_kln[k, l, :N_k[k]] = E / (R * paramss[l]['T'])
    elif isinstance(eTs[0], list):
      K = len(eTs)
      N_k = np.zeros(K, dtype=int)

      for k in range(K):
        for c in range(len(eTs[k])):
          N_k[k] += len(eTs[k][c][probe_key])
      u_kln = np.zeros([K, L, N_k.max()], np.float)

      for k in range(K):
        E_base = 0.0
        C = len(eTs[k])
        if addMM:
          E_base += np.concatenate([eTs[k][c]['MM'] for c in range(C)])
        if addSite:
          E_base += np.concatenate([eTs[k][c]['site'] for c in range(C)])
        for l in range(L):
          E = 1. * E_base
          for scalable in scalables:
            if scalable in paramss[l].keys():
              E += paramss[l][scalable]*np.concatenate([eTs[k][c][scalable] \
                for c in range(C)])
          for key in ['k_angular_ext', 'k_spatial_ext', 'k_angular_int']:
            if key in paramss[l].keys():
              E += paramss[l][key]*np.concatenate([eTs[k][c][key] \
                for c in range(C)])
          if noBeta:
            u_kln[k, l, :N_k[k]] = E
          else:
            u_kln[k, l, :N_k[k]] = E / (R * paramss[l]['T'])

    if (K == 1) and (L == 1):
      return u_kln.ravel()
    else:
      return (u_kln, N_k)

  def _load_programs(self, phases):
    # Find the necessary programs, downloading them if necessary
    programs = []
    for phase in phases:
      for (prefix,program) in [('NAMD','namd'), \
          ('sander','sander'), ('gbnsr6','gbnsr6'), ('APBS','apbs')]:
        if phase.startswith(prefix) and not program in programs:
          programs.append(program)
      if phase.find('ALPB') > -1:
        if not 'elsize' in programs:
          programs.append('elsize')
        if not 'ambpdb' in programs:
          programs.append('ambpdb')
    if 'apbs' in programs:
      for program in ['ambpdb', 'molsurf']:
        if not program in programs:
          programs.append(program)
    for program in programs:
      self.args.FNs[program] = path_tools.findPaths([program])[program]
    path_tools.loadModules(programs)

  def _postprocess(self,
      conditions=[('original',0, 0,'R'), ('BC',-1,-1,'L'), \
                  ('CD',   -1,-1,'L'), ('CD',-1,-1,'RL')],
      phases=None,
      readOnly=False, redo_CD=False, debug=DEBUG):
    """
    Obtains the NAMD energies of all the conditions using all the phases.
    Saves both MMTK and NAMD energies after NAMD energies are estimated.

    state == -1 means the last state
    cycle == -1 means all cycles

    """
    # Clear evaluators to save memory
    self.system.clear_evaluators()

    if phases is None:
      phases = list(set(self.args.params['BC']['phases'] + \
        self.args.params['CD']['phases']))

    updated_processes = []

    # Identify incomplete calculations
    incomplete = []
    for (p, state, cycle, moiety) in conditions:
      # Check that the values are legitimate
      if not p in ['BC', 'CD', 'original']:
        raise Exception("Type should be in ['BC', 'CD', 'original']")
      if not moiety in ['R', 'L', 'RL']:
        raise Exception("Species should in ['R','L', 'RL']")
      if p != 'original' and self.data[p].protocol == []:
        continue
      if state == -1:
        state = len(self.data[p].protocol) - 1
      if cycle == -1:
        cycles = range(self.data[p].cycle)
      else:
        cycles = [cycle]

      # Check for completeness
      for c in cycles:
        for phase in phases:
          label = moiety + phase

          # Skip postprocessing
          # if the function is NOT being rerun in redo_CD mode
          # and one of the following:
          # the function is being run in readOnly mode,
          # the energies are already in memory.
          if (not (redo_CD and p=='CD')) and \
            (readOnly \
            or (p == 'original' and \
                (label in self.args.original_Es[state][c].keys()) and \
                (self.args.original_Es[state][c][label] is not None)) \
            or (p != 'original' and \
                ('MM' in self.data[p].Es[state][c].keys()) and \
                (label in self.data[p].Es[state][c].keys()) and \
                (len(self.data[p].Es[state][c]['MM'])==\
                 len(self.data[p].Es[state][c][label])))):
            pass
          else:
            incomplete.append((p, state, c, moiety, phase))

    if incomplete == []:
      return True

    del p, state, c, moiety, phase, cycles, label

    self._load_programs([val[-1] for val in incomplete])

    # Write trajectories and queue calculations
    m = multiprocessing.Manager()
    task_queue = m.Queue()
    time_per_snap = m.dict()
    for (p, state, c, moiety, phase) in incomplete:
      if moiety + phase not in time_per_snap.keys():
        time_per_snap[moiety + phase] = m.list()

    # Decompress prmtop and inpcrd files
    decompress = (self.args.FNs['prmtop'][moiety].endswith('.gz')) or \
                 (self.args.FNs['inpcrd'][moiety].endswith('.gz'))
    if decompress:
      for key in ['prmtop', 'inpcrd']:
        if self.args.FNs[key][moiety].endswith('.gz'):
          import shutil
          shutil.copy(self.args.FNs[key][moiety],
                      self.args.FNs[key][moiety] + '.BAK')
          os.system('gunzip -f ' + self.args.FNs[key][moiety])
          os.rename(self.args.FNs[key][moiety] + '.BAK',
                    self.args.FNs[key][moiety])
          self.args.FNs[key][moiety] = self.args.FNs[key][moiety][:-3]

    toClean = []

    for (p, state, c, moiety, phase) in incomplete:
      # Identify the configurations
      if (moiety == 'R'):
        if not 'receptor' in self.data['CD'].confs.keys():
          continue
        confs = [self.data['CD'].confs['receptor']]
      else:
        confs = self.data[p].confs['samples'][state][c]

      # Identify the file names
      if p == 'original':
        prefix = p
      else:
        prefix = '%s%d_%d' % (p, state, c)

      p_dir = {
        'BC': self.args.dir['BC'],
        'original': self.args.dir['CD'],
        'CD': self.args.dir['CD']
      }[p]

      if phase.startswith('NAMD'):
        traj_FN = os.path.join(p_dir, '%s.%s.dcd' % (prefix, moiety))
      elif phase.startswith('sander'):
        traj_FN = os.path.join(p_dir, '%s.%s.mdcrd' % (prefix, moiety))
      elif phase.startswith('gbnsr6'):
        traj_FN = os.path.join(p_dir, '%s.%s%s' % (prefix, moiety, phase),
                               'in.crd')
      elif phase.startswith('OpenMM'):
        traj_FN = None
      elif phase in ['APBS_PBSA']:
        traj_FN = os.path.join(p_dir, '%s.%s.pqr' % (prefix, moiety))
      outputname = os.path.join(p_dir, '%s.%s%s' % (prefix, moiety, phase))

      # Writes trajectory
      self._write_traj(traj_FN, confs, moiety)
      if (traj_FN is not None) and (not traj_FN in toClean):
        toClean.append(traj_FN)

      # Queues the calculations
      task_queue.put((confs, moiety, phase, traj_FN, outputname, debug, \
              (p,state,c,moiety+phase)))

    # Start postprocessing
    # self.log.set_lock('CD' if 'CD' in [loc[0] for loc in incomplete] else 'BC')
    self.log.tee("\n>>> Postprocessing, starting at " + \
      time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "\n")
    self.log.recordStart('postprocess')

    done_queue = m.Queue()
    processes = [multiprocessing.Process(target=self._energy_worker, \
        args=(task_queue, done_queue, time_per_snap)) \
        for p in range(self.args.cores)]
    for process in range(self.args.cores):
      task_queue.put('STOP')
    for process in processes:
      process.start()
    for process in processes:
      process.join()
    results = []
    while not done_queue.empty():
      results.append(done_queue.get())
    for process in processes:
      process.terminate()

    # Clean up files
    if not debug:
      for FN in toClean:
        if os.path.isfile(FN):
          os.remove(FN)

    # Clear decompressed files
    if decompress:
      for key in ['prmtop', 'inpcrd']:
        if os.path.isfile(self.args.FNs[key][moiety] + '.gz'):
          os.remove(self.args.FNs[key][moiety])
          self.args.FNs[key][moiety] = self.args.FNs[key][moiety] + '.gz'

    # Store energies
    updated_energy_dicts = []
    for (E, (p, state, c, label), wall_time) in results:
      if p == 'original':
        self.args.original_Es[state][c][label] = E
        updated_energy_dicts.append(self.args.original_Es[state][c])
      else:
        self.data[p].Es[state][c][label] = E
        updated_energy_dicts.append(self.data[p].Es[state][c])
      if not p in updated_processes:
        updated_processes.append(p)
    for d in updated_energy_dicts:
      self._combine_MM_and_solvent(d)

    # Print time per snapshot
    for key in time_per_snap.keys():
      if len(time_per_snap[key]) > 0:
        mean_time_per_snap = np.mean(time_per_snap[key])
        if not np.isnan(mean_time_per_snap):
          self.log.tee("  an average of %.5g s per %s snapshot"%(\
            mean_time_per_snap, key))
        else:
          self.log.tee("  time per snapshot in %s: "%(key) + \
            ', '.join(['%.5g'%t for t in time_per_snap[key]]))
      else:
        self.log.tee("  no snapshots postprocessed in %s" % (key))

    # Save data
    if 'original' in updated_processes:
      for phase in phases:
        if (self.args.params['CD']['receptor_'+phase] is None) and \
           (self.args.original_Es[0][0]['R'+phase] is not None):
          self.args.params['CD']['receptor_'+phase] = \
            self.args.original_Es[0][0]['R'+phase]
      self.save('CD', keys=['progress'])
    if 'BC' in updated_processes:
      self.save('BC')
    if ('CD' in updated_processes) or ('original' in updated_processes):
      self.save('CD')

    if len(updated_processes) > 0:
      self.log.tee("\nElapsed time for postprocessing: " + \
        HMStime(self.log.timeSince('postprocess')))
      self.log.clear_lock('CD' if 'CD' in updated_processes else 'BC')
      return len(incomplete) == len(results)

  def _energy_worker(self, input, output, time_per_snap):
    for args in iter(input.get, 'STOP'):
      (confs, moiety, phase, traj_FN, outputname, debug, reference) = args
      (p, state, c, label) = reference
      nsnaps = len(confs)

      # Make sure there is enough time remaining
      if len(time_per_snap[moiety + phase]) > 0:
        if not self.log.isTimeForTask(nsnaps * np.array(time_per_snap[moiety + phase])):
          return

      # Calculate the energy
      self.log.recordStart('energy')
      for program in ['NAMD', 'sander', 'gbnsr6', 'OpenMM', 'APBS']:
        if phase.startswith(program):
          E = np.array(getattr(self, '_%s_Energy' % program)(*args))
          break
      wall_time = self.log.timeSince('energy')

      if not np.isinf(E).any():
        self.log.tee("  postprocessed %s, state %d, cycle %d, %s in %s"%(\
          p,state,c,label,HMStime(wall_time)))

        # Store output and timings
        output.put((E, reference, wall_time))

        time_per_snap_list = time_per_snap[moiety + phase]
        time_per_snap_list.append(wall_time / nsnaps)
        time_per_snap[moiety + phase] = time_per_snap_list
      else:
        self.log.tee("  error in postprocessing %s, state %d, cycle %d, %s in %s"%(\
          p,state,c,label,HMStime(wall_time)))
        return

  def _NAMD_Energy(self,
                   confs,
                   moiety,
                   phase,
                   dcd_FN,
                   outputname,
                   debug=DEBUG,
                   reference=None):
    """
    Uses NAMD to calculate the energy of a set of configurations
    Units are the MMTK standard, kJ/mol
    """
    # NAMD ENERGY FIELDS:
    # 0. TS 1. BOND 2. ANGLE 3. DIHED 4. IMPRP 5. ELECT 6. VDW 7. BOUNDARY
    # 8. MISC 9. KINETIC 10. TOTAL 11. TEMP 12. POTENTIAL 13. TOTAL3 14. TEMPAVG
    # The saved fields are energyFields=[1, 2, 3, 4, 5, 6, 8, 12],
    # and thus the new indicies are
    # 0. BOND 1. ANGLE 2. DIHED 3. IMPRP 4. ELECT 5. VDW 6. MISC 7. POTENTIAL

    # Run NAMD
    import AlGDock.NAMD
    energyCalc = AlGDock.NAMD.NAMD(\
      prmtop=self.args.FNs['prmtop'][moiety], \
      inpcrd=self.args.FNs['inpcrd'][moiety], \
      fixed={'R':self.args.FNs['fixed_atoms']['R'], \
             'L':None, \
             'RL':self.args.FNs['fixed_atoms']['RL']}[moiety], \
      solvent={'NAMD_OBC':'GBSA', 'NAMD_Gas':'Gas'}[phase], \
      useCutoff=(phase=='NAMD_OBC'), \
      namd_command=self.args.FNs['namd'])
    E = energyCalc.energies_PE(\
      outputname, dcd_FN, energyFields=[1, 2, 3, 4, 5, 6, 8, 12], \
      keepScript=debug, write_energy_pkl_gz=False)

    return np.array(E, dtype=float) * MMTK.Units.kcal / MMTK.Units.mol

  def _sander_Energy(self, confs, moiety, phase, AMBER_mdcrd_FN, \
      outputname=None, debug=DEBUG, reference=None):
    self.args.dir['out'] = os.path.dirname(os.path.abspath(AMBER_mdcrd_FN))
    script_FN = '%s%s.in' % ('.'.join(AMBER_mdcrd_FN.split('.')[:-1]), phase)
    out_FN = '%s%s.out' % ('.'.join(AMBER_mdcrd_FN.split('.')[:-1]), phase)

    script_F = open(script_FN, 'w')
    script_F.write('''Calculating energies with sander
&cntrl
  imin=5,    ! read trajectory in for analysis
  ntx=1,     ! input is read formatted with no velocities
  irest=0,
  ntb=0,     ! no periodicity and no PME
  idecomp=0, ! no decomposition
  ntc=1,     ! No SHAKE
  cut=9999., !''')
    if phase == 'sander_Gas':
      script_F.write("""
  ntf=1,     ! Complete interaction is calculated
/
""")
    elif phase == 'sander_PBSA':
      fillratio = 4.0 if moiety == 'L' else 2.0
      script_F.write('''
  ntf=7,     ! No bond, angle, or dihedral forces calculated
  ipb=2,     ! Default PB dielectric model
  inp=2,     ! non-polar from cavity + dispersion
/
&pb
  radiopt=0, ! Use atomic radii from the prmtop file
  fillratio=%d,
  sprob=1.4,
  cavity_surften=0.0378, ! (kcal/mol) Default in MMPBSA.py
  cavity_offset=-0.5692, ! (kcal/mol) Default in MMPBSA.py
/
''' % fillratio)
    else:
      if phase.find('ALPB') > -1 and moiety.find('R') > -1:
        script_F.write("\n  alpb=1,")
        script_F.write("\n  arad=%.2f," % self.elsize)
      key = phase.split('_')[-1]
      igb = {'HCT': 1, 'OBC1': 2, 'OBC2': 5, 'GBn': 7, 'GBn2': 8}[key]
      script_F.write('''
  ntf=7,     ! No bond, angle, or dihedral forces calculated
  igb=%d,     !
  gbsa=2,    ! recursive surface area algorithm (for postprocessing)
/
''' % (igb))
    script_F.close()

    os.chdir(self.args.dir['out'])
    import subprocess
    args_list = [self.args.FNs['sander'], '-O','-i',script_FN,'-o',out_FN, \
      '-p',self.args.FNs['prmtop'][moiety],'-c',self.args.FNs['inpcrd'][moiety], \
      '-y', AMBER_mdcrd_FN, '-r',script_FN+'.restrt']
    if debug:
      print ' '.join(args_list)
    p = subprocess.Popen(args_list)
    p.wait()

    F = open(out_FN, 'r')
    dat = F.read().strip().split(' BOND')
    F.close()

    dat.pop(0)
    if len(dat) > 0:
      # For the different models, all the terms are the same except for
      # EGB/EPB (every model is different)
      # ESURF versus ECAVITY + EDISPER
      # EEL (ALPB versus not)
      E = np.array([
        rec[:rec.find('\nminimization')].replace('1-4 ', '1-4').split()[1::3]
        for rec in dat
      ],
                   dtype=float) * MMTK.Units.kcal / MMTK.Units.mol
      if phase == 'sander_Gas':
        E = np.hstack((E, np.sum(E, 1)[..., None]))
      else:
        # Mark as nan to add the Gas energies later
        E = np.hstack((E, np.ones((E.shape[0], 1)) * np.nan))

      if not debug and os.path.isfile(script_FN):
        os.remove(script_FN)
      if os.path.isfile(script_FN + '.restrt'):
        os.remove(script_FN + '.restrt')

      if not debug and os.path.isfile(out_FN):
        os.remove(out_FN)
    else:
      E = np.array([np.inf] * 11)

    os.chdir(self.args.dir['start'])
    return E
    # AMBER ENERGY FIELDS:
    # For Gas phase:
    # 0. BOND 1. ANGLE 2. DIHEDRAL 3. VDWAALS 4. EEL
    # 5. HBOND 6. 1-4 VWD 7. 1-4 EEL 8. RESTRAINT
    # For GBSA phases:
    # 0. BOND 1. ANGLE 2. DIHEDRAL 3. VDWAALS 4. EEL
    # 5. EGB 6. 1-4 VWD 7. 1-4 EEL 8. RESTRAINT 9. ESURF
    # For PBSA phase:
    # 0. BOND 1. ANGLE 2. DIHEDRAL 3. VDWAALS 4. EEL
    # 5. EPB 6. 1-4 VWD 7. 1-4 EEL 8. RESTRAINT 9. ECAVITY 10. EDISPER

  def _get_elsize(self):
    # Calculates the electrostatic size of the receptor for ALPB calculations
    # Writes the coordinates in AMBER format
    pqr_FN = os.path.join(self.args.dir['CD'], 'receptor.pqr')
    if not os.path.isdir(self.args.dir['CD']):
      os.system('mkdir -p ' + self.args.dir['CD'])

    import AlGDock.IO
    IO_crd = AlGDock.IO.crd()
    factor = 1.0 / MMTK.Units.Ang
    IO_crd.write(self.args.FNs['inpcrd']['R'], factor*self.data['CD'].confs['receptor'], \
      'title', trajectory=False)

    # Converts the coordinates to a pqr file
    inpcrd_F = open(self.args.FNs['inpcrd']['R'], 'r')
    cdir = os.getcwd()
    import subprocess
    try:
      p = subprocess.Popen(\
        [self.args.FNs['ambpdb'], \
         '-p', os.path.relpath(self.args.FNs['prmtop']['R'], cdir), \
         '-pqr'], \
        stdin=inpcrd_F, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      (stdoutdata_ambpdb, stderrdata_ambpdb) = p.communicate()
      p.wait()
    except OSError:
      os.system('ls -ltr')
      print 'Command: ' + ' '.join([os.path.relpath(self.args.FNs['ambpdb'], cdir), \
         '-p', os.path.relpath(self.args.FNs['prmtop']['R'], cdir), \
         '-pqr'])
      print 'stdout:\n' + stdoutdata_ambpdb
      print 'stderr:\n' + stderrdata_ambpdb
    inpcrd_F.close()

    pqr_F = open(pqr_FN, 'w')
    pqr_F.write(stdoutdata_ambpdb)
    pqr_F.close()

    # Runs the pqr file through elsize
    p = subprocess.Popen(\
      [self.args.FNs['elsize'], os.path.relpath(pqr_FN, cdir)], \
      stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdoutdata_elsize, stderrdata_elsize) = p.communicate()
    p.wait()

    for FN in [pqr_FN]:
      if os.path.isfile(FN):
        os.remove(FN)
    try:
      elsize = float(stdoutdata_elsize.strip())
    except ValueError:
      print 'Command: ' + ' '.join([os.path.relpath(self.args.FNs['elsize'], cdir), \
       os.path.relpath(pqr_FN, cdir)])
      print stdoutdata_elsize
      print 'Error with elsize'
    return elsize

  def _gbnsr6_Energy(self,
                     confs,
                     moiety,
                     phase,
                     inpcrd_FN,
                     outputname,
                     debug=DEBUG,
                     reference=None):
    """
    Uses gbnsr6 (part of AmberTools)
    to calculate the energy of a set of configurations
    """
    # Prepare configurations for writing to crd file
    factor = 1.0 / MMTK.Units.Ang
    if (moiety.find('R') > -1):
      receptor_0 = factor * self.data['CD'].confs['receptor'][:self.top_RL.
                                                              L_first_atom, :]
      receptor_1 = factor * self.data['CD'].confs['receptor'][self.top_RL.
                                                              L_first_atom:, :]

    if not isinstance(confs, list):
      confs = [confs]

    if (moiety.find('R') > -1):
      if (moiety.find('L') > -1):
        full_confs = [np.vstack((receptor_0, \
          conf[self.top.prmtop_atom_order_L,:]/MMTK.Units.Ang, \
          receptor_1)) for conf in confs]
      else:
        full_confs = [factor * self.data['CD'].confs['receptor']]
    else:
      full_confs = [conf[self.top.prmtop_atom_order_L,:]/MMTK.Units.Ang \
        for conf in confs]

    # Set up directory
    inpcrdFN = os.path.abspath(inpcrd_FN)
    gbnsr6_dir = os.path.dirname(inpcrd_FN)
    os.system('mkdir -p ' + gbnsr6_dir)
    os.chdir(gbnsr6_dir)
    cdir = os.getcwd()

    # Write gbnsr6 script
    chagb = 0 if phase.find('Still') > -1 else 1
    alpb = 1 if moiety.find(
      'R') > -1 else 0  # ALPB ineffective with small solutes
    gbnsr6_in_FN = moiety + 'gbnsr6.in'
    gbnsr6_in_F = open(gbnsr6_in_FN, 'w')
    gbnsr6_in_F.write("""gbnsr6
&cntrl
  inp=1
/
&gb
  alpb=%d,
  chagb=%d
/
""" % (alpb, chagb))
    gbnsr6_in_F.close()

    args_list = [self.args.FNs['gbnsr6'], \
      '-i', os.path.relpath(gbnsr6_in_FN, cdir), \
      '-o', 'stdout', \
      '-p', os.path.relpath(self.args.FNs['prmtop'][moiety], cdir), \
      '-c', os.path.relpath(inpcrd_FN, cdir)]
    if debug:
      print ' '.join(args_list)

    # Write coordinates, run gbnsr6, and store energies
    import subprocess
    import AlGDock.IO
    IO_crd = AlGDock.IO.crd()

    E = []
    for full_conf in full_confs:
      # Writes the coordinates in AMBER format
      IO_crd.write(inpcrd_FN, full_conf, 'title', trajectory=False)

      # Runs gbnsr6
      import subprocess
      p = subprocess.Popen(args_list, \
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      (stdoutdata, stderrdata) = p.communicate()
      p.wait()

      recs = stdoutdata.strip().split(' BOND')
      if len(recs) > 1:
        rec = recs[1]
        E.append(rec[:rec.find('\n -----')].replace('1-4 ',
                                                    '1-4').split()[1::3])
      else:
        self.log.tee("  error has occured in gbnsr6 after %d snapshots" %
                     len(E))
        self.log.tee("  prmtop was " + self.args.FNs['prmtop'][moiety])
        self.log.tee("  --- stdout:")
        self.log.tee(stdoutdata)
        self.log.tee("  --- stderr:")
        self.log.tee(stderrdata)

    E = np.array(E, dtype=float) * MMTK.Units.kcal / MMTK.Units.mol
    E = np.hstack((E, np.ones((E.shape[0], 1)) * np.nan))

    os.chdir(self.args.dir['start'])
    if not debug:
      os.system('rm -rf ' + gbnsr6_dir)
    return E
    # For gbnsr6 phases:
    # 0. BOND 1. ANGLE 2. DIHED 3. 1-4 NB 4. 1-4 EEL
    # 5. VDWAALS 6. EELEC 7. EGB 8. RESTRAINT 9. ESURF

  def _setup_OpenMM(self, moiety, phase):
    if not hasattr(self, '_OpenMM_sims'):
      self._OpenMM_sims = {}
    key = moiety + phase
    if not key in self._OpenMM_sims.keys():
      import simtk.openmm
      import simtk.openmm.app as OpenMM_app
      prmtop = OpenMM_app.AmberPrmtopFile(self.args.FNs['prmtop'][moiety])
      inpcrd = OpenMM_app.AmberInpcrdFile(self.args.FNs['inpcrd'][moiety])
      OMM_system = prmtop.createSystem(nonbondedMethod=OpenMM_app.NoCutoff, \
        constraints=None, implicitSolvent={
          'OpenMM_Gas':None,
          'OpenMM_GBn':OpenMM_app.GBn,
          'OpenMM_GBn2':OpenMM_app.GBn2,
          'OpenMM_HCT':OpenMM_app.HCT,
          'OpenMM_OBC1':OpenMM_app.OBC1,
          'OpenMM_OBC2':OpenMM_app.OBC2}[phase])
      # Set receptor atom mass to zero to facilitate future minimization
      if moiety == 'R':
        for i in range(OMM_system.getNumParticles()):
          OMM_system.setParticleMass(i, 0)
      elif moiety == 'RL':
        for i in range(self.top_RL.L_first_atom) + \
            range(self.top_RL.L_first_atom + self.top.universe.numberOfAtoms(), \
            OMM_system.getNumParticles()):
          OMM_system.setParticleMass(i, 0)
      dummy_integrator = simtk.openmm.LangevinIntegrator(300*simtk.unit.kelvin, \
        1/simtk.unit.picosecond, 0.002*simtk.unit.picoseconds)
      self._OpenMM_sims[key] = OpenMM_app.Simulation(prmtop.topology, \
        OMM_system, dummy_integrator)

  def _OpenMM_Energy(self, confs, moiety, phase, traj_FN=None, \
      outputname=None, debug=DEBUG, reference=None):
    key = moiety + phase
    self._setup_OpenMM(moiety, phase)  # Set up the simulation

    # Prepare the conformations by combining with the receptor if necessary
    if (moiety.find('R') > -1):
      receptor_0 = self.data['CD'].confs['receptor'][:self.top_RL.
                                                     L_first_atom, :]
      receptor_1 = self.data['CD'].confs['receptor'][self.top_RL.
                                                     L_first_atom:, :]
    if not isinstance(confs, list):
      confs = [confs]
    if (moiety.find('R') > -1):
      if (moiety.find('L') > -1):
        confs = [np.vstack((receptor_0, \
          conf[self.top.prmtop_atom_order_L,:], \
          receptor_1)) for conf in confs]
      else:
        confs = [self.data['CD'].confs['receptor']]
    else:
      confs = [conf[self.top.prmtop_atom_order_L, :] for conf in confs]

    import simtk.unit
    # Calculate the energies
    E = []
    for conf in confs:
      self._OpenMM_sims[key].context.setPositions(conf)
      s = self._OpenMM_sims[key].context.getState(getEnergy=True)
      E.append(
        [0.,
         s.getPotentialEnergy() / simtk.unit.kilojoule * simtk.unit.mole])
    return np.array(E, dtype=float) * MMTK.Units.kJ / MMTK.Units.mol

  def _APBS_Energy(self,
                   confs,
                   moiety,
                   phase,
                   pqr_FN,
                   outputname,
                   debug=DEBUG,
                   reference=None):
    """
    Uses APBS to calculate the solvation energy of a set of configurations
    Units are the MMTK standard, kJ/mol
    """
    # Prepare configurations for writing to crd file
    factor = 1.0 / MMTK.Units.Ang
    if (moiety.find('R') > -1):
      receptor_0 = factor * self.data['CD'].confs['receptor'][:self.top_RL.
                                                              L_first_atom, :]
      receptor_1 = factor * self.data['CD'].confs['receptor'][self.top_RL.
                                                              L_first_atom:, :]

    if not isinstance(confs, list):
      confs = [confs]

    if (moiety.find('R') > -1):
      if (moiety.find('L') > -1):
        full_confs = [np.vstack((receptor_0, \
          conf[self.top.prmtop_atom_order_L,:]/MMTK.Units.Ang, \
          receptor_1)) for conf in confs]
      else:
        full_confs = [factor * self.data['CD'].confs['receptor']]
    else:
      full_confs = [conf[self.top.prmtop_atom_order_L,:]/MMTK.Units.Ang \
        for conf in confs]

    # Write coordinates, run APBS, and store energies
    apbs_dir = os.path.abspath(pqr_FN)[:-4]
    os.system('mkdir -p ' + apbs_dir)
    os.chdir(apbs_dir)
    pqr_FN = os.path.join(apbs_dir, 'in.pqr')

    import subprocess
    import AlGDock.IO
    IO_crd = AlGDock.IO.crd()

    E = []
    for full_conf in full_confs:
      # Writes the coordinates in AMBER format
      inpcrd_FN = pqr_FN[:-4] + '.crd'
      IO_crd.write(inpcrd_FN, full_conf, 'title', trajectory=False)

      # Converts the coordinates to a pqr file
      inpcrd_F = open(inpcrd_FN, 'r')
      cdir = os.getcwd()
      p = subprocess.Popen(\
        [os.path.relpath(self.args.FNs['ambpdb'], cdir), \
         '-p', os.path.relpath(self.args.FNs['prmtop'][moiety], cdir), \
         '-pqr'], \
        stdin=inpcrd_F, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      (stdoutdata_ambpdb, stderrdata_ambpdb) = p.communicate()
      p.wait()
      inpcrd_F.close()

      pqr_F = open(pqr_FN, 'w')
      pqr_F.write(stdoutdata_ambpdb)
      pqr_F.close()

      # Writes APBS script
      apbs_in_FN = moiety + 'apbs-mg-manual.in'
      apbs_in_F = open(apbs_in_FN, 'w')
      apbs_in_F.write('READ\n  mol pqr {0}\nEND\n'.format(pqr_FN))

      for sdie in [80.0, 1.0]:
        if moiety == 'L':
          min_xyz = np.array([min(full_conf[a, :]) for a in range(3)])
          max_xyz = np.array([max(full_conf[a, :]) for a in range(3)])
          mol_range = max_xyz - min_xyz
          mol_center = (min_xyz + max_xyz) / 2.

          def roundUpDime(x):
            return (np.ceil((x.astype(float) - 1) / 32) * 32 + 1).astype(int)

          focus_spacing = 0.5
          focus_dims = roundUpDime(mol_range * LFILLRATIO / focus_spacing)
          args = zip(['mdh'], [focus_dims], [mol_center], [focus_spacing])
        else:
          args = zip(['mdh', 'focus'], self._apbs_grid['dime'],
                     self._apbs_grid['gcent'], self._apbs_grid['spacing'])
        for (bcfl, dime, gcent, grid) in args:
          apbs_in_F.write('''ELEC mg-manual
  bcfl {0} # multiple debye-huckel boundary condition
  chgm spl4 # quintic B-spline charge discretization
  dime {1[0]} {1[1]} {1[2]}
  gcent {2[0]} {2[1]} {2[2]}
  grid {3} {3} {3}
  lpbe # Linearized Poisson-Boltzmann
  mol 1
  pdie 1.0
  sdens 10.0
  sdie {4}
  srad 1.4
  srfm smol # Smoothed dielectric and ion-accessibility coefficients
  swin 0.3
  temp 300.0
  calcenergy total
END
'''.format(bcfl, dime, gcent, grid, sdie))
      apbs_in_F.write('quit\n')
      apbs_in_F.close()

      # Runs APBS
      #      TODO: Control the number of threads. This doesn't seem to do anything.
      #      if self.args.cores==1:
      #        os.environ['OMP_NUM_THREADS']='1'
      p = subprocess.Popen([self.args.FNs['apbs'], apbs_in_FN], \
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      (stdoutdata, stderrdata) = p.communicate()
      p.wait()

      apbs_energy = [float(line.split('=')[-1][:-7]) \
        for line in stdoutdata.split('\n') \
        if line.startswith('  Total electrostatic energy')]
      if moiety == 'L' and len(apbs_energy) == 2:
        polar_energy = apbs_energy[0] - apbs_energy[1]
      elif len(apbs_energy) == 4:
        polar_energy = apbs_energy[1] - apbs_energy[3]
      else:
        # An error has occured in APBS
        polar_energy = np.inf
        self.log.tee("  error has occured in APBS after %d snapshots" % len(E))
        self.log.tee("  prmtop was " + self.args.FNs['prmtop'][moiety])
        self.log.tee("  --- ambpdb stdout:")
        self.log.tee(stdoutdata_ambpdb)
        self.log.tee("  --- ambpdb stderr:")
        self.log.tee(stderrdata_ambpdb)
        self.log.tee("  --- APBS stdout:")
        self.log.tee(stdoutdata)
        self.log.tee("  --- APBS stderr:")
        self.log.tee(stderrdata)

      # Runs molsurf to calculate Connolly surface
      apolar_energy = np.inf
      p = subprocess.Popen([self.args.FNs['molsurf'], pqr_FN, '1.4'], \
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      (stdoutdata, stderrdata) = p.communicate()
      p.wait()

      for line in stdoutdata.split('\n'):
        if line.startswith('surface area ='):
          apolar_energy = float(line.split('=')[-1]) * \
            0.0072 * MMTK.Units.kcal/MMTK.Units.mol

      if debug:
        molsurf_out_FN = moiety + 'molsurf-mg-manual.out'
        molsurf_out_F = open(molsurf_out_FN, 'w')
        molsurf_out_F.write(stdoutdata)
        molsurf_out_F.close()
      else:
        for FN in [inpcrd_FN, pqr_FN, apbs_in_FN, 'io.mc']:
          os.remove(FN)

      E.append([polar_energy, apolar_energy, np.nan])

      if np.isinf(polar_energy) or np.isinf(apolar_energy):
        break

    os.chdir(self.args.dir['start'])
    if not debug:
      os.system('rm -rf ' + apbs_dir)
    return np.array(E, dtype=float) * MMTK.Units.kJ / MMTK.Units.mol

  def _get_APBS_grid_spacing(self, RFILLRATIO=RFILLRATIO):
    factor = 1.0 / MMTK.Units.Ang

    def roundUpDime(x):
      return (np.ceil((x.astype(float) - 1) / 32) * 32 + 1).astype(int)

    self.system.setParams({'MM': True, 'ELE': 1})
    gd = self._forceFields['ELE'].grid_data
    focus_dims = roundUpDime(gd['counts'])
    focus_center = factor * (gd['counts'] * gd['spacing'] / 2. + gd['origin'])
    focus_spacing = factor * gd['spacing'][0]

    min_xyz = np.array([
      min(factor * self.data['CD'].confs['receptor'][a, :]) for a in range(3)
    ])
    max_xyz = np.array([
      max(factor * self.data['CD'].confs['receptor'][a, :]) for a in range(3)
    ])
    mol_range = max_xyz - min_xyz
    mol_center = (min_xyz + max_xyz) / 2.

    # The full grid spans RFILLRATIO times the range of the receptor
    # and the focus grid, whatever is larger
    full_spacing = 1.0
    full_min = np.minimum(mol_center - mol_range/2.*RFILLRATIO, \
                          focus_center - focus_dims*focus_spacing/2.*RFILLRATIO)
    full_max = np.maximum(mol_center + mol_range/2.*RFILLRATIO, \
                          focus_center + focus_dims*focus_spacing/2.*RFILLRATIO)
    full_dims = roundUpDime((full_max - full_min) / full_spacing)
    full_center = (full_min + full_max) / 2.

    self._apbs_grid = {\
      'dime':[full_dims, focus_dims], \
      'gcent':[full_center, focus_center], \
      'spacing':[full_spacing, focus_spacing]}

  def _combine_MM_and_solvent(self, E, toParse=None):
    if toParse is None:
      toParse = [k for k in E.keys() \
        if (E[k] is not None) and (len(np.array(E[k]).shape)==2)]
    for key in toParse:
      if np.isnan(E[key][:, -1]).all():
        E[key] = E[key][:, :-1]
        if key.find('sander') > -1:
          prefix = key.split('_')[0][:-6]
          for c in [0, 1, 2, 6, 7]:
            E[key][:, c] = E[prefix + 'sander_Gas'][:, c]
        elif key.find('gbnsr6') > -1:
          prefix = key.split('_')[0][:-6]
          for (gbnsr6_ind, sander_ind) in [(0, 0), (1, 1), (2, 2), (3, 6),
                                           (5, 3)]:
            E[key][:, gbnsr6_ind] = E[prefix + 'sander_Gas'][:, sander_ind]
        elif key.find('APBS_PBSA'):
          prefix = key[:-9]
          totalMM = np.transpose(np.atleast_2d(E[prefix + 'NAMD_Gas'][:, -1]))
          E[key] = np.hstack((E[key], totalMM))
        E[key] = np.hstack((E[key], np.sum(E[key], 1)[..., None]))

  def _write_traj(self, traj_FN, confs, moiety, \
      title='', factor=1.0/MMTK.Units.Ang):
    """
    Writes a trajectory file
    """

    if traj_FN is None:
      return
    if traj_FN.endswith('.pqr'):
      return
    if traj_FN.endswith('.crd'):
      return
    if os.path.isfile(traj_FN):
      return

    traj_dir = os.path.dirname(os.path.abspath(traj_FN))
    if not os.path.isdir(traj_dir):
      os.system('mkdir -p ' + traj_dir)

    import AlGDock.IO
    if traj_FN.endswith('.dcd'):
      IO_dcd = AlGDock.IO.dcd(self.top.molecule,
        ligand_atom_order = self.top.prmtop_atom_order_L, \
        receptorConf = self.data['CD'].confs['receptor'], \
        ligand_first_atom = self.top_RL.L_first_atom)
      IO_dcd.write(traj_FN,
                   confs,
                   includeReceptor=(moiety.find('R') > -1),
                   includeLigand=(moiety.find('L') > -1))
    elif traj_FN.endswith('.mdcrd'):
      if (moiety.find('R') > -1):
        receptor_0 = factor * self.data['CD'].confs[
          'receptor'][:self.top_RL.L_first_atom, :]
        receptor_1 = factor * self.data['CD'].confs['receptor'][
          self.top_RL.L_first_atom:, :]

      if not isinstance(confs, list):
        confs = [confs]
      if (moiety.find('R') > -1):
        if (moiety.find('L') > -1):
          confs = [np.vstack((receptor_0, \
            conf[self.top.prmtop_atom_order_L,:]/MMTK.Units.Ang, \
            receptor_1)) for conf in confs]
        else:
          confs = [factor * self.data['CD'].confs['receptor']]
      else:
        confs = [conf[self.top.prmtop_atom_order_L,:]/MMTK.Units.Ang \
          for conf in confs]

      import AlGDock.IO
      IO_crd = AlGDock.IO.crd()
      IO_crd.write(traj_FN, confs, title, trajectory=True)
      self.log.tee("  wrote %d configurations to %s" % (len(confs), traj_FN))
    else:
      raise Exception('Unknown trajectory type')

  def _clear_f_RL(self):
    # stats_RL will include internal energies, interaction energies,
    # the cycle by which the bound state is equilibrated,
    # the mean acceptance probability between replica exchange neighbors,
    # and the rmsd, if applicable
    phase_f_RL_keys = \
      [phase+'_solv' for phase in self.args.params['CD']['phases']]

    # Initialize variables as empty lists
    stats_RL = [('u_K_'+FF,[]) \
      for FF in ['ligand','sampled']+self.args.params['CD']['phases']]
    stats_RL += [('Psi_'+FF,[]) \
      for FF in ['grid']+self.args.params['CD']['phases']]
    stats_RL += [(item,[]) \
      for item in ['equilibrated_cycle','cum_Nclusters','mean_acc','rmsd']]
    self.stats_RL = dict(stats_RL)
    self.stats_RL['protocol'] = self.data['CD'].protocol
    # Free energy components
    self.f_RL = dict([(key,[]) \
      for key in ['grid_MBAR'] + phase_f_RL_keys])
    # Binding PMF estimates
    self.B = {'MMTK_MBAR': []}
    for phase in self.args.params['CD']['phases']:
      for method in ['min_Psi', 'mean_Psi', 'EXP', 'MBAR']:
        self.B[phase + '_' + method] = []

    # Store empty list
    if self.args.params['CD']['pose'] == -1:
      f_RL_FN = os.path.join(self.args.dir['CD'], 'f_RL.pkl.gz')
    else:
      f_RL_FN = os.path.join(self.args.dir['CD'], \
        'f_RL_pose%03d.pkl.gz'%self.args.params['CD']['pose'])
    if hasattr(self, 'run_type') and (not self.log.run_type.startswith('timed')):
      self.log.tee(
        write_pkl_gz(f_RL_FN, (self.f_L, self.stats_RL, self.f_RL, self.B)))

  def save(self, p, keys=['progress', 'data']):
    """Saves results

    Parameters
    ----------
    p : str
      The process, either 'BC' or 'CD'
    keys : list of str
      Save the progress, the data, or both
    """
    if 'progress' in keys:
      self.log.tee(self.args.save_pkl_gz(p, self.data[p]))
    if 'data' in keys:
      self.log.tee(self.data[p].save_pkl_gz())

  def __del__(self):
    if (not DEBUG) and len(self.args.toClear) > 0:
      print "\n>>> Clearing files"
      for FN in self.args.toClear:
        if os.path.isfile(FN):
          os.remove(FN)
          print '  removed ' + os.path.relpath(FN, self.args.dir['start'])

if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(
    description=
    'Molecular docking with adaptively scaled alchemical interaction grids')

  for key in arguments.keys():
    parser.add_argument('--' + key, **arguments[key])
  args = parser.parse_args()

  if args.run_type in ['render_docked', 'render_intermediates']:
    from AlGDock.BindingPMF_plots import BPMF_plots
    self = BPMF_plots(**vars(args))
  else:
    self = BPMF(**vars(args))
