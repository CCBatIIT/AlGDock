#!/usr/bin/env python

# TODO: Free energy of external confinement for poseBPMFs

import os
import cPickle as pickle
import gzip
import copy
import sys


import AlGDock.IO 
from AlGDock.IO import load_pkl_gz
from AlGDock.IO import write_pkl_gz
from AlGDock.logger import NullDevice

import time
import numpy as np

from collections import OrderedDict

from AlGDock import dictionary_tools

try:
  import MMTK
  import MMTK.Units
  from MMTK.ParticleProperties import Configuration
  from MMTK.ForceFields import ForceField
except ImportError:
  MMTK = None

try:
  import openmm
  import openmm.unit as unit
  from openmm.app import AmberPrmtopFile, AmberInpcrdFile, Simulation, NoCutoff
  from openmm import *
  import parmed as pmd # need to install ParmEd

except ImportError:
  openmm = None

import Scientific
try:
  from Scientific._vector import Vector
except:
  from Scientific.Geometry.VectorModule import Vector

import pymbar.timeseries

import multiprocessing
from multiprocessing import Process
import arguments
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
    """Parses the input arguments and runs the requested calculation"""

    #         mod_path = os.path.join(os.path.dirname(a.__file__), 'BindingPMF.py')
    #         print """###########
    # # AlGDock #
    # ###########
    # Molecular docking with adaptively scaled alchemical interaction grids
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

    if not 'max_time' in kwargs.keys():
      kwargs['max_time'] = None
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

    self.run(kwargs['run_type'])

  def _setup(self):
    """Creates an MMTK InfiniteUniverse and adds the ligand"""
    if MMTK:
      from AlGDock.topology import TopologyMMTK
      self.top = TopologyMMTK(self.args)
      self.top_RL = TopologyMMTK(self.args, includeReceptor=True)
    else:
      from AlGDock.topology import TopologyUsingOpenMM
      self.top = TopologyUsingOpenMM(self.args)
      self.top_RL = TopologyUsingOpenMM(self.args, includeReceptor=True)
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
              if MMTK:
                self.top.universe.setConfiguration(
                Configuration(self.top.universe, conf))
                coms.append(np.array(self.top.universe.centerOfMass()))

              else:
                self.top.OMM_simulaiton.context.setPositions(conf)
                structure = pmd.openmm.load_topology(self.top.molecule, self.top.OMM_system)
                coms.append(np.array(structure.center_of_mass) * 0.1) # Convert from Å to nm

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
          if MMTK:
            self.top.universe.setConfiguration(
            Configuration(self.top.universe, confs[-1]))
          else:
            self.top.OMM_simulaiton.context.setPositions(confs[-1])

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
      if MMTK:
        lig_crd = complex_crd[self.top_RL.L_first_atom:self.top_RL.L_first_atom + \
        self.top.universe.numberOfAtoms(),:]
        self.data['CD'].confs['receptor'] = np.vstack(\
        (complex_crd[:self.top_RL.L_first_atom,:],\
         complex_crd[self.top_RL.L_first_atom + self.top.universe.numberOfAtoms():,:]))
      else:
        natoms = sum(1 for _ in self.top.OMM_simulaiton.topology.atoms())
        lig_crd = complex_crd[self.top_RL.L_first_atom:self.top_RL.L_first_atom + natoms,:]
        self.data['CD'].confs['receptor'] = np.vstack(\
        (complex_crd[:self.top_RL.L_first_atom,:],\
         complex_crd[self.top_RL.L_first_atom + natoms:,:]))

    elif self.args.FNs['inpcrd']['L'] is not None:
      self.data['CD'].confs['receptor'] = None
      if os.path.isfile(self.args.FNs['inpcrd']['L']):
        lig_crd = IO_crd.read(self.args.FNs['inpcrd']['L'], multiplier=0.1)
    else:
      lig_crd = None

    if lig_crd is not None:
      if MMTK:
        self.data['CD'].confs['ligand'] = lig_crd[self.top.inv_prmtop_atom_order_L, :]
        self.top.universe.setConfiguration(\
        Configuration(self.top.universe,self.data['CD'].confs['ligand']))
      else:
        self.data['CD'].confs['ligand'] = lig_crd
        self.top.OMM_simulaiton.context.setPositions(self.data['CD'].confs['ligand'])
      if MMTK:
        if self.top_RL.universe is not None:
          self.top_RL.universe.setConfiguration(\
          Configuration(self.top_RL.universe, \
          np.vstack((self.data['CD'].confs['receptor'],self.data['CD'].confs['ligand']))))
      else:
        if self.top_RL.OMM_simulaiton is not None:
          self.top.OMM_simulaiton.context.setPositions(self.data['CD'], \
          np.vstack((self.data['CD'].confs['receptor'],self.data['CD'].confs['ligand'])))

    if self.args.params['CD']['rmsd'] is not False:
      if self.args.params['CD']['rmsd'] is True:
        if lig_crd is not None:
          rmsd_crd = lig_crd[self.top.inv_prmtop_atom_order_L, :]
        else:
          raise Exception('Reference structure for rmsd calculations unknown')
      else:
        if MMTK:
          rmsd_crd = IO_crd.read(self.args.params['CD']['rmsd'], \
          natoms=self.top.universe.numberOfAtoms(), multiplier=0.1)
          rmsd_crd = rmsd_crd[self.top.inv_prmtop_atom_order_L, :]
        else:
          n_atoms = sum(1 for _ in self.top.OMM_simulaiton.topology.atoms())
          rmsd_crd = IO_crd.read(self.args.params['CD']['rmsd'], natoms=n_atoms, multiplier=0.1)
      self.data['CD'].confs['rmsd'] = rmsd_crd

      self.get_rmsds.set_ref_configuration(self.data['CD'].confs['rmsd'])

    # If configurations are being rescored, start with a docked structure
    (confs, Es) = self._get_confs_to_rescore(site=False, minimize=False)
    if len(confs) > 0:
      if MMTK:
        self.top.universe.setConfiguration(
        Configuration(self.top.universe, confs[-1]))
      else:
        self.top.OMM_simulaiton.context.setPositions(confs[-1])

    from AlGDock.simulation_iterator import SimulationIterator
    self.iterator = SimulationIterator(self.args, self.top, self.system)

    # Load progress
    from AlGDock.postprocessing import Postprocessing
    Postprocessing(self.args, self.log, self.top, self.top_RL, self.system, self.data, self.save).run(readOnly=True)

    self.calc_f_L(readOnly=True)
    self.calc_f_RL(readOnly=True)

    if self.args.random_seed > 0:
      np.random.seed(self.args.random_seed)

  def run(self, run_type):
    from AlGDock.postprocessing import Postprocessing

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
      Postprocessing(self.args, self.log, self.top, self.top_RL, self.system, self.data, self.save).run([('BC', -1, -1, 'L')])
      self.calc_f_L()
    elif run_type == 'initial_CD':
      self.initial_CD()
    elif run_type == 'CD':  # Sample the CD process
      self.sim_process('CD')
      Postprocessing(self.args, self.log, self.top, self.top_RL, self.system, self.data, self.save).run()
      self.calc_f_RL()
      # self.targeted_FEP()
    elif run_type == 'timed':  # Timed replica exchange sampling
      BC_complete = self.sim_process('BC')
      if BC_complete:
        pp_complete = Postprocessing(self.args, self.log, self.top, self.top_RL, self.system, self.data, self.save).run([('BC', -1, -1, 'L')])
        if pp_complete:
          self.calc_f_L()
          CD_complete = self.sim_process('CD')
          if CD_complete:
            pp_complete = Postprocessing(self.args, self.log, self.top, self.top_RL, self.system, self.data, self.save).run()
            if pp_complete:
              self.calc_f_RL()
              # self.targeted_FEP()
    elif run_type == 'timed_BC':  # Timed BC only
      BC_complete = self.sim_process('BC')
      if BC_complete:
        pp_complete = Postprocessing(self.args, self.log, self.top, self.top_RL, self.system, self.data, self.save).run([('BC', -1, -1, 'L')])
        if pp_complete:
          self.calc_f_L()
    elif run_type == 'timed_CD':  # Timed CD only
      CD_complete = self.sim_process('CD')
      if CD_complete:
        pp_complete = Postprocessing(self.args, self.log, self.top, self.top_RL, self.system, self.data, self.save).run()
        if pp_complete:
          self.calc_f_RL()
          # self.targeted_FEP()
    elif run_type == 'postprocess':  # Postprocessing
      Postprocessing(self.args, self.log, self.top, self.top_RL, self.system, self.data, self.save).run()
    elif run_type == 'redo_postprocess':
      Postprocessing(self.args, self.log, self.top, self.top_RL, self.system, self.data, self.save).run(redo_CD=True)
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
      Postprocessing(self.args, self.log, self.top, self.top_RL, self.system, self.data, self.save).run([('BC', -1, -1, 'L')])
      self.calc_f_L()
      self.sim_process('CD')
      Postprocessing(self.args, self.log, self.top, self.top_RL, self.system, self.data, self.save).run()
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
      MBAR = self.run_MBAR(u_kln, N_k)[0]
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
    from AlGDock.postprocessing import Postprocessing
    pp_complete = Postprocessing(self.args, self.log, self.top, self.top_RL, self.system, self.data, self.save).run([('BC', -1, -1, 'L')])
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
      MBAR = self.run_MBAR(u_kln, N_k)[0]
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
    from AlGDock.postprocessing import Postprocessing
    pp_complete = Postprocessing(self.args, self.log, self.top, self.top_RL, self.system, self.data, self.save).run()
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
    self.data['CD'].confs['starting_poses'] = None
    from AlGDock.postprocessing import Postprocessing
    pp_complete = Postprocessing(self.args, self.log, self.top, self.top_RL, self.system, self.data, self.save).run([('original', 0, 0, 'R')])

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
        from AlGDock.replica_exchange import ReplicaExchange
        ReplicaExchange(self.args, self.log, self.top, self.system,
                      self.iterator, self.data, self.save, self._u_kln).run(process)
        self.SIRS(process)
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
        from AlGDock.replica_exchange import ReplicaExchange
        ReplicaExchange(self.args, self.log, self.top, self.system,
                      self.iterator, self.data, self.save, self._u_kln).run('BC')
        self.SIRS(process)
        cycle_times.append(self.log.timeSince('repX cycle'))
        if not self.log.isTimeForTask(cycle_times):
          return False
        E_MM = []
        for k in range(len(self.data['BC'].Es[0])):
          E_MM += list(self.data['BC'].Es[0][k]['MM'])

    # Clear evaluators to save memory
    self.system.clear_evaluators()

    return True  # The process has completed

  def SIRS(self, process):
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

    protocol = self.data[process].protocol
    Es_repX = [[copy.deepcopy(self.data[process].Es[k][-1])] for k in range(len(protocol))]
    (u_kln, N_k) = self._u_kln(Es_repX, protocol)

    # This is a more direct way to get the weights
    from pymbar.utils import kln_to_kn
    u_kn = kln_to_kn(u_kln, N_k=N_k)

    from pymbar.utils import logsumexp
    log_denominator_n = logsumexp(f_k - u_kn.T, b=N_k, axis=1)
    logW = f_k - u_kn.T - log_denominator_n[:, np.newaxis]
    W_nl = np.exp(logW)
    for k in range(len(protocol)):
      W_nl[:, k] = W_nl[:, k] / np.sum(W_nl[:, k])

    # This is for conversion to 2 indicies: state and snapshot
    cum_N_state = np.cumsum([0] + list(N_k))

    def linear_index_to_snapshot_index(ind):
      state_index = list(ind < cum_N_state).index(True) - 1
      nis_index = ind - cum_N_state[state_index]
      return (state_index, nis_index)

    # Selects new replica exchange snapshots
    confs_repX = self.data[process].confs['last_repX']
    self.data[process].confs['replicas'] = []
    for k in range(len(protocol)):
      (s,n) = linear_index_to_snapshot_index(\
        np.random.choice(range(W_nl.shape[0]), size = 1, p = W_nl[:,k])[0])
      self.data[process].confs['replicas'].append(np.copy(confs_repX[s][n]))

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
    weights = self.run_MBAR(u_kln, N_k, augmented=True)[1][:, -1]
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

  def run_MBAR(self, u_kln, N_k, augmented=False):
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

  def _u_kln(self, eTs, protocol, noBeta=False):
    """
    Computes a reduced potential energy matrix.  k is the sampled state.  l is the state for which energies are evaluated.

    Input:
    eT is a
      -dictionary (of mapped energy terms) of numpy arrays (over states)
      -list (over states) of dictionaries (of mapped energy terms) of numpy arrays (over configurations), or a
      -list (over states) of lists (over cycles) of dictionaries (of mapped energy terms) of numpy arrays (over configurations)
    protocol is a list of thermodynamic states
    noBeta means that the energy will not be divided by RT

    Output: u_kln or (u_kln, N_k)
    u_kln is the matrix (as a numpy array)
    N_k is an array of sample sizes
    """
    L = len(protocol)

    addMM = ('MM' in protocol[0].keys()) and (protocol[0]['MM'])
    addSite = ('site' in protocol[0].keys()) and (protocol[0]['site'])
    probe_keys = ['MM','k_angular_ext','k_spatial_ext','k_angular_int'] + \
      scalables
    probe_key = [key for key in protocol[0].keys() if key in probe_keys][0]

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
          if scalable in protocol[l].keys():
            E += protocol[l][scalable] * eTs[scalable]
        for key in ['k_angular_ext', 'k_spatial_ext', 'k_angular_int']:
          if key in protocol[l].keys():
            E += protocol[l][key] * eTs[key]
        if noBeta:
          u_kln.append(E)
        else:
          u_kln.append(E / (R * protocol[l]['T']))
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
            if scalable in protocol[l].keys():
              E += protocol[l][scalable] * eTs[k][scalable]
          for key in ['k_angular_ext', 'k_spatial_ext', 'k_angular_int']:
            if key in protocol[l].keys():
              E += protocol[l][key] * eTs[k][key]
          if noBeta:
            u_kln[k, l, :N_k[k]] = E
          else:
            u_kln[k, l, :N_k[k]] = E / (R * protocol[l]['T'])
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
            if scalable in protocol[l].keys():
              E += protocol[l][scalable]*np.concatenate([eTs[k][c][scalable] \
                for c in range(C)])
          for key in ['k_angular_ext', 'k_spatial_ext', 'k_angular_int']:
            if key in protocol[l].keys():
              E += protocol[l][key]*np.concatenate([eTs[k][c][key] \
                for c in range(C)])
          if noBeta:
            u_kln[k, l, :N_k[k]] = E
          else:
            u_kln[k, l, :N_k[k]] = E / (R * protocol[l]['T'])

    if (K == 1) and (L == 1):
      return u_kln.ravel()
    else:
      return (u_kln, N_k)

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

  for key in arguments.args.keys():
    parser.add_argument('--' + key, **arguments.args[key])
  args = parser.parse_args()

  if args.run_type in ['render_docked', 'render_intermediates']:
    from AlGDock.BindingPMF_plots import BPMF_plots
    self = BPMF_plots(**vars(args))
  else:
    self = BPMF(**vars(args))
