#!/usr/bin/env python

import os # Miscellaneous operating system interfaces
from os.path import abspath
from os.path import exists
from os.path import isfile
from os.path import join
import cPickle as pickle
import gzip
import copy

import time
import numpy as np

import MMTK
import MMTK.Units
from MMTK.ParticleProperties import Configuration

import Scientific
from Scientific._vector import Vector  # @UnresolvedImport

import AlGDock as a
import pymbar.timeseries

from multiprocessing import Process
import multiprocessing

# Constants
R = 8.3144621 * MMTK.Units.J / MMTK.Units.mol / MMTK.Units.K
T_HIGH = 600.0 * MMTK.Units.K
T_TARGET = 300.0 * MMTK.Units.K
RT_HIGH = R * T_HIGH
RT_TARGET = R * T_TARGET

term_map = {
  'cosine dihedral angle':'MM',
  'electrostatic/pair sum':'MM',
  'harmonic bond':'MM',
  'harmonic bond angle':'MM',
  'Lennard-Jones':'MM',
  'site':'site',
  'sLJr':'sLJr',
  'sELE':'sELE',
  'sLJa':'sLJa',
  'LJr':'LJr',
  'LJa':'LJa',
  'ELE':'ELE',
  'electrostatic':'misc'}

allowed_phases = ['Gas','GBSA','NAMD_Gas','NAMD_GBSA']

########################
# Auxilliary functions #
########################

def random_rotate():
  """
  Return a random rotation matrix
  """
  u = np.random.uniform(size=3)

  # Random quaternion
  q = np.array([np.sqrt(1-u[0])*np.sin(2*np.pi*u[1]),
               np.sqrt(1-u[0])*np.cos(2*np.pi*u[1]),
               np.sqrt(u[0])*np.sin(2*np.pi*u[2]),
               np.sqrt(u[0])*np.cos(2*np.pi*u[2])])
  
  # Convert the quaternion into a rotation matrix 
  rotMat = np.array([[q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3],
                     2*q[1]*q[2] - 2*q[0]*q[3],
                     2*q[1]*q[3] + 2*q[0]*q[2]],
                    [2*q[1]*q[2] + 2*q[0]*q[3],
                     q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3],
                     2*q[2]*q[3] - 2*q[0]*q[1]],
                    [2*q[1]*q[3] - 2*q[0]*q[2],
                     2*q[2]*q[3] + 2*q[0]*q[1],
                     q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3]]])
  return rotMat

#  def choice(N, weights=None):
#    """
#    Weighted random number selection of N values 
#    """
#    weights /= np.sum(weights)
#    cumweights = np.cumsum(weights)
#    beginIndicies = np.zeros(N, dtype=int)
#    randVals = np.random.uniform(size=N)
#
#    for (bI_idx,randVal) in zip(range(N),randVals):
#      beginIndicies[bI_idx] = sum(cumweights < randVal)
#    beginIndicies.sort()
#    return beginIndicies

def merge_dictionaries(dicts, required_consistency=[]):
  """
  Merges a list of dictionaries, giving priority to items in descending order.
  Items in the required_consistency list must be consistent with one another.
  """
  merged = {}
  for a in range(len(dicts)): # Loop over all dictionaries,
                              # giving priority to the first
    if not isinstance(dicts[a],dict):
      continue
    for key in dicts[a].keys():
      if isinstance(dicts[a][key],dict): # item is a dictionary
        merged[key] = merge_dictionaries(
          [dicts[n][key] for n in range(len(dicts)) if key in dicts[n].keys()],
          required_consistency=required_consistency)
      elif (key not in merged.keys()):
        # Merged dictionary will contain value from
        # first dictionary where key appears
        merged[key] = dicts[a][key]
        # Check for consistency with other dictionaries
        for b in (range(a) + range(a+1,len(dicts))):
          if isinstance(dicts[b],dict) and (key in dicts[b].keys()) \
              and (dicts[a][key] is not None) and (dicts[b][key] is not None):
            if (isinstance(dicts[b][key],np.ndarray)):
              inconsistent_items = (dicts[b][key]!=dicts[a][key]).any()
            else:
              inconsistent_items = (dicts[b][key]!=dicts[a][key])
            if inconsistent_items:
              if key in required_consistency:
                print 'Dictionary items for %s are inconsistent:'%key
                print dicts[a][key]
                print dicts[b][key]
                raise Exception('Items must be consistent!')
      elif (merged[key] is None): # Replace None
        merged[key] = dicts[a][key]
  return merged

def convert_dictionary_relpath(d, relpath_o=None, relpath_n=None):
  """
  Converts all file names in a dictionary from one relative path to another.
  If relpath_o is None, nothing is joined to the original path.
  If relpath_n is None, an absolute path is used.
  """
  converted = {}
  for key in d.keys():
    if d[key] is None:
      pass
    elif isinstance(d[key],dict):
      converted[key] = convert_dictionary_relpath(d[key],
        relpath_o = relpath_o, relpath_n = relpath_n)
    elif isinstance(d[key],str):
      if relpath_o is not None:
        p = abspath(join(relpath_o,d[key]))
      else:
        p = abspath(d[key])
      if exists(p): # Only save file names for existent paths
        if relpath_n is not None:
          converted[key] = os.path.relpath(p,relpath_n)
        else:
          converted[key] = p
  return converted

def HMStime(s):
  """
  Given the time in seconds, an appropriately formatted string.
  """
  if s<60.:
    return '%.3f s'%s
  elif s<3600.:
    return '%d:%.3f'%(int(s/60%60),s%60)
  else:
    return '%d:%d:%.3f'%(int(s/3600),int(s/60%60),s%60)

class NullDevice():
  """
  A device to suppress output
  """
  def write(self, s):
    pass

##############
# Main Class #
##############

class BindingPMF:
  def __init__(self,
      # Arguments are either constant or run-dependent.
      # Constant arguments are stored in the dir_cool and dir_dock directories
      # when they are first used.  In subsequent instances of Dock, the stored
      # values will be used.
      # Run-dependent arguments can vary from instance to instance.
      #
      # Directory and file names
      #   Run-dependent
      dir_dock=None,
      dir_cool=None,
      namd=None,
      vmd=None,
      sander=None,
      #   Stored in both dir_cool and dir_dock
      ligand_database=None,
      forcefield=None, frcmodList=None,
      ligand_prmtop=None, ligand_inpcrd=None,
      #   Stored in dir_dock
      receptor_prmtop=None, receptor_inpcrd=None,
      receptor_fixed_atoms=None,
      complex_prmtop=None, complex_inpcrd=None,
      complex_fixed_atoms=None,
      dir_grid=None, grid_LJr=None, grid_LJa=None, grid_ELE=None,
      # Arguments - simulation parameters and constants
      #   Run-dependent
      cool_repX_cycles=None,
      dock_repX_cycles=None,
      run_type=None,
      nodes=None,
      #   Defaults for dir_cool and dir_dock
      protocol=None, no_protocol_refinement=None, therm_speed=None,
      sampler=None, MCMC_moves=1,
      seeds_per_state=None, steps_per_seed=None,
      repX_cycles=None,
      min_repX_acc=None,
      sweeps_per_cycle=None, steps_per_sweep=None,
      keep_intermediate=None,
      phases=None,
      #   Stored in dir_cool
      cool_protocol=None, cool_therm_speed=None,
      cool_sampler=None,
      cool_seeds_per_state=None, cool_steps_per_seed=None,
      cool_sweeps_per_cycle=None, cool_steps_per_sweep=None,
      cool_keep_intermediate=None,
      #   Stored in dir_dock
      dock_protocol=None, dock_therm_speed=None,
      dock_sampler=None,
      dock_seeds_per_state=None, dock_steps_per_seed=None,
      dock_sweeps_per_cycle=None, dock_steps_per_sweep=None,
      dock_keep_intermediate=None,
      score=None, rmsd=None,
      site=None, site_center=None, site_direction=None, # Site parameters
      site_max_X=None, site_max_R=None, # Site parameters
      site_density=None,
      do_calc_random_dock_stats=None,
      receptor_Gas=None,
      receptor_GBSA=None,
      receptor_NAMD_Gas=None,
      receptor_NAMD_GBSA=None): # Energy values
    """Parses the input arguments and runs the requested docking calculation"""
    
    mod_path = join(os.path.dirname(a.__file__),'HREX.py')
    print """
###########
# AlGDock #
###########
Molecular docking with adaptively scaled alchemical interaction grids

version {0}
in {1}
last modified {2}

available CPUs: {3}
    """.format(a.__version__, mod_path, \
      time.ctime(os.path.getmtime(mod_path)), multiprocessing.cpu_count())
      
    self._do_multiprocess = (multiprocessing.cpu_count()>1) and \
      ((nodes is None) or ((nodes is not None) and (nodes>1)))
    
    self.confs = {'cool':{}, 'dock':{}}
    
    print '*** Directories ***'
    self.dir = {}
    self.dir['start'] = os.getcwd()
    
    if dir_dock is not None:
      self.dir['dock'] = abspath(dir_dock)
    else:
      self.dir['dock'] = abspath('.')
    
    if dir_cool is not None:
      self.dir['cool'] = abspath(dir_cool)
    else:
      self.dir['cool'] = self.dir['dock'] # Default that may be
                                          # overwritten by stored directory
    
    # Load previously stored file names and arguments
    FNs = {}
    args = {}
    for p in ['dock','cool']:
      params = self._load_progress(p)
      if params is not None:
        (fn_dict, arg_dict) = params
        FNs[p] = convert_dictionary_relpath(fn_dict,
          relpath_o=self.dir[p], relpath_n=None)
        args[p] = arg_dict
        if (p=='dock') and (dir_cool is None) and \
           ('dir_cool' in FNs[p].keys()) and \
           (FNs[p]['dir_cool'] is not None):
          self.dir['cool'] = FNs[p]['dir_cool']
      else:
        FNs[p] = {}
        args[p] = {}
  
    print self.dir
  
    # Set up file name dictionary
    print '\n*** Files ***'

    for p in ['cool','dock']:
      if p in FNs.keys():
        if FNs[p]!={}:
          print 'previous stored in %s directory:'%p
          print FNs[p]

    if dir_grid is None:
      dir_grid = ''
    FNs['new'] = {
      'ligand_database':a.findPath([ligand_database]),
      'forcefield':a.findPath([forcefield,'../Data/gaff.dat'] + \
                               a.search_paths['gaff.dat']),
      'frcmodList':frcmodList,
      'prmtop':{'L':a.findPath([ligand_prmtop]),
                'R':a.findPath([receptor_prmtop]),
                'RL':a.findPath([complex_prmtop])},
      'inpcrd':{'L':a.findPath([ligand_inpcrd]),
                'R':a.findPath([receptor_inpcrd]),
                'RL':a.findPath([complex_inpcrd])},
      'fixed_atoms':{'R':a.findPath([receptor_fixed_atoms]),
                     'RL':a.findPath([complex_fixed_atoms])},
      'grids':{'LJr':a.findPath([grid_LJr,
                               join(dir_grid,'LJr.nc'),
                               join(dir_grid,'LJr.dx'),
                               join(dir_grid,'LJr.dx.gz')]),
               'LJa':a.findPath([grid_LJa,
                               join(dir_grid,'LJa.nc'),
                               join(dir_grid,'LJa.dx'),
                               join(dir_grid,'LJa.dx.gz')]),
               'ELE':a.findPath([grid_ELE,
                               join(dir_grid,'electrostatic.nc'),
                               join(dir_grid,'electrostatic.dx'),
                               join(dir_grid,'electrostatic.dx.gz'),
                               join(dir_grid,'pbsa.nc'),
                               join(dir_grid,'pbsa.dx'),
                               join(dir_grid,'pbsa.dx.gz')])},
      'dir_cool':self.dir['cool'],
      'namd':a.findPath([namd] + a.search_paths['namd']),
      'vmd':a.findPath([vmd] + a.search_paths['vmd']),
      'sander':a.findPath([vmd] + a.search_paths['sander'])}

    if not (FNs['cool']=={} and FNs['dock']=={}):
      print 'from arguments and defaults:'
      print FNs['new']
      print '\nto be used:'

    self._FNs = merge_dictionaries(
      [FNs[src] for src in ['new','cool','dock']],
      required_consistency=['L','R','RL','ligand_database'])
  
    # Default: a force field modification is in the same directory as the ligand
    if (self._FNs['frcmodList'] is None):
      dir_lig = os.path.dirname(self._FNs['prmtop']['L'])
      frcmod = a.findPath([\
        abspath(join(dir_lig, \
          os.path.basename(self._FNs['prmtop']['L'])[:-7]+'.frcmod')),\
        abspath(join(dir_lig,'lig.frcmod')),\
        abspath(join(dir_lig,'ligand.frcmod'))])
      if frcmod is not None:
        self._FNs['frcmodList'] = [frcmod]
    elif isinstance(self._FNs['frcmodList'],str):
      self._FNs['frcmodList'] = [self._FNs['frcmodList']]

    print self._FNs
    
    args['default_cool'] = {
        'protocol':'Adaptive',
        'no_protocol_refinement':False,
        'therm_speed':1.0,
        'sampler':'HMC',
        'seeds_per_state':50,
        'steps_per_seed':20000,
        'repX_cycles':20,
        'min_repX_acc':0.3,
        'sweeps_per_cycle':100,
        'steps_per_sweep':2500,
        'phases':['NAMD_Gas','NAMD_GBSA'],
        'keep_intermediate':False}

    args['default_dock'] = dict(args['default_cool'].items() + {
      'site':None, 'site_center':None, 'site_direction':None,
      'site_max_X':None, 'site_max_R':None,
      'site_density':50.,
      'MCMC_moves':1,
      'do_calc_random_dock_stats':False,
      'rmsd':False,
      'score':False}.items() + \
      [('receptor_'+phase,None) for phase in allowed_phases])

    # Store passed arguments in dictionary
    namespace = locals()
    for p in ['cool','dock']:
      args['new_'+p] = {}
      for key in args['default_'+p].keys():
        specific_key = p + '_' + key
        if (specific_key in namespace.keys()) and \
           (namespace[specific_key] is not None):
          # Use the specific key if it exists
          args['new_'+p][key] = namespace[specific_key]
        elif (key in ['site_center', 'site_direction'] +
                     ['receptor_'+phase for phase in allowed_phases]) and \
             (namespace[key] is not None):
          # Convert these to arrays of floats
          args['new_'+p][key] = np.array(namespace[key], dtype=float)
        else:
          # Use the general key
          args['new_'+p][key] = namespace[key]

    self.params = {}
    for p in ['cool','dock']:
      self.params[p] = merge_dictionaries(
        [args[src] for src in ['new_'+p,'default_'+p,p]])

    # Allow updates of receptor energies
    for phase in allowed_phases:
      if args['new_dock']['receptor_'+phase] is not None:
        self.params['dock']['receptor_'+phase] = args['new_dock']['receptor_'+phase]

    # Allow updates of score parameter
    if args['new_dock']['score'] is not None:
      self.params['dock']['score'] = args['new_dock']['score']

    print '\n*** Simulation parameters and constants ***'
    for p in ['cool','dock']:
      print 'for %s:'%p
      print self.params[p]

    # Variables dependent on the parameters
    self.original_Es = [[{}]]
    for phase in allowed_phases:
      if self.params['dock']['receptor_'+phase] is not None:
        self.original_Es[0][0]['R'+phase] = \
          self.params['dock']['receptor_'+phase]
      else:
        self.original_Es[0][0]['R'+phase] = None
        
    self._scalables = ['sLJr','sELE','LJr','LJa','ELE']

    self.lambda_full = {'T':T_HIGH,'MM':True,'site':True}
    self.lambda_scalables = {'T':T_HIGH}
    for scalable in self._scalables:
        self.lambda_full[scalable] = 1
        self.lambda_scalables[scalable] = 1

    # Check for existence of required files
    do_dock = (hasattr(args,'run_type') and \
              (args.run_type in ['pose_energies','random_dock',\
                'initial_dock','dock','all']))
              
    for FN in [self._FNs['ligand_database'],
               self._FNs['forcefield'],
               self._FNs['prmtop']['L'],
               self._FNs['inpcrd']['L']]:
      if (FN is None) or (not isfile(FN)):
        raise Exception('Required file %s is missing!'%FN)

    for FN in [self._FNs['prmtop']['RL'],
               self._FNs['inpcrd']['RL'],
               self._FNs['fixed_atoms']['RL'],
               self._FNs['grids']['LJr'],
               self._FNs['grids']['LJa'],
               self._FNs['grids']['ELE']]:
      if (FN is None) or (not isfile(FN)):
        if do_dock:
          raise Exception('Missing file is required for docking!')
        else:
          print 'Missing file is required for docking!'

    if ((self._FNs['inpcrd']['RL'] is None) and \
        (self._FNs['inpcrd']['R'] is None)):
        if do_dock:
          raise Exception('Receptor coordinates needed for docking')
        else:
          print 'Receptor coordinates needed for docking'
        
    print '\n*** Setting up the simulation ***'
    self._setup_universe(do_dock = do_dock)

    self.run_type = run_type
    self.cycles_run = 0
    if run_type=='pose_energies':
      (unique_confs, unique_conf_Es) = self._get_dock_score_confs()
      if self.params['dock']['score'] is True:
        prefix = 'xtal'
      else:
        prefix = os.path.basename(self.params['dock']['score']).split('.')[0]
      E = self._calc_E(confs=unique_confs, E=unique_conf_Es, type='all', \
        prefix=prefix)
      self._write_pkl_gz(join(self.dir['dock'],prefix+'.pkl.gz'),E)
    elif run_type=='one_step':
      # Does one of the following:
      # 1. Initial cooling
      # 2. Cooling replica exchange
      # 3. Initial docking
      # 4. One cycle of docking replica exchange
      self.cool()
      if self.cycles_run==0:
        self.dock()
    elif run_type=='initial_cool':
      self.initial_cool()
    elif run_type=='cool':
      self.cool()
    elif run_type=='random_dock':
      self.initial_dock(randomOnly=True)
    elif run_type=='initial_dock':
      self.initial_dock()
    elif run_type=='dock':
      self.dock()
    elif run_type=='all':
      self.cool()
      self.dock()
    elif run_type=='store_params':
      self._save_progress('cool')
      self._save_progress('dock')
    elif run_type=='free_energies':
      self.calc_f_L()
      self.calc_f_RL()
    elif run_type=='postprocess':
      self._postprocess()
      self.calc_f_L()
      self.calc_f_RL()
    elif run_type=='redo_postprocess':
      self._postprocess(redo_dock=True)
      self.calc_f_L()
      self.calc_f_RL()
    elif run_type=='clear_intermediates':
      for process in ['cool','dock']:
        print 'Clearing intermediates for '+process
        for state_ind in range(1,len(self.confs[process]['samples'])-1):
          for cycle_ind in range(len(self.confs[process]['samples'][state_ind])):
            self.confs[process]['samples'][state_ind][cycle_ind] = []
        self._save_progress(process)

  def _setup_universe(self, do_dock=True):
    """Creates an MMTK InfiniteUniverse and adds the ligand"""

    MMTK.Database.molecule_types.directory = \
      os.path.dirname(self._FNs['ligand_database'])

    # Force fields
    from MMTK.ForceFields import Amber12SBForceField

    self._forceFields = {}
    self._forceFields['gaff'] = Amber12SBForceField(
      parameter_file=self._FNs['forcefield'],mod_files=self._FNs['frcmodList'])

    if (self.params['dock']['site']=='Sphere') and \
       (self.params['dock']['site_center'] is not None) and \
       (self.params['dock']['site_max_R'] is not None):
      from AlGDock.ForceFields.Sphere.Sphere import SphereForceField
      self._forceFields['site'] = SphereForceField(
        center=self.params['dock']['site_center'],
        max_R=self.params['dock']['site_max_R'], name='site')
    elif (self.params['dock']['site']=='Cylinder') and \
         (self.params['dock']['site_center'] is not None) and \
         (self.params['dock']['site_direction'] is not None):
      from AlGDock.ForceFields.Cylinder.Cylinder import CylinderForceField
      self._forceFields['site'] = CylinderForceField(
        origin=self.params['dock']['site_center'],
        direction=self.params['dock']['site_direction'],
        max_X=self.params['dock']['site_max_X'],
        max_R=self.params['dock']['site_max_R'], name='site')
    else:
      if do_dock:
        raise Exception('Binding site type not recognized!')
      else:
        print 'Binding site type not recognized!'
  
    # Set up the system
    import sys
    original_stderr = sys.stderr
    sys.stderr = NullDevice()

    self.molecule = MMTK.Molecule(\
      os.path.basename(self._FNs['ligand_database']))

    sys.stderr = original_stderr

    # Helpful variables for referencing and indexing atoms in the molecule
    self.molecule.heavy_atoms = [ind for (atm,ind) in zip(self.molecule.atoms,range(self.molecule.numberOfAtoms())) if atm.type.name!='hydrogen']
    self.molecule.nhatoms = len(self.molecule.heavy_atoms)

    self.molecule.prmtop_atom_order = np.array([atom.number \
      for atom in self.molecule.prmtop_order], dtype=int)
    self.molecule.inv_prmtop_atom_order = np.zeros(shape=len(self.molecule.prmtop_atom_order), dtype=int)
    for i in range(len(self.molecule.prmtop_atom_order)):
      self.molecule.inv_prmtop_atom_order[self.molecule.prmtop_atom_order[i]] = i
    
    # Randomly rotate the molecule and translate it into the binding site
    from Scientific.Geometry.Transformation import Rotation
    self.molecule.applyTransformation(Rotation(random_rotate()))
    if 'site' in self._forceFields.keys():
      self.molecule.translateTo(Vector(self._forceFields['site'].randomPoint()))
    else:
      print 'Molecule not translated into binding site'

    self.universe = MMTK.Universe.InfiniteUniverse()
    self.universe.addObject(self.molecule)
    self._set_universe_force_field({'MM':True, 'T':T_HIGH})
    self._ligand_natoms = self.universe.numberOfAtoms()

    from Integrators.VelocityVerlet.VelocityVerlet \
      import VelocityVerletIntegrator
    from Integrators.HamiltonianMonteCarlo.HamiltonianMonteCarlo \
      import HamiltonianMonteCarloIntegrator
    from NUTS import NUTSIntegrator  # @UnresolvedImport

    # Samplers may accept the following options:
    # steps - number of MD steps
    # T - temperature in K
    # delta_t - MD time step
    # normalize - normalizes configurations
    # adapt - uses an adaptive time step
    self.sampler = {}
    self.sampler['init'] = NUTSIntegrator(self.universe)
    for p in ['cool', 'dock']:
      if self.params[p]['sampler'] == 'HMC':
        self.sampler[p] = HamiltonianMonteCarloIntegrator(self.universe)
      elif self.params[p]['sampler'] == 'NUTS':
        self.sampler[p] = NUTSIntegrator(self.universe)
      elif self.params[p]['sampler'] == 'VV':
        self.sampler[p] = VelocityVerletIntegrator(self.universe)
      else:
        raise Exception('Unrecognized sampler!')

    # Determine ligand atomic index
    if (self._FNs['prmtop']['R'] is not None) and \
       (self._FNs['prmtop']['RL'] is not None):
      import AlGDock.IO
      IO_prmtop = AlGDock.IO.prmtop()
      prmtop_R = IO_prmtop.read(self._FNs['prmtop']['R'])
      prmtop_RL = IO_prmtop.read(self._FNs['prmtop']['RL'])
      ligand_ind = [ind for (r,ind) in \
        zip(prmtop_RL['RESIDUE_LABEL'],range(len(prmtop_RL['RESIDUE_LABEL']))) \
        if r not in prmtop_R['RESIDUE_LABEL']]
      if len(ligand_ind)==0:
        raise Exception('Ligand not found in complex prmtop')
      elif len(ligand_ind) > 1:
        raise Exception('Ligand residue label is ambiguous')
      self._ligand_first_atom = prmtop_RL['RESIDUE_POINTER'][ligand_ind[0]] - 1
    else:
      if do_dock:
        raise Exception('Missing AMBER prmtop files for receptor')
      else:
        print 'Missing AMBER prmtop files for receptor'

    # Read the reference ligand and receptor coordinates
    import AlGDock.IO
    IO_crd = AlGDock.IO.crd()
    if self._FNs['inpcrd']['R'] is not None:
      if isfile(self._FNs['inpcrd']['L']):
        lig_crd = IO_crd.read(self._FNs['inpcrd']['L'], multiplier=0.1)
      self.confs['receptor'] = IO_crd.read(self._FNs['inpcrd']['R'], multiplier=0.1)
    elif self._FNs['inpcrd']['RL'] is not None:
      complex_crd = IO_crd.read(self._FNs['inpcrd']['RL'], multiplier=0.1)
      lig_crd = complex_crd[self._ligand_first_atom:self._ligand_first_atom + self._ligand_natoms,:]
      self.confs['receptor'] = np.vstack(\
        (complex_crd[:self._ligand_first_atom,:],\
         complex_crd[self._ligand_first_atom + self._ligand_natoms:,:]))
    elif self._FNs['inpcrd']['L'] is not None:
      if isfile(self._FNs['inpcrd']['L']):
        lig_crd = IO_crd.read(self._FNs['inpcrd']['L'], multiplier=0.1)
    self.confs['ligand'] = lig_crd[self.molecule.inv_prmtop_atom_order,:]
    
    if self.params['dock']['rmsd'] is not False:
      if self.params['dock']['rmsd'] is True:
        rmsd_crd = self.confs['ligand']
      else:
        rmsd_crd = IO_crd.read(self.params['dock']['rmsd'], \
          natoms=self.universe.numberOfAtoms(), multiplier=0.1)
        rmsd_crd = rmsd_crd[self.molecule.inv_prmtop_atom_order,:]
      self.confs['rmsd'] = rmsd_crd[self.molecule.heavy_atoms,:]
    
    # Load progress
    self._postprocess(readOnly=True)
    self.calc_f_L(readOnly=True)
    self.calc_f_RL(readOnly=True)

  ###########
  # Cooling #
  ###########
  def initial_cool(self):
    """
    Cools the ligand from T_HIGH to T_TARGET
    
    Intermediate thermodynamic states are chosen such that
    thermodynamic length intervals are approximately constant.
    Configurations from each state are subsampled to seed the next simulation.
    """

    if (len(self.cool_protocol)>0) and (self.cool_protocol[-1]['crossed']):
      return # Initial cooling is already complete
    
    self._set_lock('cool')
    self.tee("\n>>> Initial cooling of the ligand "+\
             "from %d K to %d K"%(T_HIGH,T_TARGET))
    cool_start_time = time.time()

    # Set up the force field
    self.cool_protocol = [{'MM':True, 'T':T_HIGH, \
                          'delta_t':1.5*MMTK.Units.fs,
                          'a':0.0, 'crossed':False}]
    self._set_universe_force_field(self.cool_protocol[-1])

    # Minimize and ramp the temperature from 0 to the desired high temperature
    from MMTK.Minimization import ConjugateGradientMinimizer # @UnresolvedImport
    minimizer = ConjugateGradientMinimizer(self.universe)
    for rep in range(50):
      x_o = np.copy(self.universe.configuration().array)
      e_o = self.universe.energy()
      minimizer(steps = 10)
      e_n = self.universe.energy()
      if np.isnan(e_n) or (e_o-e_n)<1000:
        self.universe.configuration().array = x_o
        break

    for T in np.linspace(0.,T_HIGH,33)[1:]:
      self.sampler['init'](steps = 500, T=T,\
                           delta_t=self.delta_t, steps_per_trial = 100)
    self.universe.normalizePosition()

    # Run at high temperature
    state_start_time = time.time()
    conf = self.universe.configuration().array
    (confs, Es_MM, tau_ac) = self._initial_sim_state('cool', \
      [conf for n in range(self.params['cool']['seeds_per_state'])], \
      normalize=True)
    self.confs['cool']['replicas'] = [confs[np.argmin(Es_MM)]]
    self.confs['cool']['samples'] = [[confs]]
    self.cool_Es = [[{'MM':Es_MM}]]
    self.cool_protocol[-1]['delta_t'] = self.delta_t
    self.cool_protocol[-1]['tau_ac'] = tau_ac
    self.tee("  generated %d configurations "%len(confs) + \
             "(tau_ac=%f, "%self.cool_protocol[-1]['tau_ac'] + \
             "dt=%f ps) "%self.cool_protocol[-1]['delta_t'] + \
             "at %d K "%self.cool_protocol[-1]['T'] + \
             "in " + HMStime(time.time()-state_start_time))
    tL_tensor = Es_MM.std()/(R*T_HIGH*T_HIGH)

    # Main loop for initial cooling:
    # choose new temperature, randomly select seeds, simulate
    T = T_HIGH
    while (not self.cool_protocol[-1]['crossed']):
      # Choose new temperature
      To = self.cool_protocol[-1]['T']
      crossed = False
      if tL_tensor>0:
        T = To - self.params['cool']['therm_speed']/tL_tensor
        if T < T_TARGET:
          T = T_TARGET
          crossed = True
      else:
        raise Exception('No variance in configuration energies')
      lambda_n = {'T':T, 'a':(T_HIGH-T)/(T_HIGH-T_TARGET), 'MM':True,
                  'crossed':crossed}
      self.cool_protocol.append(lambda_n)
      self._set_universe_force_field(lambda_n)

      # Randomly select seeds for new trajectory
      w = np.exp(-Es_MM/R*(1/T-1/To))
      seedIndicies = np.random.choice(len(Es_MM),
        size = self.params['cool']['seeds_per_state'],
        p = w/sum(w))

      # Simulate and store data
      confs_o = confs
      Es_MM_o = Es_MM

      state_start_time = time.time()
      (confs, Es_MM, tau_ac) = self._initial_sim_state('cool', \
        [confs[ind] for ind in seedIndicies], normalize=True)
      self.cool_protocol[-1]['delta_t'] = self.delta_t
      self.cool_protocol[-1]['tau_ac'] = tau_ac
      self.tee("  generated %d configurations "%len(confs) + \
               "(tau_ac=%f, "%self.cool_protocol[-1]['tau_ac'] + \
               "dt=%f ps) "%self.cool_protocol[-1]['delta_t'] + \
               "at %d K "%self.cool_protocol[-1]['T'] + \
               "in " + (HMStime(time.time()-state_start_time)))

      # Estimate the mean replica exchange acceptance rate
      # between the previous and new state
      (u_kln,N_k) = self._u_kln([[{'MM':Es_MM_o}],[{'MM':Es_MM}]], \
                                self.cool_protocol[-2:])
      N = min(N_k)
      acc = np.exp(-u_kln[0,1,:N]-u_kln[1,0,:N]+u_kln[0,0,:N]+u_kln[1,1,:N])
      mean_acc = np.mean(np.minimum(acc,np.ones(acc.shape)))

      if mean_acc<self.params['cool']['min_repX_acc']:
        # If the acceptance probability is too low,
        # reject the state and restart
        self.cool_protocol.pop()
        confs = confs_o
        Es_MM = Es_MM_o
        tL_tensor = tL_tensor*1.25 # Use a smaller step
        self.tee("  rejected new state, as estimated replica exchange acceptance rate of %f is too low"%mean_acc)
      elif (mean_acc>0.95) and (not crossed):
        # If the acceptance probability is too low,
        # reject the previous state and restart
        self.confs['cool']['replicas'][-1] = confs[np.argmin(Es_MM)]
        self.cool_protocol.pop()
        self.cool_protocol[-1] = lambda_n
        tL_tensor = Es_MM.std()/(R*T*T) # Metric tensor for the thermodynamic length
        self.tee("  rejected previous state, as estimated replica exchange acceptance rate of %f is too high"%mean_acc)
      else:
        self.confs['cool']['replicas'].append(confs[np.argmin(Es_MM)])
        self.confs['cool']['samples'].append([confs])
        if len(self.confs['cool']['samples'])>2 and \
            (not self.params['cool']['keep_intermediate']):
          self.confs['cool']['samples'][-2] = []
        self.cool_Es.append([{'MM':Es_MM}])
        tL_tensor = Es_MM.std()/(R*T*T) # Metric tensor for the thermodynamic length
        self.tee("  estimated replica exchange acceptance rate is %f\n"%mean_acc)

    # Save data
    self._cool_cycle += 1
    self._cool_total_cycle += 1
    self._save_progress('cool')
    self.tee("\nElapsed time for initial cooling: " + \
      HMStime(time.time()-cool_start_time))
    self._clear_lock('cool')

  def cool(self):
    """
    Samples different ligand configurations 
    at thermodynamic states between T_HIGH and T_TARGET
    """
    self._sim_process('cool')

  def calc_f_L(self, readOnly=False, redo=False):
    """
    Calculates ligand-specific free energies:
    1. solvation free energy of the ligand using single-step 
       free energy perturbation
    2. reduced free energy of cooling the ligand from T_HIGH to T_TARGET
    """
    # Initialize variables as empty lists or by loading data
    f_L_FN = join(self.dir['cool'],'f_L.pkl.gz')
    if redo:
      if isfile(f_L_FN):
        os.remove(f_L_FN)
      dat = None
    else:
      dat = self._load_pkl_gz(f_L_FN)
    if (dat is not None):
      (self.f_L_solvation, self.f_cool, self.cool_equilibrated_cycle, \
       self.cool_mean_acc) = dat
    else:
      for var in ['f_L_solvation', 'f_cool', \
                  'cool_equilibrated_cycle', 'cool_mean_acc']:
        setattr(self,var,[])
    if readOnly:
      return
    
    if self.cool_protocol==[]:
      return # Initial cooling is incomplete
      
    K = len(self.cool_protocol)

    # Make sure all the energies are available
    for c in range(self._cool_cycle):
      if len(self.cool_Es[-1][c].keys())==0:
        self.tee("  skipping the cooling free energy calculation")
        return

    # Estimate cycle at which simulation has equilibrated
    u_Ks = [self._u_kln([self.cool_Es[-1][c]],[self.cool_protocol[-1]]) \
      for c in range(self._cool_cycle)]
    mean_u_Ks = np.array([np.mean(u_K) for u_K in u_Ks])
    std_u_Ks = np.array([np.std(u_K) for u_K in u_Ks])
    for c in range(len(self.cool_equilibrated_cycle), self._cool_cycle):
      nearMean = (mean_u_Ks - mean_u_Ks[c])<std_u_Ks[c]
      for (item, ind) in zip(nearMean, range(len(nearMean))):
        if item:
          self.cool_equilibrated_cycle.append(ind)
          break

    updated = False
    
    # Calculate solvation free energies that have not already been calculated,
    # in units of RT
    for c in range(len(self.f_L_solvation), self._cool_cycle):
      fromCycle = self.cool_equilibrated_cycle[c]
      toCycle = c + 1
      
      # Arbitrarily, solvation is the
      # 'forward' direction and desolvation the 'reverse'
      u_F = {}
      for key in ['L'+phase for phase in self.params['cool']['phases']]+['MM']:
        u_F[key] = np.concatenate(\
          [self.cool_Es[-1][n][key] for n in range(fromCycle,toCycle)])

      # Single step FEP
      f_L_solvation = {}
      for phase in self.params['cool']['phases']:
        du_F = (u_F['L'+phase][:,-1] - u_F['MM'])/RT_TARGET
        min_du_F = min(du_F)
        f_L_solvation[phase] = -np.log(np.exp(-du_F+min_du_F).mean()) + min_du_F

      self.f_L_solvation.append(f_L_solvation)
      self.tee("  calculated %s solvation free energy of %f RT "%(\
                  self.params['cool']['phases'][-1],
                  self.f_L_solvation[-1][self.params['cool']['phases'][-1]])+\
               "using cycles %d to %d"%(fromCycle, toCycle-1))
      updated = True

    # Calculate cooling free energies that have not already been calculated,
    # in units of RT
    for c in range(len(self.f_cool), self._cool_cycle):
      fromCycle = self.cool_equilibrated_cycle[c]
      toCycle = c + 1

      # Cooling free energy
      cool_Es = []
      for cool_Es_state in self.cool_Es:
        cool_Es.append(cool_Es_state[fromCycle:toCycle])
      (u_kln,N_k) = self._u_kln(cool_Es,self.cool_protocol)
      self.f_cool.append(self._run_MBAR(u_kln,N_k))

      # Average acceptance probabilities
      cool_mean_acc = np.zeros(K-1)
      for k in range(0, K-1):
        (u_kln,N_k) = self._u_kln(cool_Es[k:k+2],self.cool_protocol[k:k+2])
        N = min(N_k)
        acc = np.exp(-u_kln[0,1,:N]-u_kln[1,0,:N]+u_kln[0,0,:N]+u_kln[1,1,:N])
        cool_mean_acc[k] = np.mean(np.minimum(acc,np.ones(acc.shape)))
      self.cool_mean_acc.append(cool_mean_acc)

      self.tee("  calculated cooling free energy of %f RT "%(\
                  self.f_cool[-1][-1])+\
               "using cycles %d to %d"%(fromCycle, toCycle-1))
      updated = True

    if updated:
      self._write_pkl_gz(f_L_FN, \
        (self.f_L_solvation, self.f_cool, \
         self.cool_equilibrated_cycle, self.cool_mean_acc))

  ###########
  # Docking #
  ###########
  def random_dock(self):
    # Select samples from the first cooling state and make sure there are enough
    E_MM = []
    confs = []
    for k in range(len(self.cool_Es[0])):
      E_MM += list(self.cool_Es[0][k]['MM'])
      confs += list(self.confs['cool']['samples'][0][k])
    while len(E_MM)<self.params['dock']['seeds_per_state']:
      self.tee("More samples from high temperature ligand simulation needed")
      self._replica_exchange('cool')
      E_MM = []
      confs = []
      for k in range(len(self.cool_Es[0])):
        E_MM += list(self.cool_Es[0][k]['MM'])
        confs += list(self.confs['cool']['samples'][0][k])

    random_dock_inds = np.array(np.linspace(0,len(E_MM), \
      self.params['dock']['seeds_per_state'],endpoint=False),dtype=int)
    cool0_Es_MM = [E_MM[ind]  for ind in random_dock_inds]
    cool0_confs = [confs[ind] for ind in random_dock_inds]

    # Do the random docking
    self.tee("\n>>> Initial docking")

    # Set up the force field with full interaction grids
    self._set_universe_force_field(self.lambda_scalables)
  
    lambda_o = {'T':T_HIGH, 'MM':True, 'site':True, \
                'crossed':False, 'a':0.0}
    for scalable in self._scalables:
      lambda_o[scalable] = 0
    self.dock_protocol = [lambda_o]

    # Either loads or generates the random translations and rotations for the first state of docking
    if not (hasattr(self,'_random_trans') and hasattr(self,'_random_rotT')):
      self._max_n_trans = 10000
      # Default density of points is 50 per nm**3
      self._n_trans = max(min(np.int(np.ceil(self._forceFields['site'].volume*self.params['dock']['site_density'])),self._max_n_trans),5)
      self._random_trans = np.ndarray((self._max_n_trans), dtype=Vector)
      for ind in range(self._max_n_trans):
        self._random_trans[ind] = Vector(self._forceFields['site'].randomPoint())
      self._max_n_rot = 100
      self._n_rot = 100
      self._random_rotT = np.ndarray((self._max_n_rot,3,3))
      for ind in range(self._max_n_rot):
        self._random_rotT[ind,:,:] = np.transpose(random_rotate())
    else:
      self._max_n_trans = self._random_trans.shape[0]
      self._n_rot = self._random_rotT.shape[0]

    # Get interaction energies.
    # Loop over configurations, random rotations, and random translations
    E = {}
    for term in (['MM','site']+self._scalables):
      # Large array creation may cause MemoryError
      E[term] = np.zeros((self.params['dock']['seeds_per_state'], \
        self._max_n_rot,self._n_trans))
    self.tee("  allocated memory for interaction energies")

    converged = False
    n_trans_o = 0
    n_trans_n = self._n_trans
    while not converged:
      for c in range(self.params['dock']['seeds_per_state']):
        E['MM'][c,:,:] = cool0_Es_MM[c]
        for i_rot in range(self._n_rot):
          conf_rot = Configuration(self.universe,\
            np.dot(cool0_confs[c], self._random_rotT[i_rot,:,:]))
          for i_trans in range(n_trans_o, n_trans_n):
            self.universe.setConfiguration(conf_rot)
            self.universe.translateBy(self._random_trans[i_trans])
            eT = self.universe.energyTerms()
            for (key,value) in eT.iteritems():
              E[term_map[key]][c,i_rot,i_trans] += value
      E_c = {}
      for term in E.keys():
        # Large array creation may cause MemoryError
        E_c[term] = np.ravel(E[term][:,:self._n_rot,:n_trans_n])
      self.tee("  allocated memory for %d translations"%n_trans_n)
      (u_kln,N_k) = self._u_kln([E_c],\
        [lambda_o,self._next_dock_state(E=E_c, lambda_o=lambda_o)])
      du = u_kln[0,1,:] - u_kln[0,0,:]
      bootstrap_reps = 50
      f_grid0 = np.zeros(bootstrap_reps)
      for b in range(bootstrap_reps):
        du_b = du[np.random.randint(0, len(du), len(du))]
        f_grid0[b] = -np.log(np.exp(-du_b+min(du_b)).mean()) + min(du_b)
      f_grid0_std = f_grid0.std()
      converged = f_grid0_std<0.1
      if not converged:
        self.tee("  with %s translations "%n_trans_n + \
                 "the predicted free energy difference is %f (%f)"%(\
                 f_grid0.mean(),f_grid0_std))
        if n_trans_n == self._max_n_trans:
          break
        n_trans_o = n_trans_n
        n_trans_n = min(n_trans_n + 25, self._max_n_trans)
        for term in (['MM','site']+self._scalables):
          # Large array creation may cause MemoryError
          # Has not been tested
          E[term] = np.dstack(E[term], \
            np.zeros((self.params['dock']['seeds_per_state'], \
            self._max_n_rot,25)))

    if self._n_trans != n_trans_n:
      self._n_trans = n_trans_n
      
    self.tee("  %d ligand configurations "%len(cool0_Es_MM) + \
             "were randomly docked into the binding site using "+ \
             "%d translations and %d rotations "%(n_trans_n,self._n_rot))
    self.tee("  the predicted free energy difference between the" + \
             " first and second docking states is " + \
             "%f (%f)"%(f_grid0.mean(),f_grid0_std))

    ravel_start_time = time.time()
    for term in E.keys():
      E[term] = np.ravel(E[term][:,:self._n_rot,:self._n_trans])
    self.tee("  raveled energy terms in " + \
      HMStime(time.time()-ravel_start_time))

    stats_start_time = time.time()
    if self.params['dock']['do_calc_random_dock_stats']:
      self.calc_random_dock_stats(E)
      self.tee("  calculated random dock stats in " + \
        HMStime(time.time()-stats_start_time))

    return (cool0_confs, E)

  def initial_dock(self, randomOnly=False):
    """
      Docks the ligand into the receptor
      
      The first state of docking is sampled by randomly placing configurations 
      from the high temperature ligand simulation into the binding site.
      Intermediate thermodynamic states are chosen such that
      thermodynamic length intervals are approximately constant.
      Configurations from each state are subsampled to seed the next simulation.
    """
    
    if (len(self.dock_protocol)>0) and (self.dock_protocol[-1]['crossed']):
      return # Initial docking already complete

    self._set_lock('dock')
    dock_start_time = time.time()

    if self.dock_protocol==[]:
      (cool0_confs, E) = self.random_dock()
      self.tee("  random docking complete in " + \
               HMStime(time.time()-dock_start_time))
      if randomOnly:
        self._clear_lock('dock')
        return
    else:
      # Continuing from a previous docking instance
      self.tee("\n>>> Initial docking, continued")
      confs = self.confs['dock']['samples'][-1][0]
      E = self.dock_Es[-1][0]

    lambda_o = self.dock_protocol[-1]

    # Main loop for initial docking:
    # choose new thermodynamic variables,
    # randomly select seeds,
    # simulate
    rejectStage = 0
    while (not self.dock_protocol[-1]['crossed']):
      # Determine next value of the protocol
      lambda_n = self._next_dock_state(E = E, lambda_o = lambda_o, \
          pow = rejectStage)
      self.dock_protocol.append(lambda_n)
      if len(self.dock_protocol)>100:
        raise Exception('Too many replicas!')

      # Randomly select seeds for new trajectory
      u_o = self._u_kln([E],[lambda_o])
      u_n = self._u_kln([E],[lambda_n])
      du = u_n-u_o
      w = np.exp(-du+min(du))
      seedIndicies = np.random.choice(len(u_o), \
        size = self.params['dock']['seeds_per_state'], \
        p=w/sum(w))

      if len(self.dock_protocol)==2: # Cooling state 0 configurations, randomly oriented
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
        seeds = [confs[ind] for ind in seedIndicies]
      self.confs['dock']['seeds'] = seeds

      # Store old data
      confs_o = confs
      E_o = E

      # Simulate
      sim_start_time = time.time()
      self._set_universe_force_field(lambda_n)
      (confs, Es_tot, tau_ac) = self._initial_sim_state('dock', seeds, \
        normalize=False, doMC=True)
      
      # Store data
      self.dock_protocol[-1]['delta_t'] = self.delta_t
      self.dock_protocol[-1]['tau_ac'] = tau_ac

      self.tee("  generated %d configurations "%len(confs) + \
               "(tau_ac=%f, "%self.dock_protocol[-1]['tau_ac'] + \
               "dt=%f ps) "%self.dock_protocol[-1]['delta_t'] + \
               "with progress %f "%self.dock_protocol[-1]['a'] + \
               "in " + HMStime(time.time()-sim_start_time))

      # Get state energies
      E = self._calc_E(confs)

      if len(self.dock_protocol)>2:
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
        elif (mean_acc>0.95) and (not lambda_n['crossed']):
          # If the acceptance probability is too high,
          # reject the previous state and restart
          self.confs['dock']['replicas'][-1] = confs[np.argmin(Es_tot)]
          self.dock_protocol.pop()
          self.dock_protocol[-1] = lambda_n
          self.dock_protocol[-1]['delta_t'] = self.delta_t
          rejectStage = 0
          lambda_o = lambda_n
          self.tee("  rejected previous state, as estimated replica exchange acceptance rate of %f is too high"%mean_acc)
        else:
          # Store data and continue with initialization
          self.confs['dock']['replicas'].append(confs[np.argmin(Es_tot)])
          self.confs['dock']['samples'].append([confs])
          if (not self.params['dock']['keep_intermediate']):
            self.confs['dock']['samples'][-2] = []
          self.dock_Es.append([E])
          self.dock_protocol[-1]['delta_t'] = self.delta_t
          rejectStage = 0
          lambda_o = lambda_n
          self.tee("  the estimated replica exchange acceptance rate is %f"%mean_acc)
      else:
        # Store data and continue with initialization (first time)
        self.confs['dock']['replicas'].append(confs[np.argmin(Es_tot)])
        self.confs['dock']['samples'].append([confs])
        self.dock_Es.append([E])
        self.dock_protocol[-1]['delta_t'] = self.delta_t
        rejectStage = 0
        lambda_o = lambda_n
      
      self._save_progress('dock')
      self.tee("")

    K = len(self.dock_protocol)
    self.tee("  %d states in the docking process sampled in %s"%(K,\
      HMStime(time.time()-dock_start_time)))
      
    self._dock_cycle += 1
    self._dock_total_cycle += 1
    self._save_progress('dock')
    self.tee("\nElapsed time for initial docking: " + \
      HMStime(time.time()-dock_start_time))
    self._clear_lock('dock')

  def dock(self):
    """
    Docks the ligand into the binding site
    by simulating at thermodynamic states
    between decoupled and fully interacting and
    between T_HIGH and T_TARGET
    """
    self._sim_process('dock')

  def calc_f_RL(self, readOnly=False, redo=False):
    """
    Calculates the binding potential of mean force
    """
    if self.dock_protocol==[]:
      return # Initial docking is incomplete

    # Initialize variables as empty lists or by loading data
    f_RL_FN = join(self.dir['dock'],'f_RL.pkl.gz')
    if redo:
      if isfile(f_RL_FN):
        os.remove(f_RL_FN)
      dat = None
    else:
      dat = self._load_pkl_gz(f_RL_FN)
    if (dat is not None):
      (self.f_cool, self.f_L_solvation, self.Psi, self.f_grid, \
       self.B, self.dock_equilibrated_cycle, self.dock_mean_acc) = dat
    else:
      for var in ['Psi', 'f_grid', 'B', 'dock_equilibrated_cycle', 'dock_mean_acc']:
        setattr(self,var,[])
    if readOnly:
      return

    # Make sure all the energies are available
    for c in range(self._dock_cycle):
      if len(self.dock_Es[-1][c].keys())==0:
        self.tee("  skipping the binding PMF calculation")
        return
    for var in ['f_L_solvation', 'f_cool', \
                'cool_equilibrated_cycle', 'cool_mean_acc']:
      if (not hasattr(self,var)) or (getattr(self,var)==[]):
        self.tee("  skipping the binding PMF calculation")
        return

    K = len(self.dock_protocol)
  
    # Calculate interaction energies
    for c in range(self._dock_cycle):
      while len(self.Psi) <= c:
        self.Psi.append({})
      if not 'MM' in self.Psi[c].keys():
        self.Psi[c]['MM'] = self.dock_Es[-1][c]['MM']/RT_TARGET
      if not 'grid' in self.Psi[c].keys():
        self.Psi[c]['grid'] = \
          (self.dock_Es[-1][c]['LJr'] + \
          self.dock_Es[-1][c]['LJa'] + \
          self.dock_Es[-1][c]['ELE'])/RT_TARGET
      for phase in self.params['dock']['phases']:
        if not phase in self.Psi[c].keys():
          self.Psi[c][phase] = \
            (self.dock_Es[-1][c]['RL'+phase][:,-1] - \
            self.dock_Es[-1][c]['L'+phase][:,-1] - \
            self.params['dock']['receptor_'+phase][-1])/RT_TARGET

    # Estimate cycle at which simulation has equilibrated
    u_Ks = [self._u_kln([self.dock_Es[-1][c]],[self.dock_protocol[-1]]) for c in range(self._dock_cycle)]
    mean_u_Ks = np.array([np.mean(u_K) for u_K in u_Ks])
    std_u_Ks = np.array([np.std(u_K) for u_K in u_Ks])
    for c in range(len(self.dock_equilibrated_cycle), self._dock_cycle):
      nearMean = (mean_u_Ks - mean_u_Ks[c])<std_u_Ks[c]
      for (item, ind) in zip(nearMean, range(len(nearMean))):
        if item:
          self.dock_equilibrated_cycle.append(ind)
          break

    # Calculate docking free energies that have not already been calculated
    updated = False
    for c in range(len(self.f_grid), self._dock_cycle):
      fromCycle = self.dock_equilibrated_cycle[c]
      toCycle = c + 1
      
      # Extract energies between fromCycle and toCycle
      dock_Es = [Es[fromCycle:toCycle] for Es in self.dock_Es]
      
      # Use MBAR for the grid scaling free energy estimate
      (u_kln,N_k) = self._u_kln(dock_Es,self.dock_protocol)
      f_grid = self._run_MBAR(u_kln,N_k)
      self.f_grid.append(f_grid)

      # Binding PMF
      B = {'grid':(-self.f_cool[-1][-1] + f_grid[-1])}
      u_K = self._u_kln([dock_Es[-1]],[self.dock_protocol[-1]])

# How to estimate the change in potential energy for GBSA:
# RLGBSAflex = RLGBSAfixedR - RGBSAfixedR + RGBSAflex
#                                           ^ This is assumed to cancel out
# RLMMTK = RMMTK + LMMTK + PsiMMTK
#          ^ This is assumed to be zero
# du = (RLGBSAfixedR - RGBSAfixedR + RGBSAflex) - (LMMTK + PsiMMTK)
#
# in code this is,
# du = \
#    (np.concatenate([self.dock_Es[-1][c]['RL'+phase][:,-1] \
#                      for c in range(fromCycle,toCycle)]) - \
#        self.params['dock']['receptor_'+phase])/RT_TARGET - \
#    u_K
#
# Alternatively,
# RLGBSAflex = RGBSAflex + LGBSA + PsiGBSA
#              ^ this has been missing
#            = RGBSAflex + LGBSA + (RLGBSAflex - RGBSAflex - LGBSA)
# du = (RGBSAflex + LGBSA + PsiGBSA) - (RMMTK + LMMTK + PsiMMTK)
#
# in code this is,
# du = self.dock_Es[-1][c]['L'+phase][:,-1]/RT_TARGET+self.Psi[c][phase] - u_K
# which gives the same numbers for c = 0
#
# Alternatively, assuming LGBSA~LMMTK
# du = PsiGBSA - PsiMMTK
# du = np.concatenate([self.Psi[c][phase] - self.Psi[c]['grid'] for c in range(fromCycle,toCycle)])

      for phase in self.params['dock']['phases']:
        du = np.concatenate([self.dock_Es[-1][c]['RL'+phase][:,-1] \
                            for c in range(fromCycle,toCycle)] - \
             self.params['dock']['receptor_'+phase][-1])/RT_TARGET - u_K
        min_du = min(du)
        # Complex solvation
        B_RL_solv = -np.log(np.exp(-du+min_du).mean()) + min_du
        weight = np.exp(-du+min_du)
        weight = weight/sum(weight)
        Psi = np.concatenate([self.Psi[c][phase] \
          for c in range(fromCycle,toCycle)])
        min_Psi = min(Psi)
        
        B[phase+'_RL_solv'] = B_RL_solv
        B[phase+'_MBAR_based'] = -self.f_L_solvation[-1][phase] - \
           self.f_cool[-1][-1] + f_grid[-1] + B_RL_solv
        B[phase+'_mean_Psi'] = np.mean(Psi)
        B[phase+'_min_Psi'] = min_Psi
        B[phase+'_inverse_FEP'] = \
          np.log(sum(weight*np.exp(Psi-min_Psi))) + min_Psi
      
      self.B.append(B)

      # Average acceptance probabilities
      dock_mean_acc = np.zeros(K-1)
      for k in range(0, K-1):
        (u_kln,N_k) = self._u_kln(dock_Es[k:k+2],self.dock_protocol[k:k+2])
        N = min(N_k)
        acc = np.exp(-u_kln[0,1,:N]-u_kln[1,0,:N]+u_kln[0,0,:N]+u_kln[1,1,:N])
        dock_mean_acc[k] = np.mean(np.minimum(acc,np.ones(acc.shape)))
      self.dock_mean_acc.append(dock_mean_acc)

      self.tee("  calculated %s binding PMF of %f RT using cycles %d to %d"%(\
        self.params['dock']['phases'][-1], B[self.params['dock']['phases'][-1]+'_MBAR_based'], fromCycle, toCycle-1))
      updated = True

    if updated:
      self._write_pkl_gz(f_RL_FN, (self.f_cool, self.f_L_solvation, self.Psi, \
        self.f_grid, self.B, self.dock_equilibrated_cycle, self.dock_mean_acc))

  ##############################
  # Analysis and Visualization #
  ##############################
  def calc_random_dock_stats(self, E):
    """
    Examines the convergence of the interaction energy 
    in the first docking state, looking at how the minimum, 
    mean, standard deviation, and other statistics vary as
    a function of the number of rotations and translations.
    """
    
    # Initiate variables
    interval = 5
    keys = ['min','mean','std','lambda_n','f_grid0','f_grid']

    r_range = range(interval, self._n_rot, interval)
    t_range = range(interval, self._n_trans, interval)
    
    random_dock_stats = {}
    for key in keys:
      random_dock_stats[key] = np.zeros((len(r_range),len(t_range)))

    Psi = np.reshape(E['ELE'],(-1,self._n_rot,self._n_trans)) + \
          np.reshape(E['LJr'],(-1,self._n_rot,self._n_trans)) + \
          np.reshape(E['LJa'],(-1,self._n_rot,self._n_trans))
    w = np.exp((-1/RT_TARGET+1/RT_HIGH)*np.reshape(E['MM'] - min(E['MM']),(-1,self._n_rot,self._n_trans)))
    
    for r in range(len(r_range)):
      for t in range(len(t_range)):
        Psi_c = np.ravel(Psi[:,:r_range[r],:t_range[t]])
        Psi_c_std = Psi_c.std()
        lambda_c = (self.params['dock']['therm_speed']*RT_HIGH)/Psi_c_std
        random_dock_stats['min'][r][t] = Psi_c.min()
        random_dock_stats['mean'][r][t] = Psi_c.mean()
        random_dock_stats['std'][r][t] = Psi_c_std
        random_dock_stats['lambda_n'][r][t] = lambda_c

        du = lambda_c*Psi_c/RT_HIGH
        min_du = min(du)
        random_dock_stats['f_grid0'][r][t] = -RT_HIGH*(np.log(np.mean(np.exp(-du+min_du))) + min_du)

        w_c = np.ravel(w[:,:r_range[r],:t_range[t]])

        du = Psi_c/RT_TARGET
        min_du = min(du)
        random_dock_stats['f_grid'][r][t] = -RT_TARGET*(np.log(np.sum(w_c*np.exp(-du+min_du))/np.sum(w_c)) + min_du)
    
    self._write_pkl_gz(join(self.dir['dock'],'random_dock_stats.pkl.gz'), random_dock_stats)
  
  ######################
  # Internal Functions #
  ######################

  def _set_universe_force_field(self,lambda_n):
    """
    Sets the force field parameters of the universe to values appropriate for the given lambda_n dictionary.
    The elements in the dictionary lambda_n can be:
      MM - True, to turn on the Generalized AMBER force field
      site - True, to turn on the binding site
      sLJr - scaling of the soft Lennard-Jones repulsive grid
      sLJa - scaling of the soft Lennard-Jones attractive grid
      sELE - scaling of the soft electrostatic grid
      LJr - scaling of the Lennard-Jones repulsive grid
      LJa - scaling of the Lennard-Jones attractive grid
      ELE - scaling of the electrostatic grid
      T - the temperature in K
    """

    self.T = lambda_n['T']
    self.RT = R*lambda_n['T']
    
    if 'delta_t' in lambda_n.keys():
      self.delta_t = lambda_n['delta_t']
    else:
      self.delta_t = 1.5*MMTK.Units.fs

    fflist = []
    if ('MM' in lambda_n.keys()) and lambda_n['MM']:
      fflist.append(self._forceFields['gaff'])
    if ('site' in lambda_n.keys()) and lambda_n['site']:
      fflist.append(self._forceFields['site'])
    for scalable in self._scalables:
      if (scalable in lambda_n.keys()) and lambda_n[scalable]>0:
        # Load the force field if it has not been loaded
        if not scalable in self._forceFields.keys():
          import time
          loading_start_time = time.time()
          grid_FN = self._FNs['grids'][{'sLJr':'LJr','sLJa':'LJa','sELE':'ELE',
            'LJr':'LJr','LJa':'LJa','ELE':'ELE'}[scalable]]
          grid_scaling_factor = 'scaling_factor_' + \
            {'sLJr':'LJr','sLJa':'LJa','sELE':'electrostatic', \
             'LJr':'LJr','LJa':'LJa','ELE':'electrostatic'}[scalable]
          if scalable=='LJr':
            from ForceFields.Grid.TrilinearISqrtGrid import TrilinearISqrtGridForceField
            self._forceFields[scalable] = TrilinearISqrtGridForceField(grid_FN,
              lambda_n[scalable], grid_scaling_factor,
              grid_name=scalable, max_val=-1)
          else:
            if scalable=='sLJr':
              max_val = 10.0
            elif scalable=='sELE':
              # The maximum value is set so that the electrostatic energy
              # less than or equal to the Lennard-Jones repulsive energy
              # for every heavy atom at every grid point
              scaling_factors_ELE = np.array([ \
                self.molecule.getAtomProperty(a, 'scaling_factor_electrostatic') \
                  for a in self.molecule.atomList()],dtype=float)
              scaling_factors_LJr = np.array([ \
                self.molecule.getAtomProperty(a, 'scaling_factor_LJr') \
                  for a in self.molecule.atomList()],dtype=float)
              scaling_factors_ELE = scaling_factors_ELE[scaling_factors_LJr>10]
              scaling_factors_LJr = scaling_factors_LJr[scaling_factors_LJr>10]
              max_val = min(abs(scaling_factors_LJr*10.0/scaling_factors_ELE))
            else:
              max_val = -1
              
            from ForceFields.Grid.TrilinearGrid import TrilinearGridForceField
            self._forceFields[scalable] = TrilinearGridForceField(grid_FN,
              lambda_n[scalable], grid_scaling_factor,
              grid_name=scalable, max_val=max_val)
          self.tee('  %s grid loaded from %s in %s'%(scalable, grid_FN, HMStime(time.time()-loading_start_time)))

        # Set the force field strength to the desired value
        self._forceFields[scalable].strength = lambda_n[scalable]
        fflist.append(self._forceFields[scalable])

    compoundFF = fflist[0]
    for ff in fflist[1:]:
      compoundFF += ff
    self.universe.setForceField(compoundFF)

  def _MC_translate_rotate(self, trials=25):
    """
    Conducts Monte Carlo translation and rotation moves.
    """

    step_size = min(\
      1.0*MMTK.Units.Ang, \
      self._forceFields['site'].max_R*MMTK.Units.nm/10)

    acc = 0
    xo = self.universe.copyConfiguration()
    eo = self.universe.energy()
    com = self.universe.centerOfMass().array
    
    for c in range(trials):
      step = np.random.randn(3)*step_size
      xn = np.dot((xo.array - com), random_rotate()) + com + step
      self.universe.setConfiguration(Configuration(self.universe,xn))
      en = self.universe.energy()
      if ((en<eo) or (np.random.random()<np.exp(-(en-eo)/self.RT))):
        acc += 1
        xo = self.universe.copyConfiguration()
        eo = en
        com += step
      else:
        self.universe.setConfiguration(xo)

    return acc

  def _initial_sim_state(self, process, seeds, normalize=False, doMC=False):
    """
    Initializes a state, returning the configurations and potential energy.
    """
    
    all_confs = []
    all_potEs = []
    Hs = []

    for seed in seeds:
      self.universe.setConfiguration(Configuration(self.universe,seed))
      if doMC and process=='dock' and self.params['dock']['MCMC_moves']>0:
        self._MC_translate_rotate()
      (confs_s, potEs_s, Hs_s, self.delta_t) = self.sampler['init'](
        steps = self.params[process]['steps_per_seed'],
        T=self.T, delta_t=self.delta_t,
        steps_per_trial = self.params[process]['steps_per_seed']/10,
        normalize=normalize, adapt=True)
      # Restrict the range of the time step
      self.delta_t = min(max(self.delta_t, 0.25*MMTK.Units.fs), 2.5*MMTK.Units.fs)
      all_confs.append(confs_s)
      all_potEs.append(np.array(potEs_s))
      Hs.append(Hs_s)
    Hs = np.mean(np.array(Hs))

    import sys
    original_stderr = sys.stderr
    sys.stderr = NullDevice()

    # Decorrelation
    if len(seeds)>1:
      tau_ac = pymbar.timeseries.integratedAutocorrelationTimeMultiple(all_potEs)
    else:
      tau_ac = pymbar.timeseries.integratedAutocorrelationTime(all_potEs)
    # stride = max(int(np.ceil(1+2*tau_ac)),1)
    
    sys.stderr = original_stderr

    confs = []
    potEs = []
    for s in range(len(seeds)):
      nsnaps = len(all_confs[s])
      # uncorr_indicies = range(min(stride-1, nsnaps-1), nsnaps, stride)
      # Only use the last configuration from each seed
      uncorr_indicies = [nsnaps-1]
      confs.extend([all_confs[s][i] for i in uncorr_indicies])
      potEs.extend([all_potEs[s][i] for i in uncorr_indicies])

    return (confs, np.array(potEs), tau_ac)
    
  def _replica_exchange(self, process):
    """
    Performs a cycle of replica exchange
    """
    if not process in ['dock','cool']:
      raise Exception('Process must be dock or cool')

    self._set_lock(process)

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

    # Setting the force field will load grids
    # before multiple processes are spawned
    for k in range(K):
      self._set_universe_force_field(lambdas[k])
    
    storage = {}
    for var in ['confs','state_inds','energies']:
      storage[var] = []
    
    cycle_start_time = time.time()
    
    # Do replica exchange
    MC_time = 0
    state_inds = range(K)
    inv_state_inds = range(K)
    for sweep in range(self.params[process]['sweeps_per_cycle']):
      E = {}
      for term in terms:
        E[term] = np.zeros(K, dtype=float)
      if process=='dock':
        MC_start_time = time.time()
        E['acc_MC'] = np.zeros(K, dtype=float)
      Ht = np.zeros(K, dtype=int)
      # Sample within each state
      results = []
      if self._do_multiprocess:
        result_queue = multiprocessing.Queue()
        processes = []
        for k in range(K):
            p = multiprocessing.Process(target=self._sample_state, \
              args=(result_queue, process, lambdas[state_inds[k]], \
                    Configuration(self.universe, confs[k])))
            processes.append(p)
        for pro in processes:
            pro.start()
        for pro in processes:
            pro.join()
        for k in range(K):
            results.append(result_queue.get())
        for pro in processes:
            pro.terminate()
      else:
        for k in range(K):
          results.append(self._sample_state(None, process, \
            lambdas[state_inds[k]], Configuration(self.universe, confs[k])))
      # Store results
      for k in range(K):
        if 'acc_MC' in results[k].keys():
          E['acc_MC'][k] = results[k]['acc_MC']
          MC_time += results[k]['MC_time']
        confs[k] = results[k]['confs']
        if process == 'cool':
            E['MM'][k] = results[k]['E_MM']
        Ht[k] += results[k]['Ht']
      if process=='dock':
        E = self._calc_E(confs, E) # Get energies
        # Get rmsd values
        if self.params['dock']['rmsd'] is not False:
          E['rmsd'] = np.array([np.sqrt(((confs[k][self.molecule.heavy_atoms,:] - \
            self.confs['rmsd'])**2).sum()/self.molecule.nhatoms) for k in range(K)])
      # Calculate u_ij (i is the replica, and j is the configuration),
      #    a list of arrays
      (u_ij,N_k) = self._u_kln(E, [lambdas[state_inds[c]] for c in range(K)])
      # Do the replica exchange
      for attempt in range(1000):
        for (t1,t2) in pairs_to_swap:
          a = inv_state_inds[t1]
          b = inv_state_inds[t2]
          ddu = -u_ij[a][b]-u_ij[b][a]+u_ij[a][a]+u_ij[b][b]
          if (ddu>0) or (np.random.uniform()<np.exp(ddu)):
            u_ij[a],u_ij[b] = u_ij[b],u_ij[a]
            state_inds[a],state_inds[b] = state_inds[b],state_inds[a]
            inv_state_inds[state_inds[a]],inv_state_inds[state_inds[b]] = \
              inv_state_inds[state_inds[b]],inv_state_inds[state_inds[a]]
      # Store data in local variables
      storage['confs'].append(list(confs))
      storage['state_inds'].append(list(state_inds))
      storage['energies'].append(copy.deepcopy(E))

    Ht = Ht/self.params[process]['sweeps_per_cycle']
    storage['Ht'] = Ht

    # Estimate relaxation time from empirical state transition matrix
    state_inds = np.array(storage['state_inds'])
    Nij = np.zeros((K,K),dtype=int)
    for (i,j) in zip(state_inds[:-1,:],state_inds[1:,:]):
      for k in range(K):
        Nij[j[k],i[k]] += 1
    N = (Nij+Nij.T)
    Tij = np.array(N,dtype=float)/sum(N,1)
    (eval,evec)=np.linalg.eig(Tij)
    tau2 = 1/(1-eval[1])
    
    # Estimate relaxation time from autocorrelation
    tau_ac = pymbar.timeseries.integratedAutocorrelationTimeMultiple(state_inds.T)
    # There will be at least 25 and up to sweeps_per_cycle saved samples
    # max(int(np.ceil((1+2*tau_ac)/50)),1) is the minimum stride,
    # which is based due to 25 samples per autocorrelation time.
    # max(self.params['dock']['sweeps_per_cycle']/25)
    # is the maximum stride to have 25 samples, if possible.
    stride = min(max(int(np.ceil((1+2*tau_ac)/25.)),1), \
                 max(int(np.ceil(self.params['dock']['sweeps_per_cycle']/25.)),1))

    store_indicies = np.array(\
      range(min(stride-1,self.params[process]['sweeps_per_cycle']-1), self.params[process]['sweeps_per_cycle'], stride), dtype=int)
    nsaved = len(store_indicies)

    self.tee("  storing %d configurations for %d replicas"%(nsaved, len(confs)) + \
      " in cycle %d"%cycle + \
      " (tau2=%f, tau_ac=%f)"%(tau2,tau_ac) + \
      " with %s for MC"%(HMStime(MC_time)) + \
      " in " + HMStime(time.time()-cycle_start_time))

    # Get indicies for storing global variables
    inv_state_inds = np.zeros((nsaved,K),dtype=int)
    for snap in range(nsaved):
      state_inds = storage['state_inds'][store_indicies[snap]]
      for state in range(K):
        inv_state_inds[snap][state_inds[state]] = state

    # Reorder energies and replicas for storage
    if process=='dock':
      terms.append('acc_MC') # Make sure to save the acceptance probability
      if self.params['dock']['rmsd'] is not False:
        terms.append('rmsd') # Make sure to save the rmsd
    Es = []
    for state in range(K):
      E_state = {}
      if state==0:
        E_state['repXpath'] = storage['state_inds']
      for term in terms:
        E_state[term] = np.array([storage['energies'][store_indicies[snap]][term][inv_state_inds[snap][state]] for snap in range(nsaved)])
      Es.append([E_state])

    self.confs[process]['replicas'] = \
      [storage['confs'][store_indicies[-1]][inv_state_inds[-1][state]] \
       for state in range(K)]

    # If it is not the last cycle, consider
    # refinings the protocol by inserting states between neighbors with low
    # mean replica exchange probabilities

    if ((not self.params[process]['no_protocol_refinement']) and \
        ((getattr(self,'_%s_total_cycle'%process)+1) \
          < self.params[process]['repX_cycles'])):
      # Estimate mean replica exchange probabilities between neighbors
      mean_acc = np.zeros(K-1)
      for k in range(K-1):
        (u_kln,N_k) = self._u_kln(Es[k:k+2],lambdas[k:k+2])
        N = min(N_k)
        acc = np.exp(-u_kln[0,1,:N]-u_kln[1,0,:N]+u_kln[0,0,:N]+u_kln[1,1,:N])
        mean_acc[k] = np.mean(np.minimum(acc,np.ones(acc.shape)))
      insert_state = (mean_acc<self.params['dock']['min_repX_acc'])
      if insert_state.any():
        for k in reversed(range(K-1)):
          if insert_state[k]:
            self.tee("  due to exchange acceptance probability of %f, inserted state between %d and %d"%(mean_acc[k],k,k+1))
            # Duplicating the configuration of the later state
            new_conf = np.copy(self.confs[process]['replicas'][k+1])
            self.confs[process]['replicas'].insert(k+1, new_conf)
        K = len(self.confs[process]['replicas'])
        self.tee("  refined protocol to have %d states after total cycle %d"%(\
          K, getattr(self,'_%s_total_cycle'%process)))
        
        # Determine new protocol with equal intervals in thermodynamic length
        progress_o = np.array([l['a'] for l in lambdas])
        tL_tensor = np.array([self._tL_tensor(Es[k],lambdas[k],process) \
          for k in range(len(lambdas))])
        tL_intervals = abs(progress_o[1:]-progress_o[:-1])*(tL_tensor[:-1] + tL_tensor[1:])
        tL_o = np.insert(np.cumsum(tL_intervals),0,0.)
        tL_n = np.linspace(0,tL_o[-1],K)
        progress_n = np.interp(tL_n,tL_o,progress_o)
        progress_n.sort()
        
        lambdas = [self._lambda(a,process) for a in progress_n]
        setattr(self,process+'_protocol', lambdas)

        # Clear progress and save the data
        setattr(self, process+'_Es', [[] for state in range(K)])
        self.confs[process]['samples'] = [[] for state in range(K)]
        setattr(self,'_%s_cycle'%process,0)
        setattr(self,'_%s_total_cycle'%process,
          getattr(self,'_%s_total_cycle'%process)+1)
        self._save_progress(process)
        if process=='cool':
          self.calc_f_L(redo=True,readOnly=True)
        else:
          self.calc_f_RL(redo=True,readOnly=True)
        self._clear_lock(process)
        return

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

    setattr(self,'_%s_cycle'%process,cycle + 1)
    setattr(self,'_%s_total_cycle'%process,
      getattr(self,'_%s_total_cycle'%process) + 1)
    self._save_progress(process)
    self._clear_lock(process)

  def _sample_state(self, results_queue, process, lambda_k, initial_conf):
    self._set_universe_force_field(lambda_k)
    self.universe.setConfiguration(initial_conf)

    T = lambda_k['T']
    
    if 'delta_t' in lambda_k.keys():
      delta_t = lambda_k['delta_t']
    else:
      delta_t = 1.5*MMTK.Units.fs
    
    results = {}
    # Perform MCMC moves
    if (process == 'dock') and (self.params['dock']['MCMC_moves'] > 0) and \
        (lambda_k['a'] < 0.1):
      MC_start_time = time.time()
      results['acc_MC'] = self._MC_translate_rotate(25)/25.
      results['MC_time'] = (time.time() - MC_start_time)

    # Execute sampler
    (confs_k, potEs_k, Hs_k, delta_t) = self.sampler[process](\
      steps=self.params[process]['steps_per_sweep'], \
      T=lambda_k['T'], delta_t=delta_t, \
      normalize=(process=='cool'),
      adapt=False)

    # Store results in dictionary
    results['confs'] = np.copy(confs_k[-1])
    if process == 'cool':
      results['E_MM'] = potEs_k[-1]
    results['Ht'] = Hs_k

    if results_queue is not None:
      results_queue.put(results)  # put the result into the result_queue
    else:
      return results

  def _sim_process(self, process):
    """
    Simulate and analyze a cooling or docking process.
    
    As necessary, first conduct an initial cooling or docking
    and then run a desired number of replica exchange cycles.
    """
    if (getattr(self,process+'_protocol')==[]) or \
       (not getattr(self,process+'_protocol')[-1]['crossed']):
      getattr(self,'initial_'+process)()
      self.cycles_run += 1
      if self.run_type=='one_step':
        return

    # Main loop for replica exchange
    if (self.params[process]['repX_cycles'] is not None) and \
       ((getattr(self,'_%s_total_cycle'%process) < self.params[process]['repX_cycles']) or (getattr(self,'_%s_cycle'%process)==0)):
      # Score configurations from another program
      if (process=='dock') and \
         (self._dock_cycle==1) and (self._dock_total_cycle==1) and \
         (self.params['dock']['score'] is not False):
        self._get_dock_score_confs()

      self.tee("\n>>> Replica exchange sampling for the {0}ing process".format(process), process=process)
      import time
      repEx_start_time = time.time()
      start_cycle = getattr(self,'_%s_total_cycle'%process)
      while ((getattr(self,'_%s_total_cycle'%process) < self.params[process]['repX_cycles']) or (getattr(self,'_%s_cycle'%process)==0)):
        self._replica_exchange(process)
        self.cycles_run += 1
        if self.run_type=='one_step' and process=='dock':
            break
      self.tee("Elapsed time for %d total cycles of replica exchange was %s\n"%(\
         (getattr(self,'_%s_total_cycle'%process)-start_cycle), \
         HMStime((time.time()-repEx_start_time))), process=process)

    # Do additional replica exchange on the cooling process
    #   if there are not enough configurations
    if (process=='cool'):
      E_MM = []
      for k in range(len(self.cool_Es[0])):
        E_MM += list(self.cool_Es[0][k]['MM'])
      while len(E_MM)<self.params['dock']['seeds_per_state']:
        self.tee("More samples from high temperature ligand simulation needed", process='cool')
        self._replica_exchange('cool')
        E_MM = []
        for k in range(len(self.cool_Es[0])):
          E_MM += list(self.cool_Es[0][k]['MM'])

    self._set_lock(process)
    if process=='cool':
      self._postprocess(conditions=[('cool',-1,-1,'L')])
      self.calc_f_L()
    elif process=='dock':
      if not (((getattr(self,'_%s_total_cycle'%process) < self.params[process]['repX_cycles']) or (getattr(self,'_%s_cycle'%process)==0))):
        self._postprocess()
        self.calc_f_RL()
    self._clear_lock(process)

  def _get_dock_score_confs(self):
    """
    Gets configurations to score from another program.
    Returns unique configurations and corresponding energies (if applicable).
    The final results are stored in self.confs['dock']['replicas'].
    """
    self._set_lock('dock')
    minimized_confs = []
    minimized_conf_Es = []
    unique_confs = []
    if self.params['dock']['score'] is True:
      unique_confs = [self.confs['ligand']]
      self.confs['dock']['replicas'] = [np.copy(self.confs['ligand']) \
        for k in range(len(self.dock_protocol))]
      self.tee("\n>>> Rescoring default configuration")
    elif self.params['dock']['score'] is False:
      raise Exception('Error running _get_dock_score_confs with score==False')
    else:
      count = {'dock6':0, 'initial_dock':0, 'duplicated':0}
      
      if self.params['dock']['score'].endswith('.mol2') or \
         self.params['dock']['score'].endswith('.mol2.gz'):
        (confs_dock6,E_mol2) = self._read_dock6(self.params['dock']['score'])
        # Add configurations where the ligand is in the binding site
        self._set_universe_force_field({'site':True,'T':T_TARGET})
        for ind in range(len(confs_dock6)):
          self.universe.setConfiguration(Configuration(self.universe, confs_dock6[ind]))
          if self.universe.energy()<1.:
            unique_confs.append(confs_dock6[ind])
        count['dock6'] = len(unique_confs)
      else:
        raise Exception('Unrecognized file type for configurations')

      if self.confs['dock']['seeds'] is not None:
        unique_confs = unique_confs + self.confs['dock']['seeds']
        count['initial_dock'] = len(self.confs['dock']['seeds'])
      
      # Minimize each configuration
      self._set_universe_force_field(self._lambda(1.0,process='dock'))
      from MMTK.Minimization import SteepestDescentMinimizer # @UnresolvedImport
      minimizer = SteepestDescentMinimizer(self.universe)

      self.tee("\n>>> Minimizing seed configurations")
      min_start_time = time.time()
      for conf in unique_confs:
        self.universe.setConfiguration(Configuration(self.universe, conf))
        minimizer(steps = 1000)
        minimized_confs.append(np.copy(self.universe.configuration().array))
        minimized_conf_Es.append(self.universe.energy())
      self.tee("\nElapsed time for minimization: " + \
        HMStime(time.time()-min_start_time))

      # Sort minimized configurations
      minimized_conf_Es, minimized_confs = \
        (list(l) for l in zip(*sorted(zip(minimized_conf_Es, minimized_confs), \
          key=lambda p:p[0], reverse=True)))
  
      self.confs['dock']['replicas'] = minimized_confs[-len(self.dock_protocol):]
      while len(self.confs['dock']['replicas'])<len(self.dock_protocol):
        self.confs['dock']['replicas'].append(self.confs['dock']['replicas'][-1])
        count['duplicated'] += 1
      self.tee("\n>>> Rescoring configurations:")
      self.tee("  {dock6} from dock6, {initial_dock} from initial docking, and {duplicated} duplicated".format(**count))
    self._clear_lock('dock')
    return (minimized_confs, minimized_conf_Es)

  def _run_MBAR(self,u_kln,N_k):
    """
    Estimates the free energy of a transition using MBAR
    """
    import pymbar
    K = len(N_k)
    f_k_FEPF = np.zeros(K)
    f_k_FEPR = np.zeros(K)
    f_k_BAR = np.zeros(K)
    for k in range(K-1):
      w_F = u_kln[k,k+1,:N_k[k]] - u_kln[k,k,:N_k[k]]
      min_w_F = min(w_F)
      w_R = u_kln[k+1,k,:N_k[k+1]] - u_kln[k+1,k+1,:N_k[k+1]]
      min_w_R = min(w_R)
      f_k_FEPF[k+1] = -np.log(np.mean(np.exp(-w_F+min_w_F))) + min_w_F
      f_k_FEPR[k+1] = np.log(np.mean(np.exp(-w_R+min_w_R))) - min_w_R
      try:
        f_k_BAR[k+1] = pymbar.BAR(w_F, w_R, relative_tolerance=0.000001, verbose=False, compute_uncertainty=False)
      except:
        f_k_BAR[k+1] = f_k_FEPF[k+1]
    f_k_FEPF = np.cumsum(f_k_FEPF)
    f_k_FEPR = np.cumsum(f_k_FEPR)
    f_k_BAR = np.cumsum(f_k_BAR)
    try:
      f_k = pymbar.MBAR(u_kln, N_k,
        verbose = False, method = 'adaptive',
        initial_f_k = f_k_BAR,
        maximum_iterations = 20).f_k
    except:
      f_k = f_k_BAR
    return f_k

  def _u_kln(self,eTs,lambdas,noBeta=False):
    """
    Computes a reduced potential energy matrix.  k is the sampled state.  l is the state for which energies are evaluated.
    
    Input:
    eT is a 
      -dictionary (of mapped energy terms) of numpy arrays (over states)
      -list (over states) of dictionaries (of mapped energy terms) of numpy arrays (over configurations), or a
      -list (over states) of lists (over cycles) of dictionaries (of mapped energy terms) of numpy arrays (over configurations)
    lambdas is a list of thermodynamic states
    noBeta means that the energy will not be divided by RT
    
    Output: u_kln or (u_kln, N_k)
    u_kln is the matrix (as a numpy array)
    N_k is an array of sample sizes
    """
    L = len(lambdas)

    addMM = ('MM' in lambdas[0].keys()) and (lambdas[0]['MM'])
    addSite = ('site' in lambdas[0].keys()) and (lambdas[0]['site'])
    probe_key = [key for key in lambdas[0].keys() if key in (['MM'] + self._scalables)][0]
    
    if isinstance(eTs,dict):
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
        E = 1.*E_base
        for scalable in self._scalables:
          if scalable in lambdas[l].keys():
            E += lambdas[l][scalable]*eTs[scalable]
        if noBeta:
          u_kln.append(E)
        else:
          u_kln.append(E/(R*lambdas[l]['T']))
    elif isinstance(eTs[0],dict):
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
        E = 1.*E_base
        for scalable in self._scalables:
          if scalable in lambdas[l].keys():
            E += lambdas[l][scalable]*eTs[k][scalable]
        if noBeta:
          u_kln[k,l,:N_k[k]] = E
        else:
          u_kln[k,l,:N_k[k]] = E/(R*lambdas[l]['T'])
    elif isinstance(eTs[0],list):
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
          E = 1.*E_base
          for scalable in self._scalables:
            if scalable in lambdas[l].keys():
              E += lambdas[l][scalable]*np.concatenate([eTs[k][c][scalable] for c in range(C)])
          if noBeta:
            u_kln[k,l,:N_k[k]] = E
          else:
            u_kln[k,l,:N_k[k]] = E/(R*lambdas[l]['T'])

    if (K==1) and (L==1):
      return u_kln.ravel()
    else:
      return (u_kln,N_k)

  def _next_dock_state(self, E=None, lambda_o=None, pow=None):
    """
    Determines the parameters for the next docking state
    """

    if E is None:
      E = self.dock_Es[-1]

    if lambda_o is None:
      lambda_o = self.dock_protocol[-1]
    lambda_n = copy.deepcopy(lambda_o)
    
    if self.params['dock']['protocol']=='Set':
      raise Exception("Set protocol not currently supported")
    elif self.params['dock']['protocol']=='Adaptive':
      # Change grid scaling and temperature simultaneously
      tL_tensor = self._tL_tensor(E,lambda_o)
      if pow is not None:
        tL_tensor = tL_tensor*(1.25**pow)
      if tL_tensor>0:
        dL = self.params['dock']['therm_speed']/tL_tensor
        a = min(lambda_o['a'] + dL, 1)
        return self._lambda(a)
        if a == 1:
          lambda_n['crossed'] = True
      else:
        raise Exception("No variance in stage")

  def _tL_tensor(self, E, lambda_c, process='dock'):
    T = lambda_c['T']
    if process=='dock':
      # Metric tensor for the thermodynamic length
      a = lambda_c['a']
      a_sg = 1.-4.*(a-0.5)**2
      a_g = 4.*(a-0.5)**2/(1+np.exp(-100*(a-0.5)))
      da_sg_da = -8*(a-0.5)
      da_g_da = (400.*(a-0.5)**2*np.exp(-100.*(a-0.5)))/(1+np.exp(-100.*(a-0.5)))**2 + \
        (8.*(a-0.5))/(1 + np.exp(-100.*(a-0.5)))
      Psi_sg = self._u_kln([E], [{'sLJr':1,'sELE':1}], noBeta=True)
      Psi_g = self._u_kln([E], [{'LJr':1,'LJa':1,'ELE':1}], noBeta=True)
      U_RL_g = self._u_kln([E],
        [{'MM':True, 'site':True, 'T':T,\
        'sLJr':a_sg, 'sELE':a_sg, 'LJr':a_g, 'LJa':a_g, 'ELE':a_g}], noBeta=True)
      return np.abs(da_sg_da)*Psi_sg.std()/(R*T) + \
             np.abs(da_g_da)*Psi_g.std()/(R*T) + \
             np.abs(T_TARGET-T_HIGH)*U_RL_g.std()/(R*T*T)
    elif process=='cool':
      return self._u_kln([E],[{'MM':True}], noBeta=True).std()/(R*T*T)
    else:
      raise Exception("Unknown process!")

  def _lambda(self, a, process='dock', lambda_o=None):
    if process=='dock':
      if lambda_o is None:
        lambda_o = self.dock_protocol[-1]
      lambda_n = copy.deepcopy(lambda_o)
      a_sg = 1.-4.*(a-0.5)**2
      a_g = 4.*(a-0.5)**2/(1+np.exp(-100*(a-0.5)))
      if a_g<1E-10:
        a_g=0
      lambda_n['a'] = a
      lambda_n['sLJr'] = a_sg
      lambda_n['sELE'] = a_sg
      lambda_n['LJr'] = a_g
      lambda_n['LJa'] = a_g
      lambda_n['ELE'] = a_g
      lambda_n['T'] = a*(T_TARGET-T_HIGH) + T_HIGH
      lambda_n['crossed'] = (abs(a-1.0)<0.001)
    elif process=='cool':
      if lambda_o is None:
        lambda_o = self.cool_protocol[-1]
      lambda_n = copy.deepcopy(lambda_o)
      lambda_n['a'] = a
      lambda_n['T'] = T_HIGH - a*(T_HIGH-T_TARGET)
      lambda_n['crossed'] = (abs(a-1.0)<0.001)
    else:
      raise Exception("Unknown process!")
    return lambda_n

  def _postprocess(self,
      conditions=[('original',0, 0,'R'), ('cool',-1,-1,'L'), \
                  ('dock',   -1,-1,'L'), ('dock',-1,-1,'RL')],
      phases=None,
      readOnly=False, redo_dock=False, debug=False):
    """
    Obtains the NAMD energies of all the conditions using all the phases.  
    Saves both MMTK and NAMD energies after NAMD energies are estimated.
    """
    updated = []
    crd_FN_o = ''
    if phases is None:
      phases = self.params['dock']['phases']

    postprocess_start_time = time.time()

    # state == -1 means the last state
    # cycle == -1 means all cycles
    for (p, state, cycle, moiety) in conditions:
      # Check that the values are legitimate
      if not p in ['cool','dock','original']:
        raise Exception("Type should be in ['cool', 'dock', 'original']")
      if not moiety in ['R','L', 'RL']:
        raise Exception("Species should in ['R','L', 'RL']")
    
      if p!='original' and getattr(self,p+'_protocol')==[]:
        continue
    
      if state==-1:
        state = len(getattr(self,p+'_protocol'))-1

      if cycle==-1:
        cycles = range(getattr(self,'_'+p+'_cycle'))
      else:
        cycles = [cycle]
    
      for c in cycles:
        if p=='original':
          prefix = p
        else:
          prefix = '%s%d_%d'%(p, state, c)

        for phase in phases:
          p_dir = {'cool':self.dir['cool'],
                 'original':self.dir['dock'],
                 'dock':self.dir['dock']}[p]
          
          label = moiety+phase
          AMBER_mdcrd_FN = join(p_dir,'%s.%s.mdcrd'%(prefix,moiety))
          dcd_FN = join(p_dir,'%s.%s.dcd'%(prefix,moiety))
          outputname = join(p_dir,'%s.%s%s'%(prefix,moiety,phase))
          
          # Skip NAMD
          # if the function is NOT being rerun in redo mode
          # and one of the following:
          # the function is being run in readOnly mode,
          # the energies are already in memory.
          if (not (redo_dock and p=='dock')) and \
            (readOnly \
            or (p == 'original' and \
                (label in getattr(self,p+'_Es')[state][c].keys()) and \
                (getattr(self,p+'_Es')[state][c][label] is not None)) \
            or (('MM' in getattr(self,p+'_Es')[state][c].keys()) and \
                (label in getattr(self,p+'_Es')[state][c].keys()) and \
                (len(getattr(self,p+'_Es')[state][c]['MM'])==\
                 len(getattr(self,p+'_Es')[state][c][label])))):
            pass
          else:
            # Do the calculation
            # Obtain configurations
            if (moiety=='R'):
              if not 'receptor' in self.confs.keys():
                continue
              confs = self.confs['receptor']
            else:
              confs = self.confs[p]['samples'][state][c]
          
            if phase in ['NAMD_Gas','NAMD_GBSA']:
              E = self._NAMD_Energy(confs, moiety, phase, \
                dcd_FN, outputname, debug=debug)
            else:
              E = self._sander_Energy(confs, moiety, phase, \
                AMBER_mdcrd_FN, outputname, debug=debug)

            getattr(self,p+'_Es')[state][c][label] = E
          
            # Updated
            if not p in updated:
              updated.append(p)

        # Clean up
        if (not debug) and isfile(dcd_FN):
          os.remove(dcd_FN)
        if (not debug) and isfile(AMBER_mdcrd_FN):
          os.remove(AMBER_mdcrd_FN)

    # Save data
    if 'original' in updated:
      for phase in phases:
        if (self.params['dock']['receptor_'+phase] is None) and \
           (self.original_Es[0][0]['R'+phase] is not None):
          self.params['dock']['receptor_'+phase] = \
            self.original_Es[0][0]['R'+phase][-1]

    if 'cool' in updated:
      self._save_progress('cool')
    if ('dock' in updated) or ('original' in updated):
      self._save_progress('dock')

    if len(updated)>0:
      self.tee("  postprocessed data in " + \
        HMStime(time.time()-postprocess_start_time))

  def _calc_E(self, confs, E=None, type='sampling', prefix='confs'):
    """
    Calculates energies for a series of configurations
    Units are the MMTK standard, kJ/mol
    """
    if E is None:
      E = {}

    if type=='sampling' or type=='all':
      # Molecular mechanics and grid interaction energies
      self._set_universe_force_field(self.lambda_full)
      for term in (['MM','site','misc'] + self._scalables):
        E[term] = np.zeros(len(confs), dtype=float)
      for c in range(len(confs)):
        self.universe.setConfiguration(Configuration(self.universe,confs[c]))
        eT = self.universe.energyTerms()
        for (key,value) in eT.iteritems():
          E[term_map[key]][c] += value

    if type=='all':
      for phase in self.params['dock']['phases']:
        E['R'+phase] = self.params['dock']['receptor_'+phase]
        for moiety in ['L','RL']:
          dcd_FN = join(self.dir['dock'],'%s.%s.dcd'%(prefix,moiety))
          mdcrd_FN = join(self.dir['dock'],'%s.%s.mdcrd'%(prefix,moiety))
          outputname = join(self.dir['dock'],'%s.%s%s'%(prefix,moiety,phase))
          if phase in ['NAMD_Gas','NAMD_GBSA']:
            E[moiety+phase] = self._NAMD_Energy(confs, moiety, phase, dcd_FN, outputname)
          else:
            E[moiety+phase] = self._sander_Energy(confs, moiety, phase, mdcrd_FN, outputname)
    return E
  
  def _sander_Energy(self, confs, moiety, phase, AMBER_mdcrd_FN, outputname=None, debug=False):
    if not isfile(AMBER_mdcrd_FN):
      self._write_mdcrd(AMBER_mdcrd_FN, confs,
        includeReceptor=(moiety.find('R')>-1),
        includeLigand=(moiety.find('L')>-1))

    script_FN = '%s%s.in'%('.'.join(AMBER_mdcrd_FN.split('.')[:-1]),phase)
    out_FN = '%s%s.out'%('.'.join(AMBER_mdcrd_FN.split('.')[:-1]),phase)

    script_F = open(script_FN,'w')
    if phase=='GBSA':
      script_F.write('''Calculate GBSA energies
&cntrl
  imin=5,    ! read trajectory in for analysis
  ntx=1,     ! input is read formatted with no velocities
  irest=0,
  ntb=0,     ! no periodicity and no PME
  idecomp=0, ! no decomposition
  ntc=1,     ! No SHAKE
  ntf=1,     ! Complete interaction is calculated
  cut=9999., !
  igb=8,     ! Most recent AMBER GBn model, best agreement with PB
  gbsa=2,    ! recursive surface area algorithm
#  alpb=1,    ! Analytical Linearized Poisson-Boltzmann
#  arad=24.7  ! Electrostatic size
/
''')
    elif phase=='ALPB':
      # UNDER CONSTRUCTION
      # TO DO: Write something to get elsize
      #        run ambpdb to generate a pqr file
      #        run elsize
      #        save elsize as a parameter
      script_F.write('''Calculate GBSA energies
&cntrl
  imin=5,    ! read trajectory in for analysis
  ntx=1,     ! input is read formatted with no velocities
  irest=0,
  ntb=0,     ! no periodicity and no PME
  idecomp=0, ! no decomposition
  ntc=1,     ! No SHAKE
  ntf=1,     ! Complete interaction is calculated
  cut=9999., !
  igb=7,     ! Most recent AMBER GBn model, best agreement with PB
  gbsa=2,    ! recursive surface area algorithm
  alpb=1,    ! Analytical Linearized Poisson-Boltzmann
  arad=%f  ! Electrostatic size
/
'''%elsize)
    elif phase=='Gas':
      script_F.write('''Calculate Gas energies
&cntrl
  imin=5,    ! read trajectory in for analysis
  ntx=1,     ! input is read formatted with no velocities
  irest=0,
  ntb=0,     ! no periodicity and no PME
  idecomp=0, ! no decomposition
  ntc=1,     ! No SHAKE
  ntf=1,     ! Complete interaction is calculated
  cut=9999., !
/
''')
    script_F.close()
    
    os.chdir(os.path.dirname(out_FN))
    
    # Decompress prmtop and inpcrd files
    decompress = (self._FNs['prmtop'][moiety].endswith('.gz')) or \
                 (self._FNs['inpcrd'][moiety].endswith('.gz'))
    if decompress:
      for key in ['prmtop','inpcrd']:
        if self._FNs[key][moiety].endswith('.gz'):
          import shutil
          shutil.copy(self._FNs[key][moiety],self._FNs[key][moiety]+'.BAK')
          os.system('gunzip -f '+self._FNs[key][moiety])
          os.rename(self._FNs[key][moiety]+'.BAK', self._FNs[key][moiety])
          self._FNs[key][moiety] = self._FNs[key][moiety][:-3]

    import subprocess
    p = subprocess.Popen([self._FNs['sander'], '-O','-i',script_FN,'-o',out_FN, \
      '-p',self._FNs['prmtop'][moiety],'-c',self._FNs['inpcrd'][moiety], \
      '-y',AMBER_mdcrd_FN])
    p.wait()
    if not debug and isfile(script_FN):
      os.remove(script_FN)
    if isfile('restrt'):
      os.remove('restrt')
    os.chdir(self.dir['start'])

    # Clear decompressed files
    if decompress:
      for key in ['prmtop','inpcrd']:
        if isfile(self._FNs[key][moiety]+'.gz'):
          os.remove(self._FNs[key][moiety])
          self._FNs[key][moiety] = self._FNs[key][moiety] + '.gz'

    print out_FN

    F = open(out_FN,'r')
    dat = F.read().strip().split(' BOND')
    F.close()

    dat.pop(0)
    if len(dat)>0:
      E = np.array([rec[:rec.find('\nminimization')].replace('1-4 ','1-4').split()[1::3] for rec in dat],dtype=float)*MMTK.Units.kcal/MMTK.Units.mol
      E = np.hstack((E,np.sum(E,1)[...,None]))
    else:
      E = np.array([np.inf]*11)

    if not debug and isfile(out_FN):
      os.remove(out_FN)

    return E
    # AMBER ENERGY FIELDS:
    # For Gas phase:
    # 0. BOND 1. ANGLE 2. DIHEDRAL 3. VDWAALS 4. EEL 5. HBOND 6. 1-4 VWD 7. 1-4 EEL 8. RESTRAINT
    # For GBSA phase:
    # 0. BOND 1. ANGLE 2. DIHEDRAL 3. VDWAALS 4. EEL 5. EGB 6. 1-4 VWD 7. 1-4 EEL 8. RESTRAINT
    # 9. ESURF

  def _NAMD_Energy(self, confs, moiety, phase, dcd_FN, outputname, debug=False):
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
    
    # Write the configurations
    if not isfile(dcd_FN):
      self._write_dcd(dcd_FN, confs,
        includeReceptor=(moiety.find('R')>-1),
        includeLigand=(moiety.find('L')>-1))
    
    # Decompress prmtop and inpcrd files
    decompress = (self._FNs['prmtop'][moiety].endswith('.gz')) or \
                 (self._FNs['inpcrd'][moiety].endswith('.gz'))
    if decompress:
      for key in ['prmtop','inpcrd']:
        if self._FNs[key][moiety].endswith('.gz'):
          import shutil
          shutil.copy(self._FNs[key][moiety],self._FNs[key][moiety]+'.BAK')
          os.system('gunzip -f '+self._FNs[key][moiety])
          os.rename(self._FNs[key][moiety]+'.BAK', self._FNs[key][moiety])
          self._FNs[key][moiety] = self._FNs[key][moiety][:-3]
    
    # Run NAMD
    import AlGDock.NAMD
    energyCalc = AlGDock.NAMD.NAMD(\
      prmtop=self._FNs['prmtop'][moiety], \
      inpcrd=self._FNs['inpcrd'][moiety], \
      fixed={'R':self._FNs['fixed_atoms']['R'], \
             'L':None, \
             'RL':self._FNs['fixed_atoms']['RL']}[moiety], \
      solvent={'NAMD_GBSA':'GBSA', 'NAMD_Gas':'Gas'}[phase], \
      useCutoff=(phase=='NAMD_GBSA'), \
      namd_command=self._FNs['namd'])
    E = energyCalc.energies_PE(\
      outputname, dcd_FN, energyFields=[1, 2, 3, 4, 5, 6, 8, 12], \
      keepScript=debug, writeEnergyDatGZ=False)

    # Clear decompressed files
    if decompress:
      for key in ['prmtop','inpcrd']:
        if isfile(self._FNs[key][moiety]+'.gz'):
          os.remove(self._FNs[key][moiety])
          self._FNs[key][moiety] = self._FNs[key][moiety] + '.gz'

    return np.array(E, dtype=float)*MMTK.Units.kcal/MMTK.Units.mol
  
  def _write_mdcrd(self, mdcrd_FN, confs, title='', \
      includeLigand=True, includeReceptor=False, factor=1.0/MMTK.Units.Ang):
    """
    Writes an AMBER format trajectory file
    """
    n_atoms = 0
    if includeReceptor:
      receptor_0 = factor*self.confs['receptor'][:self._ligand_first_atom,:]
      receptor_1 = factor*self.confs['receptor'][self._ligand_first_atom:,:]
      n_atoms += self.confs['receptor'].shape[0]
    if includeLigand:
      n_atoms += len(self.molecule.atoms)

    if not isinstance(confs,list):
      confs = [confs]
    if includeReceptor:
      if includeLigand:
        confs = [np.vstack((receptor_0, \
          conf[self.molecule.prmtop_atom_order,:]/MMTK.Units.Ang, \
          receptor_1)) for conf in confs]
      else:
        confs = [factor*self.confs['receptor']]
    else:
      confs = [conf[self.molecule.prmtop_atom_order,:]/MMTK.Units.Ang \
        for conf in confs]
    
    import AlGDock.IO
    IO_crd = AlGDock.IO.crd()
    IO_crd.write(mdcrd_FN, confs, title, trajectory=True)
    self.tee("  wrote %d configurations to %s"%(len(confs), mdcrd_FN))

  def _read_dock6(self, mol2FN):
    # Specifically to read output from UCSF dock6
    if not isfile(mol2FN):
      raise Exception('mol2 file %s does not exist!'%mol2FN)
    if mol2FN.endswith('.mol2'):
      mol2F = open(mol2FN,'r')
    elif mol2FN.endswith('.mol2.gz'):
      mol2F = gzip.open(mol2FN,'r')
    models = mol2F.read().strip().split('########## Name:')
    mol2F.close()
    models.pop(0)
    
    crds = []
    E = {}
    for line in models[0].split('\n'):
      if line.startswith('##########'):
        label = line[11:line.find(':')].strip()
        E[label] = []
    
    for model in models:
      fields = model.split('<TRIPOS>')
      for line in fields[0].split('\n'):
        if line.startswith('##########'):
          label = line[11:line.find(':')].strip()
          E[label].append(float(line.split()[-1]))
      crds.append(np.array([l.split()[2:5] for l in fields[2].split('\n')[1:-1]],
        dtype=float)[self.molecule.inv_prmtop_atom_order,:]/10.)
    return (crds,E)

  def _write_dcd(self, dcd_file_name, confs,
      includeLigand=True, includeReceptor=False,
      factor=1.0/MMTK.Units.Ang,
      delta_t=0.1):
    """
    Writes a DCD file for a trajectory.
    If includeReceptor==True, the receptor coordinates are included.
    """
    import MMTK_DCD  # @UnresolvedImport
    from Scientific import N

    i_start = 0       # always start at frame 0
    n_savc  = 1       # save every frame
    
    if not isinstance(confs,list):
      confs = [confs]
    
    n_atoms = 0
    if includeReceptor:
      receptor_x0 = factor*self.confs['receptor'][:self._ligand_first_atom,0]
      receptor_y0 = factor*self.confs['receptor'][:self._ligand_first_atom,1]
      receptor_z0 = factor*self.confs['receptor'][:self._ligand_first_atom,2]
      receptor_x1 = factor*self.confs['receptor'][self._ligand_first_atom:,0]
      receptor_y1 = factor*self.confs['receptor'][self._ligand_first_atom:,1]
      receptor_z1 = factor*self.confs['receptor'][self._ligand_first_atom:,2]
      n_atoms += self.confs['receptor'].shape[0]
    if includeLigand:
      n_atoms += len(self.molecule.atoms)
    n_snaps = len(confs)

    fd = MMTK_DCD.writeOpenDCD(dcd_file_name, n_atoms, n_snaps,
                               i_start, n_savc, delta_t)

    if includeReceptor and includeLigand:
      for array in confs:
        array = factor*array
        x = N.concatenate((receptor_x0,N.take(array[:,0],self.molecule.prmtop_atom_order),receptor_x1)).astype(N.Float16)
        y = N.concatenate((receptor_y0,N.take(array[:,1],self.molecule.prmtop_atom_order),receptor_y1)).astype(N.Float16)
        z = N.concatenate((receptor_z0,N.take(array[:,2],self.molecule.prmtop_atom_order),receptor_z1)).astype(N.Float16)
        MMTK_DCD.writeDCDStep(fd, x, y, z)
      MMTK_DCD.writeCloseDCD(fd)
    elif includeLigand:
      for array in confs:
        array = factor*array
        x = N.take(array[:,0], self.molecule.prmtop_atom_order).astype(N.Float16)
        y = N.take(array[:,1], self.molecule.prmtop_atom_order).astype(N.Float16)
        z = N.take(array[:,2], self.molecule.prmtop_atom_order).astype(N.Float16)
        MMTK_DCD.writeDCDStep(fd, x, y, z)
      MMTK_DCD.writeCloseDCD(fd)
    else:
      x = N.concatenate((receptor_x0,receptor_x1)).astype(N.Float16)
      y = N.concatenate((receptor_y0,receptor_y1)).astype(N.Float16)
      z = N.concatenate((receptor_z0,receptor_z1)).astype(N.Float16)
      MMTK_DCD.writeDCDStep(fd, x, y, z)
      MMTK_DCD.writeCloseDCD(fd)

  def _load_pkl_gz(self, FN):
    if isfile(FN) and os.path.getsize(FN)>0:
      F = gzip.open(FN,'r')
      try:
        data = pickle.load(F)
      except:
        self.tee('  error loading '+FN)
        F.close()
        return None
      F.close()
      return data
    else:
      return None

  def _write_pkl_gz(self, FN, data):
    F = gzip.open(FN,'w')
    pickle.dump(data,F)
    F.close()
    self.tee("  wrote to "+FN)

  def _load_progress(self, p):
    progress = None
    progress_FN = join(self.dir[p],'%s_progress.pkl.gz'%p)
    for FN in [progress_FN, progress_FN+'.BAK']:
      progress = self._load_pkl_gz(FN)
      if (progress is None):
        print '  failed to load '+FN
        if isfile(FN):
          os.remove(FN)
        dirN = os.path.dirname(FN)
        if dirN!='' and (not os.path.isdir(dirN)):
          os.system('mkdir -p '+dirN)
        open(FN, 'a').close()
      else:
        print '  loaded '+FN
        break
    if progress is not None:
      params = progress[0]
      setattr(self,'%s_protocol'%p,progress[1])
      setattr(self,'_%s_cycle'%p,progress[2])
      setattr(self,'_%s_total_cycle'%p,progress[3])
      if p=='dock':
        (self._n_trans, self._max_n_trans, self._random_trans, \
         self._n_rot, self._max_n_rot, self._random_rotT) = progress[4]
      self.confs[p]['replicas'] = progress[5]
      self.confs[p]['seeds'] = progress[6]
      self.confs[p]['samples'] = progress[7]
      setattr(self,'%s_Es'%p,progress[8])
    else:
      params = None
      setattr(self,'%s_protocol'%p,[])
      setattr(self,'_%s_cycle'%p,0)
      setattr(self,'_%s_total_cycle'%p,0)
      self.confs[p]['replicas'] = None
      self.confs[p]['seeds'] = None
      self.confs[p]['samples'] = None
      setattr(self,'%s_Es'%p,None)
    return params

  def _save_progress(self, p):
    """
    Saves the protocol, 
    cycle counts,
    random orientation parameters (for docking),
    replica configurations,
    sampled configurations,
    and energies
    """
    if p=='dock':
      random_orient = (self._n_trans, self._max_n_trans, self._random_trans, \
         self._n_rot, self._max_n_rot, self._random_rotT)
    else:
      random_orient = None
    
    arg_dict = dict([tp for tp in self.params[p].items() \
                      if not tp[0] in ['repX_cycles']])
    if p=='cool':
      fn_dict = convert_dictionary_relpath({
          'ligand_database':self._FNs['ligand_database'],
          'forcefield':self._FNs['forcefield'],
          'frcmodList':self._FNs['frcmodList'],
          'prmtop':{'L':self._FNs['prmtop']['L']},
          'inpcrd':{'L':self._FNs['inpcrd']['L']}},
          relpath_o=None, relpath_n=self.dir['cool'])
      params = (fn_dict,arg_dict)
    elif p=='dock':
      fn_dict = convert_dictionary_relpath(
          dict([tp for tp in self._FNs.items() \
            if not tp[0] in ['namd','vmd','sander']]),
          relpath_o=None, relpath_n=self.dir['dock'])
    params = (fn_dict,arg_dict)
    
    progress = (params,
                getattr(self,'%s_protocol'%p),
                getattr(self,'_%s_cycle'%p),
                getattr(self,'_%s_total_cycle'%p),
                random_orient,
                self.confs[p]['replicas'],
                self.confs[p]['seeds'],
                self.confs[p]['samples'],
                getattr(self,'%s_Es'%p))
    progress_FN = join(self.dir[p],'%s_progress.pkl.gz'%p)
    if isfile(progress_FN):
      os.rename(progress_FN,progress_FN+'.BAK')
    self._write_pkl_gz(progress_FN, progress)

  def _set_lock(self, p):
    if not os.path.isdir(self.dir[p]):
      os.system('mkdir -p'+self.dir[p])
    lockFN = join(self.dir[p],'.lock')
    if isfile(lockFN):
      raise Exception(p + ' is locked')
    else:
      lockF = open(lockFN,'w')
      lockF.close()
    logFN = join(self.dir[p],p+'_log.txt')
    self.log = open(logFN,'a')

  def _clear_lock(self, p):
    lockFN = join(self.dir[p],'.lock')
    if isfile(lockFN):
      os.remove(lockFN)
    self.log.close()
    del self.log

  def tee(self, var, process=None):
    print var
    if hasattr(self,'log'):
      if isinstance(var,str):
        self.log.write(var+'\n')
      else:
        self.log.write(repr(var)+'\n')
    elif process is not None:
      self._set_lock(process)
      if isinstance(var,str):
        self.log.write(var+'\n')
      else:
        self.log.write(repr(var)+'\n')
      self._clear_lock(process)

if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(
    description='Molecular docking with adaptively scaled alchemical interaction grids')
  
  # Files and directories
  #   Run-dependent
  parser.add_argument('--dir_dock', 
    help='Directory where docking results are stored')
  parser.add_argument('--dir_cool', 
    help='Directory where cooling results are stored')
  parser.add_argument('--namd',
    help='Location of Not Another Molecular Dynamics (NAMD)')
  parser.add_argument('--vmd',
    help='Location of Visual Molecular Dynamics (VMD)')
  parser.add_argument('--sander',
    help='Location of sander (from AMBER)')
  #   Stored in both dir_cool and dir_dock
  parser.add_argument('--ligand_database', 
    help='MMTK molecule definition for ligand')
  parser.add_argument('--forcefield', 
    help='AMBER force field file')
  parser.add_argument('--frcmodList', nargs='+',
    help='AMBER force field modifications file(s)')
  parser.add_argument('--ligand_prmtop', 
    help='AMBER prmtop for the ligand')
  parser.add_argument('--ligand_inpcrd', 
    help='AMBER coordinates for the ligand')
  #   Stored in dir_dock
  parser.add_argument('--receptor_prmtop', 
    help='AMBER prmtop for the receptor')
  parser.add_argument('--receptor_inpcrd', 
    help='AMBER coordinates for the receptor')
  parser.add_argument('--receptor_fixed_atoms',
    help='PDB file with fixed atoms labeled by 1 in the occupancy column')
  parser.add_argument('--complex_prmtop', 
    help='AMBER prmtop file for the complex')
  parser.add_argument('--complex_inpcrd', 
    help='AMBER coordinates for the complex')
  parser.add_argument('--complex_fixed_atoms',
    help='PDB file with fixed atoms labeled by 1 in the occupancy column')
  parser.add_argument('--dir_grid',
    help='Directory containing potential energy grids')
  parser.add_argument('--grid_LJr', 
    help='DX file for Lennard-Jones repulsive grid')
  parser.add_argument('--grid_LJa', 
    help='DX file for Lennard-Jones attractive grid')
  parser.add_argument('--grid_ELE', 
    help='DX file for electrostatic grid')
  # Simulation settings and constants
  #   Run-dependent 
  parser.add_argument('--cool_repX_cycles', type=int,
    help='Number of replica exchange cycles for cooling')
  parser.add_argument('--dock_repX_cycles', type=int,
    help='Number of replica exchange cycles for docking')
  parser.add_argument('--run_type',
    choices=['pose_energies','one_step',\
             'initial_cool','cool','random_dock',\
             'initial_dock','dock','all', \
             'store_params','free_energies', 'postprocess',
             'redo_postprocess', 'clear_intermediates', None],
    help='Type of calculation to run')
  parser.add_argument('--nodes', type=int,
    help='Number of CPU cores to use')
  #   Defaults
  parser.add_argument('--protocol', choices=['Adaptive','Set'],
    help='Approach to determining series of thermodynamic states')
  parser.add_argument('--no_protocol_refinement', action='store_true', default=None,
    help='Does not refine the protocol during replica exchange')
  parser.add_argument('--therm_speed', type=float,
    help='Thermodynamic speed during adaptive simulation')
  parser.add_argument('--sampler',
    choices=['HMC','NUTS','VV'],
    help='Sampling method')
  parser.add_argument('--MCMC_moves', type=int,
    help='Types of MCMC moves to use')
  # For initialization
  parser.add_argument('--seeds_per_state', type=int,
    help='Number of starting configurations in each state during initialization')
  parser.add_argument('--steps_per_seed', type=int,
    help='Number of MD steps per state during initialization')
  # For replica exchange
  parser.add_argument('--repX_cycles', type=int,
    help='Number of replica exchange cycles for docking and cooling')
  parser.add_argument('--sweeps_per_cycle', type=int,
    help='Number of replica exchange sweeps per cycle')
  parser.add_argument('--steps_per_sweep', type=int,
    help='Number of MD steps per replica exchange sweep')
  parser.add_argument('--keep_intermediate', action='store_true', default=None,
    help='Keep configurations for intermediate states?')
  parser.add_argument('--min_repX_acc', type=float,
    help='Minimum value for replica exchange acceptance rate')
  # For postprocessing
  parser.add_argument('--phases', nargs='+', \
    help='Phases to use in postprocessing')
  #   Stored in dir_cool
  parser.add_argument('--cool_protocol', choices=['Adaptive','Set'],
    help='Approach to determining series of thermodynamic states')
  parser.add_argument('--cool_therm_speed', type=float,
    help='Thermodynamic speed during adaptive simulation')
  parser.add_argument('--cool_sampler',
    choices=['HMC','NUTS','VV'],
    help='Sampling method')
  # For initialization
  parser.add_argument('--cool_seeds_per_state', type=int,
    help='Number of starting configurations in each state during initialization')
  parser.add_argument('--cool_steps_per_seed', type=int,
    help='Number of MD steps per state during initialization')
  # For replica exchange
  parser.add_argument('--cool_sweeps_per_cycle', type=int,
    help='Number of replica exchange sweeps per cycle')
  parser.add_argument('--cool_steps_per_sweep', type=int,
    help='Number of MD steps per replica exchange sweep')
  parser.add_argument('--cool_keep_intermediate', action='store_true',
    default=None, help='Keep configurations for intermediate states?')
  #   Stored in dir_dock
  parser.add_argument('--dock_protocol', choices=['Adaptive','Set'],
    help='Approach to determining series of thermodynamic states')
  parser.add_argument('--dock_therm_speed', type=float,
    help='Thermodynamic speed during adaptive simulation')
  parser.add_argument('--dock_sampler',
    choices=['HMC','NUTS','VV'],
    help='Sampling method')
  # For initialization
  parser.add_argument('--dock_seeds_per_state', type=int,
    help='Number of starting configurations in each state during initialization')
  parser.add_argument('--dock_steps_per_seed', type=int,
    help='Number of MD steps per state during initialization')
  # For replica exchange
  parser.add_argument('--dock_sweeps_per_cycle', type=int,
    help='Number of replica exchange sweeps per cycle')
  parser.add_argument('--dock_steps_per_sweep', type=int,
    help='Number of MD steps per replica exchange sweep')
  parser.add_argument('--dock_keep_intermediate', action='store_true',
    default=None, help='Keep configurations for intermediate states?')
  parser.add_argument('--score', nargs='?', const=True, default=False,
    help='Start replica exchange from a configuration or set of configurations, which may be passed as an argument. The default configuration is from the ligand_inpcrd argument.')
  parser.add_argument('--rmsd', nargs='?', const=True, default=False,
    help='Calculate rmsd between snapshots and a configuration or set of configurations, which may be passed as an argument. The default configuration is from the ligand_inpcrd argument.')
  # Binding site
  parser.add_argument('--site',
    choices=['Sphere','Cylinder'], help='Type of binding site')
  parser.add_argument('--site_center', nargs=3, type=float,
    help='Position of binding site center')
  parser.add_argument('--site_direction', nargs=3, type=float,
    help='Principal axis of a cylindrical binding site')
  parser.add_argument('--site_max_X', type=float,
    help='Maximum position along principal axis in a cylindrical binding site')
  parser.add_argument('--site_max_R', type=float,
    help='Maximum radial position for a spherical or cylindrical binding site')
  parser.add_argument('--site_density', type=float,
    help='Density of center-of-mass points in the first docking state')
  # Additional calculations
  parser.add_argument('--do_calc_random_dock_stats', action='store_true',
    default=None,
    help='Calculate convergence statistics for random placement of the ligand into the active site')
  parser.add_argument('--receptor_Gas', type=float, nargs='+',
    help='Receptor potential energies in AMBER Gas implicit solvent (in units of kJ/mol)')
  parser.add_argument('--receptor_GBSA', type=float, nargs='+',
    help='Receptor potential energies in AMBER GBSA implicit solvent (in units of kJ/mol)')
  parser.add_argument('--receptor_NAMD_Gas', type=float, nargs='+',
    help='Receptor potential energies in gas phase (in units of kJ/mol)')
  parser.add_argument('--receptor_NAMD_GBSA', type=float, nargs='+',
    help='Receptor potential energies in NAMD GBSA implicit solvent (in units of kJ/mol)')
  
  args = parser.parse_args()
  self = BindingPMF(**vars(args))
