#!/usr/bin/env python

# TO DO: Revamp functions to
# -compare different Lennard-Jones grid transformations

# TO DO: Keep track of memory usage and availability
# TO DO: Rewrite initial_dock to be less memory-intensive
# TO DO: Write docstrings
# TO DO: Set up checkpoints

# Ideas from Prodock by Trosset & Scheraga
# TO DO: Grid interpolation with Bezier splines,
# which provides second derivatives
# [Oberlin, D.,Jr.; Scheraga, H. A. J Comput Chem 1998, 19, 71:85.].
# This could allow for layer sampling.

# TO DO: Scaled collective variable Monte Carlo
# for dihedral angle moves
# [Noguti, T.; Go, N. Biopolymers 1985, 24, 527:546.]

# Ideas from QXP by McMartin and Bohacek
# TO DO: Dihedral angle moves using cylindrical coordinates
# (makes energy evalulations faster?)

# Ideas from Ribodock by Morley and Afshar
# TO DO: Scale number of MC trials by number of rotatable bonds

# TO DO: Genetic Monte Carlo

import os
from os.path import join
from os.path import exists
from os.path import abspath
import pickle
import gzip
import copy

import time
import numpy as np

import MMTK
import MMTK.Units
from MMTK.ParticleProperties import Configuration

import Scientific
from Scientific_vector import Vector  # @UnresolvedImport

import AlGDock as a

# Constants
R = 8.3144621*MMTK.Units.J/MMTK.Units.mol/MMTK.Units.K
T_HIGH = 800.0*MMTK.Units.K
T_TARGET = 300.0*MMTK.Units.K
RT_HIGH = R*T_HIGH
RT_TARGET = R*T_TARGET

term_map = {
  'cosine dihedral angle':'MM',
  'electrostatic/pair sum':'MM',
  'harmonic bond':'MM',
  'harmonic bond angle':'MM',
  'Lennard-Jones':'MM',
  'site':'site',
  'sLJr':'sLJr',
  'sLJa':'sLJa',
  'LJr':'LJr',
  'LJa':'LJa',
  'ELE':'ELE',
  'electrostatic':'misc'}

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

def choice(N, weights=None):
  """
  Weighted random number selection of N values 
  """
  weights /= np.sum(weights)
  cumweights = np.cumsum(weights)
  beginIndicies = np.zeros(N, dtype=int)
  randVals = np.random.uniform(size=N)

  for (bI_idx,randVal) in zip(range(N),randVals):
    beginIndicies[bI_idx] = sum(cumweights < randVal)
  beginIndicies.sort()
  return beginIndicies

def merge_dictionaries(dicts, required_consistency=[], print_inconsistent=False):
  """
  Merges a list of dictionaries, giving priority to items in descending order.
  Items in the required_consistency list must be consistent with one another.
  """
  merged = {}
  for a in range(len(dicts)): # Loop over all dictionaries, giving priority to the first
    if not isinstance(dicts[a],dict):
      continue
    for key in dicts[a].keys():
      if dicts[a][key] is None:
        if not key in merged.keys():
          merged[key] = None
      elif isinstance(dicts[a][key],dict):
        merged[key] = merge_dictionaries(
          [dicts[n][key] for n in range(len(dicts)) if key in dicts[n].keys()],
          required_consistency=required_consistency)
      else:
        # Check for consistency with other dictionaries
        for b in (range(a) + range(a+1,len(dicts))):
          if isinstance(dicts[b],dict) and (key in dicts[b].keys()) and (dicts[b][key] is not None):
            if (isinstance(dicts[b][key],np.ndarray)):
              inconsistent_items = (dicts[b][key]!=dicts[a][key]).any()
            else:
              inconsistent_items = (dicts[b][key]!=dicts[a][key])
            if inconsistent_items:
              if print_inconsistent:
                print 'Dictionary items are inconsistent:'
                print dicts[a][key]
                print dicts[b][key]
              if key in required_consistency:
                raise Exception('Items must be consistent!')
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

class NullDevice():
  """
  A device to suppress output
  """
  def write(self, s):
    pass

##############
# Main Class #
##############

class Docker:
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
      #   Stored in both dir_cool and dir_dock
      ligand_database=None,
      forcefield=None, frcmodList=None,
      ligand_prmtop=None, ligand_inpcrd=None,
      #   Stored in dir_dock
      receptor_prmtop=None, receptor_inpcrd=None,
      receptor_fixed_atoms=None,
      complex_prmtop=None, complex_inpcrd=None,
      complex_fixed_atoms=None,
      dir_grid='', grid_LJr=None, grid_LJa=None, grid_ELE=None,
      grid_force=None, # Optional, only for NAMD grids
      # Arguments - simulation parameters and constants
      #   Run-dependent
      cool_repX_cycles=None,
      dock_repX_cycles=None,
      run_type=None,
      #   Defaults for dir_cool and dir_dock
      protocol=None, therm_speed=None,
      sampler=None,
      seeds_per_state=None, steps_per_seed=None,
      repX_cycles=None,
      sweeps_per_cycle=None, steps_per_sweep=None,
      keep_intermediate=None,
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
      legs=None,
      site=None, site_center=None, site_direction=None, # Site parameters
      site_max_X=None, site_max_R=None, # Site parameters
      site_density=None,
      do_calc_random_dock_stats=None,
      rmsd_ref=None, # RMSD calculation parameters
      receptor_GBSA=None, receptor_Gas=None): # Energy values
    """Parses the input arguments and runs the requested docking calculation"""
    
    print """
###########
# ALGDock #
###########
version {0}

Molecular docking with adaptively scaled alchemical interaction grids
    """.format(a.__version__)
    
    print '*** Directories ***'
    self.dir = {}
    if dir_dock is not None:
      self.dir['dock'] = abspath(dir_dock)
    else:
      self.dir['dock'] = abspath('.')
    
    if dir_cool is not None:
      self.dir['cool'] = abspath(dir_cool)
    else:
      self.dir['cool'] = self.dir['dock']

    for p in ['cool','dock']:
      if not os.path.isdir(self.dir[p]):
        os.makedirs(self.dir[p])

    print self.dir

    # Load previously stored file names and arguments
    FNs = {}
    args = {}
    for p in ['cool','dock']:
      store_FN = join(self.dir[p],p+'_store.pkl.gz')
      if exists(store_FN):
        (fn_dict, arg_dict) = self._load_pkl_gz(store_FN)
        FNs[p] = convert_dictionary_relpath(fn_dict,
          relpath_o=self.dir[p], relpath_n=None)
        args[p] = arg_dict
      else:
        FNs[p] = {}
        args[p] = {}
  
    # Set up file name dictionary
    print '\n*** Files ***'

    for p in ['cool','dock']:
      if p in FNs.keys():
        if FNs[p]!={}:
          print 'previous stored in %s directory:'%p
          print FNs[p]

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
                               join(dir_grid,'LJr.dx'),
                               join(dir_grid,'LJr.dx.gz')]),
               'LJa':a.findPath([grid_LJa,
                               join(dir_grid,'LJa.dx'),
                               join(dir_grid,'LJa.dx.gz')]),
               'ELE':a.findPath([grid_ELE,
                               join(dir_grid,'electrostatic.dx'),
                               join(dir_grid,'electrostatic.dx.gz'),
                               join(dir_grid,'pbsa.dx'),
                               join(dir_grid,'pbsa.dx.gz')])},
      'grid_force':a.findPath([grid_force]),
      'namd':a.findPath([namd] + a.search_paths['namd']),
      'vmd':a.findPath([vmd] + a.search_paths['vmd'])}

    if not (FNs['cool']=={} and FNs['dock']=={}):
      print 'from arguments and defaults:'
      print FNs['new']
      print '\nto be used:'

    self._FNs = merge_dictionaries(
      [FNs[src] for src in ['new','cool','dock']],
      required_consistency=['L','R','RL','ligand_database'],
      print_inconsistent=True)
  
    # Default: a force field modification is in the same directory as the ligand
    if (self._FNs['frcmodList'] is None):
      dir_ligand_db = os.path.dirname(self._FNs['ligand_database'])
      frcmod = a.findPath([abspath(join(dir_ligand_db,'lig.frcmod')),
                         abspath(join(dir_ligand_db,'ligand.frcmod'))])
      if frcmod is not None:
        self._FNs['frcmodList'] = [frcmod]

    print self._FNs

    print '\n*** Simulation parameters and constants ***'
    
    args['default_cool'] = {
        'protocol':'Adaptive',
        'therm_speed':2.0,
        'sampler':'HMC',
        'seeds_per_state':50,
        'steps_per_seed':20000,
        'repX_cycles':20,
        'sweeps_per_cycle':100,
        'steps_per_sweep':2500,
        'keep_intermediate':False}

    args['default_dock'] = dict(args['default_cool'].items() + {
      'legs':1,
      'site':None, 'site_center':None, 'site_direction':None,
      'site_max_X':None, 'site_max_R':None,
      'site_density':500.,
      'do_calc_random_dock_stats':False,
      'rmsd_ref':{'vmdcrd':None,
                  'species':'RL',
                  'selection':'resname iga'},
      'receptor_GBSA':None, 'receptor_Gas':None}.items())

    # Store passed arguments in dictionary
    if rmsd_ref is not None:
      rmsd_ref_dict = {}
      if isinstance(rmsd_ref,str):
        rmsd_ref = [rmsd_ref]
      rmsd_ref[0] = os.path.relpath(rmsd_ref[0],self.dir['dock'])
      for ind in range(len(rmsd_ref)):
        rmsd_ref_dict[['vmdcrd','species','selection'][ind]] = rmsd_ref[ind]
      rmsd_ref = rmsd_ref_dict
  
    namespace = locals()
    for p in ['cool','dock']:
      args['new_'+p] = {}
      for key in args['default_'+p].keys():
        specific_key = p + '_' + key
        if (specific_key in namespace.keys()) and \
           (namespace[specific_key] is not None):
          # Use the specific key if it exists
          args['new_'+p][key] = namespace[specific_key]
        elif (key in ['site_center', 'site_direction',
                      'receptor_Gas','receptor_GBSA']) and \
             (namespace[key] is not None):
          # Convert these to arrays of floats
          args['new_'+p][key] = np.array(namespace[key], dtype=float)
        else:
          # Use the general key
          args['new_'+p][key] = namespace[key]

    self.params = {}
    for p in ['cool','dock']:
      self.params[p] = merge_dictionaries(
        [args[src] for src in ['default_'+p,'new_'+p,p]])

    for p in ['cool','dock']:
      print 'for %s:'%p
      print self.params[p]

    # Variables dependent on the parameters
    self.original_Es = [[{}]]
    for phase in ['Gas','GBSA']:
      if self.params['dock']['receptor_'+phase] is not None:
        self.original_Es[0][0]['R'+phase] = \
          self.params['dock']['receptor_'+phase]
      
    if self.params['dock']['legs']==4:
      self._scalables = ['sLJr','sLJa','LJr','LJa','ELE']
    else:
      self._scalables = ['LJr','LJa','ELE']
  
    self.lambda_full = {'T':T_HIGH,'MM':True,'site':True}
    self.lambda_scalables = {'T':T_HIGH}
    for scalable in self._scalables:
        self.lambda_full[scalable] = 1
        self.lambda_scalables[scalable] = 1

    # Check for existence of required files
    for FN in [self._FNs['ligand_database'],
               self._FNs['forcefield'],
               self._FNs['prmtop']['L'],
               self._FNs['inpcrd']['L']]:
      if (FN is None) or (not exists(FN)):
        raise Exception('Required file %s is missing!'%FN)

    for FN in [self._FNs['prmtop']['RL'],
               self._FNs['inpcrd']['RL'],
               self._FNs['fixed_atoms']['RL'],
               self._FNs['grids']['LJr'],
               self._FNs['grids']['LJa'],
               self._FNs['grids']['ELE']]:
      if (FN is None) or (not exists(FN)):
        print 'Missing file %s is necessary for docking!'%FN
        if ('run_type' is args) and \
           (args.run_type is not None) and \
           (args.run_type in ['Dock0_0','ToDock0_0','Dock','ToDock']):
          raise Exception('File %s is required!'%FN)
  
    if (self._FNs['inpcrd']['RL'] is None) and \
       (self._FNs['inpcrd']['R'] is None):
      raise Exception('Receptor coordinates needed')

    self._setup_universe()

    if run_type=='initial_cool':
      self.initial_cool()
    elif run_type=='cool':
      self.cool()
    elif run_type=='initial_dock':
      self.initial_dock()
    elif run_type=='dock':
      self.dock()
    elif run_type=='all':
      self.cool()
      self.calc_f_L()
      self.dock()
      self.calc_f_RL()
      self.calc_RMSD()

  def _setup_universe(self):
    """Creates an MMTK InfiniteUniverse and adds the ligand"""

    MMTK.Database.molecule_types.directory = \
      os.path.dirname(self._FNs['ligand_database'])

    # Force fields
    from MMTK.ForceFields import Amber12SBForceField

    self._forceFields = {}
    self._forceFields['gaff'] = Amber12SBForceField(
      parameter_file=self._FNs['forcefield'],mod_files=self._FNs['frcmodList'])

    if self.params['dock']['site']=='Sphere':
      from AlGDock.ForceFields.Sphere.Sphere import SphereForceField
      self._forceFields['site'] = SphereForceField(
        center=self.params['dock']['site_center'],
        max_R=self.params['dock']['site_max_R'], name='site')
    elif self.params['dock']['site']=='Cylinder':
      from AlGDock.ForceFields.Cylinder.Cylinder import CylinderForceField
      self._forceFields['site'] = CylinderForceField(
        origin=self.params['dock']['site_center'],
        direction=self.params['dock']['site_direction'],
        max_X=self.params['dock']['site_max_X'],
        max_R=self.params['dock']['site_max_R'], name='site')
    else:
      raise Exception('Binding site type not recognized!')
  
    # Set up the system
    import sys
    original_stderr = sys.stderr
    sys.stderr = NullDevice()

    self.molecule = MMTK.Molecule(\
      os.path.basename(self._FNs['ligand_database']))
    
    sys.stderr = original_stderr
    
    # Randomly rotate the molecule and translate it into the binding site
    from Scientific.Geometry.Transformation import Rotation
    self.molecule.applyTransformation(Rotation(random_rotate()))
    self.molecule.translateTo(Vector(self._forceFields['site'].randomPoint()))
    
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
    self.sampler['init'] = VelocityVerletIntegrator(self.universe)
    for p in ['cool', 'dock']:
      if self.params[p]['sampler'] == 'HMC':
        self.sampler[p] = HamiltonianMonteCarloIntegrator(self.universe)
      elif self.params[p]['sampler'] == 'NUTS':
        self.sampler[p] = NUTSIntegrator(self.universe)
      elif self.params[p]['sampler'] == 'VV':
        self.sampler[p] = VelocityVerletIntegrator(self.universe)
      else:
        raise Exception('Unrecognized sampler!')

    prmtop_atom_order = np.array([atom.number \
      for atom in self.molecule.prmtop_order], dtype=int)
    inv_prmtop_atom_order = np.zeros(shape=len(prmtop_atom_order), dtype=int)
    for i in range(len(prmtop_atom_order)):
      inv_prmtop_atom_order[prmtop_atom_order[i]] = i

    # Determine ligand atomic index
    def load_prmtop(prmtopFN, varnames=['RESIDUE_LABEL','RESIDUE_POINTER']):
      """
      Loads records from an AMBER prmtop file and stores data in a dictionary
      """
      def load_record(record):
        items = []
        lines = record.split('\n')
        lines.pop(0) # Name
        FORMAT = lines.pop(0).strip()[8:-1] # Format
        if FORMAT.find('a')>-1: # Text
          w = int(FORMAT[FORMAT.find('a')+1:])
          for line in lines:
            items = items + [line[x:x+w] for x in range(0,len(line),w)]
          return np.array(items)
        elif FORMAT.find('I')>-1: # Integer
          w = int(FORMAT[FORMAT.find('I')+1:])
          for line in lines:
            items = items + [int(line[x:x+w]) for x in range(0,len(line),w)]
          return np.array(items, dtype=int)
        elif FORMAT.find('E')>-1: # Scientific
          w = int(FORMAT[FORMAT.find('E')+1:FORMAT.find('.')])
          for line in lines:
            items = items + [float(line[x:x+w]) for x in range(0,len(line),w)]
          return np.array(items, dtype=float)

      prmtopF = open(prmtopFN,'r')
      prmtopData = prmtopF.read().split('%FLAG ')
      prmtopF.close()
      
      prmtop = {}
      for record in prmtopData:
        name = record[:record.find('\n')].strip()
        if name in varnames:
          prmtop[name] = load_record(record)
      return prmtop
    
    if (self._FNs['prmtop']['R'] is not None) and \
       (self._FNs['prmtop']['RL'] is not None):
      prmtop_R = load_prmtop(self._FNs['prmtop']['R'])
      prmtop_RL = load_prmtop(self._FNs['prmtop']['RL'])
      ligand_ind = [ind for (r,ind) in \
        zip(prmtop_RL['RESIDUE_LABEL'],range(len(prmtop_RL['RESIDUE_LABEL']))) \
        if r not in prmtop_R['RESIDUE_LABEL']]
      if len(ligand_ind) > 1:
        raise Exception('Ligand residue label is ambiguous')
      self._ligand_first_atom = prmtop_RL['RESIDUE_POINTER'][ligand_ind[0]] - 1
    else:
      raise Exception('Missing AMBER prmtop files')

    # Read the receptor coordinates
    if self._FNs['inpcrd']['R'] is not None:
      self._receptor_conf = self._read_crd(self._FNs['inpcrd']['R'],split=False)
    else:
      complex_crd = self._read_crd(self._FNs['inpcrd']['RL'], split=False)
      self._receptor_conf = np.vstack(\
        (complex_crd[:self._ligand_first_atom,:],\
         complex_crd[self._ligand_first_atom + self._ligand_natoms:,:]))

    # Load progress
    self.confs = {}
    self._load_progress('cool')
    self._load_progress('dock')
  
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
    
    The following are saved:
    (1) the cooling protocol (cool_protocol.pkl.gz)
    (2) starting configurations for replica exchange (cool_replicas.crd.gz)
    """

    if (len(self.cool_protocol)>0) and (self.cool_protocol[-1]['crossed']):
      return # Initial cooling is already complete
    
    self._set_lock('cool')
    print "\n>>> Initial cooling of the ligand "+\
          "from %d K to %d K"%(T_HIGH,T_TARGET)
    cool_start_time = time.time()

    # Save the simulation parameters
    self._write_pkl_gz(join(self.dir['cool'],'cool_store.pkl.gz'),
      (convert_dictionary_relpath({
        'ligand_database':self._FNs['ligand_database'],
        'forcefield':self._FNs['forcefield'],
        'frcmodList':self._FNs['frcmodList'],
        'prmtop':{'L':self._FNs['prmtop']['L']},
        'inpcrd':{'L':self._FNs['inpcrd']['L']}},
        relpath_o=None, relpath_n=self.dir['cool']),
       dict([tp for tp in self.params['cool'].items() \
          if not tp[0] in ['repX_cycles']])))

    # Set up the force field
    self.cool_protocol = [{'MM':True, 'T':T_HIGH, \
                          'delta_t':1.5*MMTK.Units.fs, 'crossed':False}]
    self._set_universe_force_field(self.cool_protocol[-1])

    # Minimize and ramp the temperature from 0 to the desired high temperature
    from MMTK.Minimization import ConjugateGradientMinimizer # @UnresolvedImport
    minimizer = ConjugateGradientMinimizer(self.universe)
    minimizer(steps = 500)

    for T in np.linspace(0.,T_HIGH,33)[1:]:
      self.sampler['init'](steps = 500, T=T,\
                           delta_t=self.delta_t, steps_per_trial = 100)

    # Run at high temperature
    state_start_time = time.time()
    conf = self.universe.configuration().array
    (confs, Es_MM, g) = self._initial_sim_state('cool', \
      [conf for n in range(self.params['cool']['seeds_per_state'])], \
      normalize=True)
    self.confs['cool_replicas'] = [confs[np.argmin(Es_MM)]]
    self.cool_protocol[-1]['delta_t'] = self.delta_t
    self.cool_protocol[-1]['g'] = g
    print "  generated %d configurations "%len(confs) + \
          "at %d K "%self.cool_protocol[-1]['T'] + \
          "in %f s "%(time.time()-state_start_time) + \
          "(g=%f, "%self.cool_protocol[-1]['g'] + \
          "Delta t=%f)"%self.cool_protocol[-1]['delta_t']

    # Main loop for initial cooling:
    # choose new temperature, randomly select seeds, simulate
    T = T_HIGH
    while (not self.cool_protocol[-1]['crossed']):
      # Choose new temperature
      To = self.cool_protocol[-1]['T']
      dL = Es_MM.std()/(R*T*T)
      if dL>0:
        T = T - self.params['cool']['therm_speed']/dL
        if T < T_TARGET:
          lambda_n = {'T':T_TARGET, 'MM':True, 'crossed':True}
        else:
          lambda_n = {'T':T, 'MM':True, 'crossed':False}
      else:
        lambda_n = {'T':T_TARGET, 'MM':True, 'crossed':True}
      self.cool_protocol.append(lambda_n)
      self._set_universe_force_field(self.cool_protocol[-1])

      # Randomly select seeds for new trajectory
      seedIndicies = choice(self.params['cool']['seeds_per_state'], \
                             np.exp(-Es_MM/R*(1/T-1/To)))

      # Simulate and store data
      state_start_time = time.time()
      (confs, Es_MM, g) = self._initial_sim_state('cool', \
        [confs[ind] for ind in seedIndicies], normalize=True)
      self.confs['cool_replicas'].append(confs[np.argmin(Es_MM)])
      self.cool_protocol[-1]['delta_t'] = self.delta_t
      self.cool_protocol[-1]['g'] = g
      print "  generated %d configurations "%len(confs) + \
            "at %d K "%self.cool_protocol[-1]['T'] + \
            "in %f s "%(time.time()-state_start_time) + \
            "(g=%f, "%self.cool_protocol[-1]['g'] + \
            "Delta t=%f)"%self.cool_protocol[-1]['delta_t']

    self.cool_Es = [[] for s in range(len(self.cool_protocol))]

    # Save data
    self._write_crd(join(self.dir['cool'],'cool_replicas.crd.gz'),
      self.confs['cool_replicas'], title='cool_replicas')
    self._write_pkl_gz(join(self.dir['cool'],'cool_protocol.pkl.gz'), \
      self.cool_protocol)

    self._clear_lock('cool')
    print "\nElapsed time for initial cooling: "+\
          "%f s"%(time.time()-cool_start_time)

  def cool(self):
    """
    Samples different ligand configurations 
    at thermodynamic states between T_HIGH and T_TARGET
    """
    self._sim_process('cool')

  def calc_f_L(self, readOnly=False):
    """
    Calculates ligand-specific free energies:
    1. solvation free energy of the ligand using single-step 
       free energy perturbation
    2. reduced free energy of cooling the ligand from T_HIGH to T_TARGET
    """
    
    if self.cool_protocol==[]:
      return # Initial cooling is incomplete
    
    # Initialize variables by loading data or as empty lists
    f_L_FN = join(self.dir['cool'],'f_L.pkl.gz')
    if exists(f_L_FN):
      (self.f_L_solvation, self.f_cool, self.cool_equilibrated_cycle, \
       self.cool_mean_acc) = self._load_pkl_gz(f_L_FN)
    else:
      for var in ['f_L_solvation', 'f_cool', \
                  'cool_equilibrated_cycle', 'cool_mean_acc']:
        if not hasattr(self,var):
          setattr(self,var,[])
    if readOnly:
      return

    self._postprocess([('cool',-1,-1,'L')])

    # Estimate cycle at which simulation has equilibrated
    K = len(self.cool_protocol)
    u_Ks = [self._u_kln(\
      [self.cool_Es[K-1][c]],\
      [self.cool_protocol[K-1]]) for c in range(self._cool_cycle)]
    mean_u_Ks = np.array([np.mean(u_K) for u_K in u_Ks])
    std_u_Ks = np.array([np.std(u_K) for u_K in u_Ks])
    for c in range(len(self.cool_equilibrated_cycle), self._cool_cycle):
      nearMean = (mean_u_Ks - mean_u_Ks[c])<std_u_Ks[c]
      for (item, ind) in zip(nearMean, range(len(nearMean))):
        if item:
          self.cool_equilibrated_cycle.append(ind)
          break

    # Calculate cooling free energies that have not already been calculated
    updated = False
    for c in range(len(self.f_L_solvation), self._cool_cycle):
      fromCycle = self.cool_equilibrated_cycle[c]
      toCycle = c + 1
      
      u_GBSA = np.concatenate(\
        [self.cool_Es[K-1][c]['LGBSA'][:,-1] for c in range(fromCycle,toCycle)])
      u_Gas = np.concatenate(\
        [self.cool_Es[K-1][c]['LGas'][:,-1] for c in range(fromCycle,toCycle)])
      du = (u_GBSA - u_Gas)/RT_TARGET
      min_du = min(du)
      self.f_L_solvation.append(-np.log(np.exp(-du+min_du).mean()) + min_du)

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

      print "  calculated cooling free energy of %f kT "%(self.f_cool[-1][-1])+\
            "using cycles %d to %d"%(fromCycle, toCycle-1)
      updated = True

    if updated:
      self._write_pkl_gz(f_L_FN, \
        (self.f_L_solvation, self.f_cool, \
         self.cool_equilibrated_cycle, self.cool_mean_acc))

  ###########
  # Docking #
  ###########
  def initial_dock(self):
    """
      Docks the ligand into the receptor
      
      The first state of docking is sampled by randomly placing configurations 
      from the high temperature ligand simulation into the binding site.
      Intermediate thermodynamic states are chosen such that
      thermodynamic length intervals are approximately constant.
      Configurations from each state are subsampled to seed the next simulation.
      
      The following are saved:
      (1) the docking protocol (dock_protocol.pkl.gz)
      (2) starting configurations for replica exchange (dock_replicas.crd.gz)
    """
    
    if self.cool_protocol==[]:
      raise Exception('Cannot perform initial docking '+\
                      'until cooling is complete')
    
    if (len(self.dock_protocol)>0) and (self.dock_protocol[-1]['crossed']):
      return # Initial docking already complete
    
    # Load the first cooling state and make sure there are enough samples
    print "\n>>> Loading the first state of cooling"
    (cool0_confs, cool0_Es_MM) = self._load_cool0()
    while len(cool0_Es_MM)<self.params['dock']['seeds_per_state']:
      print "More samples from high temperature ligand simulation needed"
      self._replica_exchange('cool')
      (cool0_confs, cool0_Es_MM) = self._load_cool0()
    print "  using %d configurations "%len(cool0_Es_MM) +\
          "from the first state of cooling"

    self._set_lock('dock')
    print "\n>>> Initial docking"
    dock_start_time = time.time()

    # Store docking parameters
    self._write_pkl_gz(join(self.dir['dock'],'dock_store.pkl.gz'),
      (convert_dictionary_relpath(
        dict([tp for tp in self._FNs.items() if not tp[0] in ['namd','vmd']]),
        relpath_o=None, relpath_n=self.dir['dock']),
       dict([tp for tp in self.params['dock'].items() \
             if not tp[0] in ['repX_cycles']])))

    # Set up the force field with full interaction grids
    self._set_universe_force_field(self.lambda_scalables)
  
    lambda_o = {'T':T_HIGH, 'MM':True, 'site':True, \
                'sLJr_crossed':False, 'crossed':False}
    for scalable in self._scalables:
      lambda_o[scalable] = 0
    self.dock_protocol = [lambda_o]

    state_start_time = time.time()
    self._random_translations_and_rotations()

    # Get interaction energies.
    # Loop over configurations, random rotations, and random translations
    E = {}
    for term in (['MM','site']+self._scalables):
      # MemoryError sometimes happens here
      E[term] = np.zeros((self.params['dock']['seeds_per_state'], \
                          self._max_n_rot,self._max_n_trans))

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
        # MemoryError sometimes happens here
        E_c[term] = np.ravel(E[term][:,:self._n_rot,:n_trans_n])
      (u_kln,N_k) = self._u_kln([E_c],\
        [lambda_o,self._next_dock_state(E=E_c, lambda_o=lambda_o)])
      du = u_kln[0,1,:] - u_kln[0,0,:]
      bootstrap_reps = 50
      f_grid0 = np.zeros(bootstrap_reps)
      for b in range(bootstrap_reps):
        f_grid0[b] = -np.log(np.mean(np.exp(
          -du[np.random.randint(0, len(du), len(du))])))
      f_grid0_std = f_grid0.std()
      converged = f_grid0_std<0.01
      if not converged:
        print "  with %s translations "%n_trans_n + \
              "the predicted free energy difference is %f (%f)"%(f_grid0.mean(),f_grid0_std)
        if n_trans_n == self._max_n_trans:
          break
        n_trans_o = n_trans_n
        n_trans_n = min(n_trans_n + 25, self._max_n_trans)

    if self._n_trans != n_trans_n:
      # Save number translations used
      self._n_trans = n_trans_n
      self._write_pkl_gz(join(self.dir['dock'],'random_trans.pkl.gz'),
          (self._n_trans, self._max_n_trans, self._random_trans))

    print "  %d ligand configurations were randomly docked "%len(cool0_Es_MM) + \
          "into the binding site "+ \
          "using %d translations and %d rotations "%(n_trans_n,self._n_rot) + \
          "in %f s"%(time.time()-state_start_time)
    print "  the predicted free energy difference between the" + \
          " first and second docking stages is %f (%f)"%(f_grid0.mean(),f_grid0_std)

    for term in E.keys():
      E[term] = np.ravel(E[term][:,:self._n_rot,:self._n_trans])

    if self.params['dock']['do_calc_random_dock_stats']:
        self.calc_random_dock_stats(E)

    # Main loop for initial docking: choose new thermodynamic variables, randomly select seeds, simulate
    while not self.dock_protocol[-1]['crossed']:
      # Determine next value of the protocol
      lambda_n = self._next_dock_state(E = E, lambda_o = lambda_o)
      self.dock_protocol.append(lambda_n)
      
      # Use free energy perturbation to estimate the free energy difference between
      # the first and second docking stage
      if len(self.dock_protocol)==2:
        (u_kln,N_k) = self._u_kln([E], [lambda_o, lambda_n])
        self.dock_protocol[0]['f_grid0'] = \
          -np.log(np.mean(np.exp(-(u_kln[0,1,:] - u_kln[0,0,:]))))
      
      # Randomly select seeds for new trajectory
      u_o = self._u_kln([E],[lambda_o])
      u_n = self._u_kln([E],[lambda_n])
      seedIndicies = choice(self.params['dock']['seeds_per_state'], np.exp(u_o-u_n))

      if len(self.dock_protocol)==2: # Cooling state 0 configurations, randomly oriented
        # Use the lowest energy configuration in the first docking stage for replica exchange
        ind = np.argmin(u_n)
        (c,i_rot,i_trans) = np.unravel_index(ind, (self.params['dock']['seeds_per_state'], self._n_rot, self._n_trans))
        repX_conf = np.add(np.dot(cool0_confs[c], self._random_rotT[i_rot,:,:]),\
                           self._random_trans[i_trans].array)
        self.confs['dock_replicas'] = [repX_conf]
        seeds = []
        for ind in seedIndicies:
          (c,i_rot,i_trans) = np.unravel_index(ind, (self.params['dock']['seeds_per_state'], self._n_rot, self._n_trans))
          seeds.append(np.add(np.dot(cool0_confs[c], self._random_rotT[i_rot,:,:]), self._random_trans[i_trans].array))
      else: # Seeds from last state
        seeds = [confs[ind] for ind in seedIndicies]

      # Simulate
      state_start_time = time.time()
      self._set_universe_force_field(lambda_n)
      (confs, Es_tot, g) = self._initial_sim_state('dock', seeds, normalize=False)
      
      # Store data
      self.confs['dock_replicas'].append(confs[np.argmin(Es_tot)])
      self.dock_protocol[-1]['delta_t'] = self.delta_t
      print '  generated %d configurations in %f s using with the following parameters:'%(\
        len(confs), time.time()-state_start_time)
      print self.dock_protocol[-1]

      # Get state energies
      self._set_universe_force_field(self.lambda_full)
      E = {}
      for term in (['MM','site','misc'] + self._scalables):
        E[term] = np.zeros(len(confs))
      for c in range(len(confs)):
        self.universe.setConfiguration(Configuration(self.universe,confs[c]))
        eT = self.universe.energyTerms()
        for (key,value) in eT.iteritems():
          E[term_map[key]][c] += value

      lambda_o = lambda_n

    self.dock_Es = [[] for s in range(len(self.dock_protocol))]
    print "  there are %d stages in the docking process"%(len(self.dock_protocol))

    # Save data
    self._write_crd(join(self.dir['dock'], 'dock_replicas.crd.gz'), self.confs['dock_replicas'], title='dock_replicas')
    self._write_pkl_gz(join(self.dir['dock'],'dock_protocol.pkl.gz'),self.dock_protocol)

    self._clear_lock('dock')
    print "\nElapsed time for initial docking: %f s"%(time.time()-dock_start_time)

  def dock(self):
    """
    Docks the ligand into the binding site
    by simulating at thermodynamic states
    between decoupled and fully interacting and
    between T_HIGH and T_TARGET
    """
    self._sim_process('dock')

  def calc_f_RL(self, readOnly=False):
    """
    Calculates the binding potential of mean force
    """
    if self.dock_protocol==[]:
      return # Initial docking is incomplete

    # Initialize variables by loading data or as empty lists
    f_RL_FN = join(self.dir['dock'],'f_RL.pkl.gz')
    if exists(f_RL_FN):
      dat = self._load_pkl_gz(f_RL_FN)
      if dat is not None:
        (self.f_grid, self.B, self.dock_equilibrated_cycle, self.dock_mean_acc) = dat
      else:
        for var in ['f_grid', 'B', 'dock_equilibrated_cycle', 'dock_mean_acc']:
          setattr(self,var,[])
    else:
      for var in ['f_grid', 'B', 'dock_equilibrated_cycle', 'dock_mean_acc']:
        setattr(self,var,[])
    if readOnly:
      return
  
    self._postprocess()
    self.calc_Psi()

    # Estimate cycle at which simulation has equilibrated
    K = len(self.dock_protocol)
    u_Ks = [self._u_kln([self.dock_Es[K-1][c]],[self.dock_protocol[K-1]]) for c in range(self._dock_cycle)]
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
      
      f_grid = [0, self.dock_protocol[0]['f_grid0']]

      # Extract energies between fromCycle and toCycle
      dock_Es = []
      for dock_Es_state in self.dock_Es:
        dock_Es.append(dock_Es_state[fromCycle:toCycle])
      
      # Use MBAR for the remaining states.
      (u_kln,N_k) = self._u_kln(\
        dock_Es[1:K],
        self.dock_protocol[1:K])
      f_grid += list(self._run_MBAR(u_kln,N_k)[1:] + f_grid[-1])

      self.f_grid.append(f_grid)

      # Binding PMF
      B = {}
      u_K = self._u_kln([dock_Es[K-1]],[self.dock_protocol[K-1]])

      for phase in ['Gas','GBSA']:
        du = \
          np.concatenate([self.dock_Es[K-1][c]['RL'+phase][:,-1] \
                            for c in range(fromCycle,toCycle)])/RT_TARGET - \
          (u_K + self.params['dock']['receptor_'+phase]/RT_TARGET)
        min_du = min(du)
        B_RL = -np.log(np.exp(-du+min_du).mean()) + min_du
        weight = np.exp(-du+min_du)
        weight = weight/sum(weight)
        Psi = np.concatenate([self.Psi[c][phase] for c in range(fromCycle,toCycle)])
        min_Psi = min(Psi)
        
        B[phase+'_MBAR_based'] = -self.f_cool[-1][-1] + f_grid[-1] + B_RL
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

      print "  calculated binding PMF of %f kT using cycles %d to %d"%(B['GBSA_MBAR_based'], fromCycle, toCycle-1)
      updated = True

    if updated:
      self._write_pkl_gz(f_RL_FN, (self.f_grid, self.B, self.dock_equilibrated_cycle, self.dock_mean_acc))

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

    self._random_translations_and_rotations()
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

  def calc_RMSD(self, rmsd_ref=None):
    """
    Uses VMD to calculate the root mean square deviation between 
    the docked conformations and the reference structure in ref_FN.
    """    
    if ((self.dock_protocol==[]) or (not self.dock_protocol[-1]['crossed'])):
      raise Exception('Initial docking not done!')
    if (self._FNs['vmd'] is None) or (not exists(self._FNs['vmd'])):
      print 'VMD not found! Cannot calculate RMSD.'
      return
    if rmsd_ref is None:
      if self.params['dock']['rmsd_ref'] is None:
        print 'Reference structure is needed for an RMSD calculation'
        return
      else:
        rmsd_ref = self.params['dock']['rmsd_ref']
    
    print "\n>>> RMSD calculation"
    
    rmsd_ref['vmdcrd'] = a.findPath([rmsd_ref['vmdcrd'],
      abspath(join(self.dir['dock'],rmsd_ref['vmdcrd']))])
    if rmsd_ref['vmdcrd'] is None:
      raise Exception('Reference structure is missing')

    rmsd_FN = join(self.dir['dock'],'rmsd.pkl.gz')
    if exists(rmsd_FN):
      self.rmsd = self._load_pkl_gz(rmsd_FN)
    else:
      self.rmsd = []
            
    ref_FN = rmsd_ref['vmdcrd']
    state = len(self.dock_protocol)-1

    # Write the DCD files
    newDCD = False
    for c in range(self._dock_cycle):
      while len(self.rmsd)<=c:
        self.rmsd.append({})
      if not ref_FN in self.rmsd[c].keys():
        crd = self._read_crd(join(self.dir['dock'],'dock%d_%d.crd.gz'%(state,c)), split=True)
        dcdFN = join(self.dir['dock'],'dock%d_%d.dcd'%(state,c))
        if not exists(dcdFN):
          self._write_dcd(dcdFN, crd)
        newDCD = True

    if newDCD:
      # Write the VMD script
      import tempfile
      scriptF = tempfile.NamedTemporaryFile(mode='w', dir=self.dir['dock'], delete=False)
      scriptF.write('''
set id_ref [mol new '''+self._FNs['prmtop'][rmsd_ref['species']]+''']
set sel_ref [atomselect $id_ref "'''+rmsd_ref['selection']+''' and not mass 1.008"]
mol addfile '''+abspath(ref_FN)+''' type crd molid $id_ref

set id_lig [mol new '''+self._FNs['prmtop']['L']+''']
set sel_lig [atomselect $id_lig "'''+rmsd_ref['selection']+''' and not mass 1.008"]
''')
      for c in range(self._dock_cycle):
        if not ref_FN in self.rmsd[c].keys():
          dcdFN = join(self.dir['dock'],'dock%d_%d.dcd'%(state,c))
          scriptF.write('''
mol addfile '''+dcdFN+''' type dcd molid $id_lig waitfor all
set fo [open "'''+join(self.dir['dock'],'rmsd%d_%d.txt'%(state,c))+'''" "w"]
set nfr [molinfo $id_lig get numframes]
for {set n 0} {$n < $nfr} {incr n} {
  $sel_lig frame $n
  set rmsd [measure rmsd $sel_ref $sel_lig]
  puts $fo "$rmsd"
}
close $fo
animate delete all
file delete '''+dcdFN+'\n')
      scriptF.write('\nquit\n')
      scriptF.close()

      # Execute VMD
      import subprocess
      vmd_out = open(join(self.dir['dock'], 'vmd.out'), 'w')
      subprocess.call([self._FNs['vmd'], '-nt', '-dispdev', 'text', '-e', scriptF.name],\
                      stdout=vmd_out)
      vmd_out.close()
      os.remove(scriptF.name)

      # Read the rmsd data
      for c in range(self._dock_cycle):
        if not ref_FN in self.rmsd[c].keys():
          temp_F = open(join(self.dir['dock'],'rmsd%d_%d.txt'%(state,c)),'r')
          self.rmsd[c][ref_FN] = np.array(
            [float(item) for item in temp_F.read().strip().split('\n')])
          temp_F.close()
          os.remove(join(self.dir['dock'],'rmsd%d_%d.txt'%(state,c)))

    self._write_pkl_gz(rmsd_FN, self.rmsd)

  def calc_Psi(self):
    """
    Calculate the interaction energy between the ligand and receptor. Also stores the ligand internal energy.
    """
    if ((self.dock_protocol==[]) or (not self.dock_protocol[-1]['crossed'])):
      raise Exception('Initial docking not done!')

    Psi_FN = join(self.dir['dock'],'Psi.pkl.gz')
    if exists(Psi_FN):
      self.Psi = self._load_pkl_gz(Psi_FN)
    else:
      self.Psi = []

    newPsi = False
    for c in range(self._dock_cycle):
      while len(self.Psi) <= c:
        self.Psi.append({})
      if not 'MM' in self.Psi[c].keys():
        self.Psi[c]['MM'] = self.dock_Es[-1][c]['MM']/RT_TARGET
        newPsi = True
      if not 'grid' in self.Psi[c].keys():
        self.Psi[c]['grid'] = \
          (self.dock_Es[-1][c]['LJr'] + \
          self.dock_Es[-1][c]['LJa'] + \
          self.dock_Es[-1][c]['ELE'])/RT_TARGET
        newPsi = True
      for phase in ['Gas','GBSA']:
        if not phase in self.Psi[c].keys():
          self.Psi[c][phase] = \
            (self.dock_Es[-1][c]['RL'+phase][:,-1] - \
            self.dock_Es[-1][c]['L'+phase][:,-1] - \
            self.params['dock']['receptor_'+phase])/RT_TARGET
          newPsi = True

    if newPsi:
      self._write_pkl_gz(Psi_FN, self.Psi)

  def plot_energies(self, process='cool', firstCycle=0, toCycle=None):
    """
    Plots timeseries and histograms of the energies for each state.
    Requires matplotlib extension for python.
    """
    try:
      import matplotlib.pyplot as plt  # @UnresolvedImport
    except:
      print 'plot_energies requires matplotlib'
      return
    
    K = len(getattr(self,process+'_protocol'))
    
    if toCycle is None:
      toCycle = getattr(self,'_%s_cycle'%process)
    Es = [getattr(self,process+'_Es')[0][firstCycle:toCycle]]
    for Es_state in getattr(self,process+'_Es'):
      Es.append(Es_state[firstCycle:toCycle])
    (u_kln,N_k) = self._u_kln(Es, getattr(self,process+'_protocol'))
    
    plt.figure(0)
    for k in range(K-1):
        plt.plot(u_kln[k,k,:N_k[k]],'.-')
    
    plt.figure(1)
    for k in range(K-1):
        plt.hist(u_kln[k,k,:N_k[k]])

  def view_samples(self, process='dock', state=None,
                   includeReceptor=False, displayComplex=False,
                   overlap=False, showGrid=True, conf_indices=None,
                   saveImage=False, execute=True, colorID=None, clearDCD=True):
    """
    Views the trajectory with VMD
    
    state - the simulation state to view
    includeReceptor - show the receptor as well as the ligand
    overlap - show ligands overlapping each other, rather than a trajectory
    showGrid - shows the Lennard-Jones repulsive grid
    stride - steps between snapshots
    """
    if (self._FNs['vmd'] is None) or (not exists(self._FNs['vmd'])):
      print 'VMD not found! Cannot view trajectory'
      return
    if includeReceptor:
      if (self._FNs['prmtop']['RL'] is None):
        raise Exception('AMBER prmtop file for complex required!')
      prmtopFN = self ._FNs['prmtop']['RL']
    else:
      if (self._FNs['prmtop']['L'] is None):
        raise Exception('AMBER prmtop file for the ligand required!')
      prmtopFN = self._FNs['prmtop']['L']
    
    # Default is last state
    if state is None:
      state = len(getattr(self,process+'_protocol'))-1

    if process=='dock':
      firstCycle = self.dock_equilibrated_cycle[-1]
    elif process=='cool':
      firstCycle = self.cool_equilibrated_cycle[-1]
    else:
      raise Exception('Invalid Process')
    
    # Convert from crd to dcd files
    for c in range(firstCycle,getattr(self,'_%s_cycle'%process)):
      if not exists(join(self.dir[process],process+'%d_%d.dcd'%(state,c))):
        crd = self._read_crd(join(self.dir[process],process+'%d_%d.crd.gz'%(state,c)), split=True)
        self._write_dcd(join(self.dir[process],process+'%d_%d.dcd'%(state,c)), crd)

    script = '# View samples\n'
    script += 'update off\n'

    script += 'set ID [mol new {' + prmtopFN + '}]\n'

    if colorID is not None:
      script += 'mol modcolor 0 $ID ColorID %d\n'%colorID

    for c in range(firstCycle,getattr(self,'_%s_cycle'%process)):
      script += 'mol addfile {' + join(self.dir[process],process+'%d_%d.dcd'%(state,c)) + '} waitfor all molid $ID\n'

    if overlap:
      script += 'set nfr [molinfo $ID get numframes]\n'
      script += 'mol drawframes $ID 0 0:[expr $nfr-1]\n'
    else:
      script += 'mol modstyle 0 0 Licorice 0.300000 10.000000 10.000000\n'
    script += 'display projection Orthographic\n'
  
    if showGrid:
      if self._FNs['grids']['LJr'][-3:]=='.gz':
        import shutil
        shutil.copyfile(self._FNs['grids']['LJr'], 'temp_grid.LJr.dx.gz')
        os.system('gunzip temp_grid.LJr.dx.gz')
        grid_LJr = 'temp_grid.LJr.dx'
      else:
        grid_LJr = self._FNs['grids']['LJr']
      
      script += 'mol addfile {'+grid_LJr+'} type dx first 0 last -1 step 1 waitfor 1 volsets {0 } 0\n'
      # Lennard-Jones Repulsive
      script += 'mol representation Isosurface 75.0 0 2 1 1 1\n'
      script += 'mol addrep 0\n'

    if displayComplex:
      script += 'set complexID [mol new {%s}]\n'%self._FNs['prmtop']['RL']
      script += 'mol addfile {%s} type crd\n'%(self.params['dock']['rmsd_ref']['vmdcrd'])
      script += '''
mol modstyle 0 $complexID NewCartoon 0.300000 10.000000 4.100000 0
mol color ColorID 0
mol representation Licorice 0.300000 10.000000 10.000000
mol selection resname iga
mol material Opaque
mol addrep $complexID
mol color Name
mol representation Lines 2.0
mol selection same residue as within 10 of resname iga
mol material Opaque
mol addrep $complexID
'''

    script += 'display resize 600 600\n'
  
    script += 'axes location off\n'
    script += 'rotate y by -90.000000\n'
    script += 'scale by 5\n'
    script += 'translate by 0 -0.5 0\n'
    script += 'display nearclip set 0.200000\n'

    script += 'draw color green\n'
    script += 'draw text {20 6 32} "Grid scaling: %.2e" size 2 thickness 2\n'%self.dock_protocol[state]['LJr']
    script += 'draw text {20 4 32} "T: %.2f K" size 2 thickness 2\n'%self.dock_protocol[state]['T']

    script += 'update\n'
    script += 'update on\n'

    if clearDCD:
      for c in range(getattr(self,'_%s_cycle'%process)):
        script += 'file delete {' + join(self.dir[process],process+'%d_%d.dcd'%(state,c)) + '}\n'
    if showGrid and self._FNs['grids']['LJr'][-3:]=='.gz':
      script += 'file delete {' + grid_LJr + '}\n'
    if saveImage:
      script += 'render TachyonInternal ' + os.path.join(self.dir['dock'],'state%d.tga'%state) + ' /usr/bin/open %s\n'
      script += 'quit\n'

    if execute:
      import tempfile
      scriptFN = tempfile.mktemp()
      scriptF = open(scriptFN,'w')
      scriptF.write(script)
      scriptF.close()
      import subprocess
      subprocess.Popen([self._FNs['vmd'], '-nt', '-e', scriptFN])

    return script
    
  ######################
  # Internal Functions #
  ######################

  def _set_universe_force_field(self,lambda_n):
    """
    Sets the force field parameters of the universe to values appropriate for the given lambda_n dictionary.
    The elements in the dictionary lambda_n can be:
      MM - True, to turn on Generalized AMBER force field
      site - True, to turn on the binding site
      sLJr - scaling of the soft Lennard-Jones repulsive grid
      sLJa - scaling of the soft Lennard-Jones attractive grid
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
          grid_FN = self._FNs['grids'][{'sLJr':'LJr','sLJa':'LJa','LJr':'LJr','LJa':'LJa','ELE':'ELE'}[scalable]]
          grid_scaling_factor = 'scaling_factor_' + \
            {'sLJr':'LJr','sLJa':'LJa','LJr':'LJr','LJa':'LJa','ELE':'electrostatic'}[scalable]
          if scalable=='LJr':
            from ForceFields.TrilinearISqrtGrid.TrilinearISqrtGrid import TrilinearISqrtGridForceField
            self._forceFields[scalable] = TrilinearISqrtGridForceField(grid_FN, lambda_n[scalable], grid_scaling_factor, grid_name=scalable, max_val=-1)
          else:
            from ForceFields.TrilinearGrid.TrilinearGrid import TrilinearGridForceField
            self._forceFields[scalable] = TrilinearGridForceField(
              grid_FN, lambda_n[scalable], grid_scaling_factor, grid_name=scalable,
              max_val=10. if scalable[0]=='s' else 10000.)
          print '  %s loaded in %f s'%(grid_FN, time.time()-loading_start_time)

        # Set the force field strength to the desired value
        self._forceFields[scalable].strength = lambda_n[scalable]
        fflist.append(self._forceFields[scalable])

    compoundFF = fflist[0]
    for ff in fflist[1:]:
      compoundFF += ff
    self.universe.setForceField(compoundFF)

  def _random_translations_and_rotations(self):
    """
    Either loads or generates the random translations and rotations for the first state of docking.
    """
    if not hasattr(self,'_random_trans'):
      random_translation_FN = join(self.dir['dock'],'random_trans.pkl.gz')
      if exists(random_translation_FN):
        (self._n_trans, self._max_n_trans, self._random_trans) = \
          self._load_pkl_gz(random_translation_FN)
      else:
        self._max_n_trans = 10000
        # Default density of points is 500 per nm**3
        self._n_trans = min(np.int(np.ceil(self._forceFields['site'].volume*self.params['dock']['site_density'])),self._max_n_trans)
        self._random_trans = np.ndarray((self._max_n_trans), dtype=Vector)
        for ind in range(self._max_n_trans):
          self._random_trans[ind] = Vector(self._forceFields['site'].randomPoint())
        self._write_pkl_gz(random_translation_FN,
          (self._n_trans, self._max_n_trans, self._random_trans))
    else:
      self._max_n_trans = self._random_trans.shape[0]

    if not hasattr(self,'_random_rotT'):
      random_rotation_FN = join(self.dir['dock'],'random_rotT.pkl.gz')
      if exists(random_rotation_FN):
        (self._n_rot, self._max_n_rot, self._random_rotT) = \
          self._load_pkl_gz(random_rotation_FN)
      else:
        self._max_n_rot = 100
        self._n_rot = 100
        self._random_rotT = np.ndarray((self._max_n_rot,3,3))
        for ind in range(self._max_n_rot):
          self._random_rotT[ind,:,:] = np.transpose(random_rotate())
        self._write_pkl_gz(random_rotation_FN,
          (self._n_rot, self._max_n_rot, self._random_rotT))
    else:
      self._n_rot = self._random_rotT.shape[0]

  def _initial_sim_state(self, process, seeds, normalize=False):
    """
    Initializes a state, returning the configurations and potential energy.
    """
    
    all_confs = []
    all_potEs = []
    Hs = []

    for seed in seeds:
      self.universe.setConfiguration(Configuration(self.universe,seed))
      (confs_s, potEs_s, Hs_s, self.delta_t) = self.sampler['init'](
        steps = self.params[process]['steps_per_seed'],
        T=self.T, delta_t=self.delta_t,
        steps_per_trial = self.params[process]['steps_per_seed']/5,
        normalize=normalize, adapt=True)
      # Restrict the range of the time step
      self.delta_t = min(max(self.delta_t, 0.3*MMTK.Units.fs), 3.*MMTK.Units.fs)
      all_confs.append(confs_s)
      all_potEs.append(np.array(potEs_s))
      Hs.append(Hs_s)
    Hs = np.mean(np.array(Hs))

    import sys
    original_stderr = sys.stderr
    sys.stderr = NullDevice()

    # Decorrelation
    import pymbar.timeseries
    if len(seeds)>1:
      g = pymbar.timeseries.statisticalInefficiencyMultiple(all_potEs)
    else:
      g = pymbar.timeseries.statisticalInefficiency(all_potEs)
    stride = int(np.ceil(g))

    sys.stderr = original_stderr

    confs = []
    potEs = []
    for s in range(len(seeds)):
      nsnaps = len(all_confs[s])
      uncorr_indicies = range(min(stride-1, nsnaps-1), nsnaps, stride)
      confs.extend([all_confs[s][i] for i in uncorr_indicies])
      potEs.extend([all_potEs[s][i] for i in uncorr_indicies])

    return (confs, np.array(potEs), g)
    
  def _replica_exchange(self, process):
    """
    Performs a cycle of replica exchange
    """
    if not process in ['dock','cool']:
      raise Exception('Process must be dock or cool')

    cycle_start_time = time.time()

    self._set_lock(process)

    if process=='cool':
      terms = ['MM']
    else:
      terms = ['MM','site','misc'] + self._scalables

    if not process+'_replicas' in self.confs.keys():
      self.confs[process+'_replicas'] = self._read_crd(join(self.dir[process],process+'_replicas.crd.gz'))
    confs = self.confs[process+'_replicas']
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
    
    storage = {}
    for var in ['confs','state_inds','energies']:
      storage[var] = []
    
    state_inds = range(K)
    E = {}
    for term in terms:
      E[term] = np.zeros(K, dtype=float)
    Ht = np.zeros(K, dtype=int)
    for sweep in range(self.params[process]['sweeps_per_cycle']):
      # Sample within each state
      for k in range(K):
        self._set_universe_force_field(lambdas[state_inds[k]])
        self.universe.setConfiguration(Configuration(self.universe, confs[k]))
        # By far, this step takes the most time within the cycle
        (confs_k, potEs_k, Hs_k, delta_t) = self.sampler[process](\
          steps = self.params[process]['steps_per_sweep'],\
          T=self.T, delta_t=self.delta_t,\
          normalize=(process=='cool'), adapt=False)
        confs[k] = np.copy(confs_k[-1])
        if process=='cool':
          E['MM'][k] = potEs_k[-1]
        Ht[k] += Hs_k
      if process=='dock':
        # Get energies
        self._set_universe_force_field(self.lambda_full)
        for term in terms:
          E[term] = np.zeros(K, dtype=float)
        for k in range(K):
          self.universe.setConfiguration(Configuration(self.universe, confs[k]))
          eT = self.universe.energyTerms()
          for (key,value) in eT.iteritems():
            E[term_map[key]][k] += value
      # Calculate u_ij (i is the replica, and j is the configuration),
      #    a list of arrays
      (u_ij,N_k) = self._u_kln(E, [lambdas[state_inds[c]] for c in range(K)])
      # Do the replica exchange
      for attempt in range(100):
        for (a,b) in pairs_to_swap:
          ddu = -u_ij[a][b]-u_ij[b][a]+u_ij[a][a]+u_ij[b][b]
          if (ddu>0) or (np.random.uniform()<np.exp(ddu)):
            u_ij[a],u_ij[b] = u_ij[b],u_ij[a]
            state_inds[a],state_inds[b] = state_inds[b],state_inds[a]
      # Store data in local variables
      storage['confs'].append(list(confs))
      storage['state_inds'].append(list(state_inds))
      storage['energies'].append(copy.deepcopy(E))

    print "  generated configurations for %d replicas"%(len(confs)) + \
      " in replica exchange cycle %d"%getattr(self,'_%s_cycle'%process) + \
      " in %f s"%(time.time()-cycle_start_time)

    Ht = Ht/self.params[process]['sweeps_per_cycle']
    storage['Ht'] = Ht

    tau2 = 1.0
    stride = 1
    store_indicies = np.array(\
      range(min(stride-1,self.params[process]['sweeps_per_cycle']-1), self.params[process]['sweeps_per_cycle'], stride), dtype=int)
    nsaved = len(store_indicies)

    # Get indicies for storing global variables
    inv_state_inds = np.zeros((nsaved,K),dtype=int)
    for snap in range(nsaved):
      state_inds = storage['state_inds'][store_indicies[snap]]
      for state in range(K):
        inv_state_inds[snap][state_inds[state]] = state

    # Store global variables and data
    for state in range(K):
      E_state = {}
      if state==0:
        E_state['repXpath'] = inv_state_inds
      for term in terms:
        E_state[term] = np.array([storage['energies'][store_indicies[snap]][term][inv_state_inds[snap][state]] for snap in range(nsaved)])
      getattr(self,process+'_Es')[state].append(E_state)

    self.confs[process+'_replicas'] = [storage['confs'][store_indicies[-1]][inv_state_inds[-1][state]] for state in range(K)]
    self._write_crd(join(self.dir[process],process+'_replicas.crd.gz'), self.confs[process+'_replicas'], title=(process+'_replicas'))

    for state in range(K):
      if self.params[process]['keep_intermediate'] or ((process=='cool') and (state==0)) or (state==(K-1)):
        confs = [storage['confs'][store_indicies[snap]][inv_state_inds[snap][state]] for snap in range(nsaved)]
        self._write_crd(
          join(self.dir[process],process+'%d_%d.crd.gz'%(state,getattr(self,'_%s_cycle'%process))),
          confs,
          title='State_%d_%d'%(state,getattr(self,'_%s_cycle'%process)))

    postprocess_start_time = time.time()
    setattr(self,'_%s_cycle'%process,getattr(self,'_%s_cycle'%process) + 1)
    if process=='cool':
      self._postprocess([('cool',-1,-1,'L')])
      self.calc_f_L()
    elif process=='dock':
      self._postprocess()
      self.calc_f_RL()

    self._clear_lock(process)

    print '  postprocessed data in %f s\n'%(time.time()-postprocess_start_time)

  def _sim_process(self, process):
    """
    Simulate and analyze a cooling or docking process.
    
    As necessary, first conduct an initial docking
    and then run a desired number of replica exchange cycles.
    """
    if getattr(self,process+'_protocol')==[]:
      getattr(self,'initial_'+process)()
    if (self.params[process]['repX_cycles'] is not None) and \
       (getattr(self,'_%s_cycle'%process) < self.params[process]['repX_cycles']):
      print "\n>>> Replica exchange sampling for the {0}ing process".format(process)
      import time
      repEx_start_time = time.time()
      start_cycle = getattr(self,'_%s_cycle'%process)
      while getattr(self,'_%s_cycle'%process) < self.params[process]['repX_cycles']:
        self._replica_exchange(process)
      print "Elapsed time for %d cycles of replica exchange was %f s"%(\
         (getattr(self,'_%s_cycle'%process)-start_cycle), \
         (time.time()-repEx_start_time))

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
        f_k_BAR[k+1] = pymbar.computeBAR(w_F, w_R, relative_tolerance=0.000001, verbose=False, compute_uncertainty=False)
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

  def _next_dock_state(self, E=None, lambda_o=None):
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
      if self.params['dock']['legs']==1:
        # Change grid scaling and temperature simultaneously
        a = lambda_o['LJr']
        T = lambda_o['T']
        Psi_g = self._u_kln([E],
          [{'LJr':1,'LJa':1,'ELE':1}], noBeta=True)
        U_RL_g = self._u_kln([E],
          [{'MM':True,'site':True,'T':lambda_o['T'],\
          'LJr':a,'LJa':a,'ELE':a}], noBeta=True)
        sigma_dl_a = Psi_g.std()/(R*T) + \
          U_RL_g.std()*np.abs(T_TARGET-T_HIGH)/(R*T*T)
        if sigma_dl_a>0:
          dL = self.params['dock']['therm_speed']/sigma_dl_a
          a = min(lambda_o['LJr'] + dL, 1)
        else:
          a = 1
        lambda_n['LJr'] = a
        lambda_n['LJa'] = a
        lambda_n['ELE'] = a
        lambda_n['T'] = a*(T_TARGET-T_HIGH) + T_HIGH
        if a == 1:
          lambda_n['crossed'] = True
      else:
        raise Exception("Multiple legs not supported")

    return lambda_n

  def _postprocess(self,
      conditions=[('original',0, 0,'R'), ('cool',-1,-1,'L'), \
                  ('dock',   -1,-1,'L'), ('dock',-1,-1,'RL')],
      phases=['GBSA','Gas'],
      readOnly=False, debug=False):
    """
    Obtains the NAMD energies of all the conditions using all the phases.  
    Saves both MMTK and NAMD energies after NAMD energies are estimated.
    """
    # NAMD ENERGY FIELDS:
    # 0. TS 1. BOND 2. ANGLE 3. DIHED 4. IMPRP 5. ELECT 6. VDW 7. BOUNDARY
    # 8. MISC 9. KINETIC 10. TOTAL 11. TEMP 12. POTENTIAL 13. TOTAL3 14. TEMPAVG
    # The saved fields are energyFields=[1, 2, 3, 4, 5, 6, 8, 12],
    # and thus the new indicies are
    # 0. BOND 1. ANGLE 2. DIHED 3. IMPRP 4. ELECT 5. VDW 6. MISC 7. POTENTIAL
    
    updated = []
    crd_FN_o = ''
    
    # state == -1 means the last state
    # cycle == -1 means all cycles
    for (p, state, cycle, moiety) in conditions:
      # Check that the values are legitimate
      if not p in ['cool','dock','original']:
        raise Exception("Type should be in ['cool', 'dock', 'original']")
      if not moiety in ['R','L', 'RL']:
        raise Exception("Species should in ['R','L', 'RL']")
    
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
          fixed = {'R':self._FNs['fixed_atoms']['R'],
                   'L':None,
                   'RL':self._FNs['fixed_atoms']['RL']}[moiety]

          label = moiety+phase
          crd_FN = join(p_dir,'%s.crd.gz'%(prefix))
          dcd_FN = join(p_dir,'%s.%s.dcd'%(prefix,moiety))
    
          # Continue if the energies are already in memory
          if readOnly or \
            (p == 'original' and (label in getattr(self,p+'_Es')[state][c].keys())) \
            or (('MM' in getattr(self,p+'_Es')[state][c].keys()) and \
                (label in getattr(self,p+'_Es')[state][c].keys()) and \
                (len(getattr(self,p+'_Es')[state][c]['MM'])==\
                 len(getattr(self,p+'_Es')[state][c][label]))):
            pass
          else:
            # Do the calculation
            # Obtain configurations
            if (moiety=='R'):
              confs = self._receptor_conf
            else:
              if crd_FN != crd_FN_o:
                confs = self._read_crd(crd_FN)
              crd_FN_o = crd_FN
            
            # Write the configurations
            if not exists(dcd_FN):
              self._write_dcd(dcd_FN, confs,
                includeReceptor=(moiety.find('R')>-1),
                includeLigand=(moiety.find('L')>-1))
            
            # Run NAMD
            import AlGDock.NAMD
            energyCalc = AlGDock.NAMD.NAMD(prmtop=self._FNs['prmtop'][moiety],
                                   inpcrd=self._FNs['inpcrd'][moiety],
                                   fixed=fixed,
                                   solvent={'GBSA':'GBSA','Gas':'vacuum'}[phase],
                                   useCutoff=(phase=='GBSA'),
                                   namd_command=self._FNs['namd'])
            E = energyCalc.energies_PE(\
              join(p_dir,'%s.%s%s'%(prefix,moiety,phase)), dcd_FN,
              energyFields=[1, 2, 3, 4, 5, 6, 8, 12],
              keepScript=debug, writeEnergyDatGZ=False)
            getattr(self,p+'_Es')[state][c][label] = \
              np.array(E, dtype=float)*MMTK.Units.kcal/MMTK.Units.mol
          
            # Updated
            if not (p, c) in updated:
              updated.append((p, c))

        # Clean up
        if (not debug) and exists(dcd_FN):
          os.remove(dcd_FN)

    # If the receptor energy is updated
    if ((self.params['dock']['receptor_Gas'] is None) or \
        (self.params['dock']['receptor_GBSA'] is None)) and \
       ((self.original_Es[0][0]['RGas'] is not None) or \
        (self.original_Es[0][0]['RGBSA'] is not None)):
      for phase in ['Gas','GBSA']:
        self.params['dock']['receptor_'+phase] = self.original_Es[0][0]['R'+phase]

      self._write_pkl_gz(join(self.dir['dock'],'dock_store.pkl.gz'),
        (convert_dictionary_relpath(
          dict([tp for tp in self._FNs.items() if not tp[0] in ['namd','vmd']]),
          relpath_o=None, relpath_n=self.dir['dock']),
         dict([tp for tp in self.params['dock'].items() if not tp[0] in ['repX_cycles']])))

    # Save progress
    for (p,c) in updated:
      print "  postprocessed {0}ing cycle {1} with NAMD".format(p,c)
      if p=='cool' or p=='dock':
        nstates = len(getattr(self,p+'_protocol'))
        energyFN = join(self.dir[p], '%s_Es_%d.pkl.gz'%(p, c))
        energies = [getattr(self, p+'_Es')[state][c] for state in range(nstates)]
        self._write_pkl_gz(energyFN, energies)

  def _load_cool0(self):
    # Load configurations
    confs = []
    c = 0
    while exists(join(self.dir['cool'],'cool0_%d.crd.gz'%c)):
      confs += self._read_crd(join(self.dir['cool'],'cool0_%d.crd.gz'%c))
      c += 1
    # Load energies
    E_MM = []
    for k in range(len(self.cool_Es[0])):
      E_MM += list(self.cool_Es[0][k]['MM'])
    # Check for consistency
    if len(E_MM)!=len(confs):
      raise Exception('Error in loading the first cooling state')
    # Indicies to be used in first docking state
    random_dock_inds = np.array(np.linspace(0,len(E_MM), \
      self.params['dock']['seeds_per_state'],endpoint=False),dtype=int)
    return ([confs[ind] for ind in random_dock_inds], \
            [E_MM[ind]  for ind in random_dock_inds])

  def _read_crd(self, crdFN, split=True):
    """ 
    Reads an AMBER or VMD format coordinate file.
    """
    if crdFN[-3:] == '.gz':
      crdF = gzip.open(crdFN,'r')
    else:
      crdF = open(crdFN,'r')
    dat = crdF.read().strip().split('\n')
    crdF.close()

    title = dat.pop(0) # Title

    if len(dat[0].split())>1:
      # VMD format (does not specify number of atoms)
      crd = []
      for line in dat:
        crd = crd + [float(x)/10.0 for x in line.split()]
      crd = np.resize(crd,(len(crd)/3,3))
      self._write_crd(crdFN, crd, title=title)
    else:
      # AMBER format
      natoms = int(dat.pop(0)) # Number of atoms
      if split and (natoms!=self.universe.numberOfAtoms()):
        print 'Incorrect number of atoms in crd file'
        return np.array([])

      w = 12
      crd = []
      for line in dat:
        crd = crd + [float(line[x:x+w])/10.0 for x in range(0,len(line),w)]
      crd = np.resize(crd,(len(crd)/3,3))

    if split:
      crd = np.vsplit(crd,crd.shape[0]/self.universe.numberOfAtoms())
      print "  read %d configurations from %s"%(len(crd), crdFN)
    return crd
    
  def _write_crd(self, crdFN, crd, title='', appendIfExists=False, vmdStyle=False):
    """
    Writes an AMBER format trajectory file
    """
    if (appendIfExists and exists(crdFN)):
      if crdFN[-3:] == '.gz':
        crdF = gzip.open(crdFN,'a')
      else:
        crdF = open(crdFN,'a')
    else:
      if exists(crdFN):
        os.rename(crdFN,crdFN+'.BAK')
      if crdFN[-3:] == '.gz':
        crdF = gzip.open(crdFN,'w')
      else:
        crdF = open(crdFN,'w')
      # Write the header
      crdF.write(title+'\n') # Title
      if not vmdStyle:
        crdF.write('%d\n'%self.universe.numberOfAtoms())

    for atom in (np.vstack(crd)/MMTK.Units.Ang):
      crdF.write('%12.7f%12.7f%12.7f\n'%(atom[0],atom[1],atom[2]))
    crdF.close()
    print "  wrote %d configurations to %s"%(len(crd), crdFN)

  def _write_dcd(self, dcd_file_name, trajectory_array,
      includeLigand=True, includeReceptor=False,
      factor=1.0/MMTK.Units.Ang,
      delta_t=0.1):
    """
    Writes a DCD file for a trajectory.
    If includeReceptor==True, the receptor coordinates are included.
    """
    import MMTK_DCD  # @UnresolvedImport
    from Scientific import N

    prmtop_atom_order = [ref.number for ref in self.molecule.prmtop_order]

    i_start = 0       # always start at frame 0
    n_savc  = 1       # save every frame
    
    n_atoms = 0
    if includeReceptor:
      receptor_x0 = factor*self._receptor_conf[:self._ligand_first_atom,0]
      receptor_y0 = factor*self._receptor_conf[:self._ligand_first_atom,1]
      receptor_z0 = factor*self._receptor_conf[:self._ligand_first_atom,2]
      receptor_x1 = factor*self._receptor_conf[self._ligand_first_atom:,0]
      receptor_y1 = factor*self._receptor_conf[self._ligand_first_atom:,1]
      receptor_z1 = factor*self._receptor_conf[self._ligand_first_atom:,2]
      n_atoms += self._receptor_conf.shape[0]
    if includeLigand:
      n_atoms += len(self.molecule.atoms)
    n_snaps = len(trajectory_array)

    fd = MMTK_DCD.writeOpenDCD(dcd_file_name, n_atoms, n_snaps,
                               i_start, n_savc, delta_t)

    if includeReceptor and includeLigand:
      for array in trajectory_array:
        array = factor*array
        x = N.concatenate((receptor_x0,N.take(array[:,0],prmtop_atom_order),receptor_x1)).astype(N.Float16)
        y = N.concatenate((receptor_y0,N.take(array[:,1],prmtop_atom_order),receptor_y1)).astype(N.Float16)
        z = N.concatenate((receptor_z0,N.take(array[:,2],prmtop_atom_order),receptor_z1)).astype(N.Float16)
        MMTK_DCD.writeDCDStep(fd, x, y, z)
      MMTK_DCD.writeCloseDCD(fd)
    elif includeLigand:
      for array in trajectory_array:
        array = factor*array
        x = N.take(array[:,0], prmtop_atom_order).astype(N.Float16)
        y = N.take(array[:,1], prmtop_atom_order).astype(N.Float16)
        z = N.take(array[:,2], prmtop_atom_order).astype(N.Float16)
        MMTK_DCD.writeDCDStep(fd, x, y, z)
      MMTK_DCD.writeCloseDCD(fd)
    else:
      x = N.concatenate((receptor_x0,receptor_x1)).astype(N.Float16)
      y = N.concatenate((receptor_y0,receptor_y1)).astype(N.Float16)
      z = N.concatenate((receptor_z0,receptor_z1)).astype(N.Float16)
      MMTK_DCD.writeDCDStep(fd, x, y, z)
      MMTK_DCD.writeCloseDCD(fd)

  def _load_pkl_gz(self, FN):
    if exists(FN):
      F = gzip.open(FN,'r')
      try:
        data = pickle.load(F)
      except:
        print 'Could not load '+FN
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
    print "  wrote to "+FN

  def _load_progress(self, p):
    setattr(self, p+'_protocol', 
            self._load_pkl_gz(join(self.dir[p], p+'_protocol.pkl.gz'))) 
    if getattr(self, p+'_protocol') is None:
      setattr(self, p+'_protocol', [])
      setattr(self, '_%s_cycle'%p, 0)
      return

    # Load energies in each cycle
    if (not hasattr(self, p+'_Es')) or (getattr(self, p+'_Es')==[]):
      nstates = len(getattr(self, p+'_protocol'))
      setattr(self, p+'_Es', [[] for state in range(nstates)])

      c = 0
      while exists(join(self.dir[p], '%s_Es_%d.pkl.gz'%(p, c))):
        cycle_Es = self._load_pkl_gz(join(self.dir[p], '%s_Es_%d.pkl.gz'%(p, c)))
        for state in range(len(cycle_Es)):
          getattr(self, p+'_Es')[state].append(cycle_Es[state])
        c += 1

    if (not getattr(self, p+'_protocol')==[]) and \
       (getattr(self, p+'_protocol')[-1]['crossed']):
      setattr(self,'_%s_cycle'%p, len(getattr(self, p+'_Es')[0]))
    else:
      setattr(self,'_%s_cycle'%p, 0)
     
  def _set_lock(self, p):
    lockFN = join(self.dir[p],'.lock')
    if exists(lockFN):
      raise Exception(p + ' is locked')
    else:
      lockF = open(lockFN,'w')
      lockF.close()

  def _clear_lock(self, p):
    lockFN = join(self.dir[p],'.lock')
    if exists(lockFN):
      os.remove(lockFN)

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
  parser.add_argument('--grid_force', 
    help='PDB file with grid force multipliers (for NAMD)')
  # Simulation settings and constants
  #   Run-dependent 
  parser.add_argument('--cool_repX_cycles', type=int,
    help='Number of replica exchange cycles for cooling')
  parser.add_argument('--dock_repX_cycles', type=int,
    help='Number of replica exchange cycles for docking')
  parser.add_argument('--run_type',
    choices=['initial_cool','cool','initial_dock','dock','all'],
    help='Type of calculation to run')
  #   Defaults
  parser.add_argument('--protocol', choices=['Adaptive','Set'],
    help='Approach to determining series of thermodynamic states')
  parser.add_argument('--therm_speed', type=float,
    help='Thermodynamic speed during adaptive simulation')
  parser.add_argument('--sampler',
    choices=['HMC','NUTS','VV'],
    help='Sampling method')
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
  parser.add_argument('--keep_intermediate', action='store_true',
    help='Keep configurations for intermediate states?')
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
    help='Keep configurations for intermediate states?')
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
  parser.add_argument('--legs', type=int, choices=[1],
    help='Number of legs in the docking thermodynamic cycle. 1: Scale LJ and ELE grids and lower the temperature simultaneously.')
  # For replica exchange
  parser.add_argument('--dock_sweeps_per_cycle', type=int,
    help='Number of replica exchange sweeps per cycle')
  parser.add_argument('--dock_steps_per_sweep', type=int,
    help='Number of MD steps per replica exchange sweep')
  parser.add_argument('--dock_keep_intermediate', action='store_true',
    help='Keep configurations for intermediate states?')
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
    help='Density of center-of-mass points in the first docking stage')
  # Additional calculations
  parser.add_argument('--do_calc_random_dock_stats', action='store_true',
    help='Calculate convergence statistics for random placement of the ligand into the active site')
  parser.add_argument('--rmsd_ref',
    help='Performs an RMSD calculation. : 1. the file name, 2. the species, and 3. a VMD selection string')
  parser.add_argument('--receptor_GBSA', type=float, nargs='+',
    help='Receptor potential energies in GBSA implicit solvent (in units of kJ/mol)')
  parser.add_argument('--receptor_Gas', type=float, nargs='+',
    help='Receptor potential energies in gas phase (in units of kJ/mol)')
  
  args = parser.parse_args()
  self = Docker(**vars(args))
