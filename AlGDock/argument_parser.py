import os
import multiprocessing

import cPickle as pickle
import gzip

import numpy as np

from collections import OrderedDict

from AlGDock import arguments
from AlGDock import dictionary_tools
from AlGDock import path_tools

from AlGDock.IO import load_pkl_gz
from AlGDock.IO import write_pkl_gz


class SimulationArguments:
  """Stores simulation arguments

  ...

  Attributes
  ----------
  cores : int
    The number of CPU cores to use
  dir : dict of str
    Directories where parameters and data are stored
  do_CD : bool
    Whether to perform simulation for states CD
  params : OrderedDict of OrderedDict
    Simulation parameters for BC and CD
  FNs : OrderedDict of OrderedDict
    Paths of input files
  original_Es : list of list of Dictionaries
    Energies of the receptor snapshot in different force fields
  random_seed : int
    The random number seed
  toClear : list
    Files to remove at the end of the calculation
  """
  def __init__(self, **kwargs):
    """Parses the input arguments

    """

    # Set undefined keywords to None
    for key in arguments.args.keys():
      if not key in kwargs.keys():
        kwargs[key] = None
    if kwargs['dir_grid'] is None:
      kwargs['dir_grid'] = ''

    # Multiprocessing options.
    # Default is to use 1 core.
    # If cores is a number, then that number (or the maximum number)
    # of cores will be used.

    # Default
    availablecores = multiprocessing.cpu_count()
    if kwargs['cores'] is None:
      self.cores = 1
    elif (kwargs['cores'] == -1):
      self.cores = availablecores
    else:
      self.cores = min(kwargs['cores'], availablecores)
    print "using %d/%d available cores" % (self.cores, availablecores)

    if kwargs['rotate_matrix'] is not None:
      self._view_args_rotate_matrix = kwargs['rotate_matrix']

    if kwargs['random_seed'] is None:
      self.random_seed = 0
    else:
      self.random_seed = kwargs['random_seed']
      print 'using random number seed of %d' % self.random_seed

    self.dir = {}
    self.dir['start'] = os.getcwd()

    if kwargs['dir_CD'] is not None:
      self.dir['CD'] = os.path.abspath(kwargs['dir_CD'])
    else:
      self.dir['CD'] = os.path.abspath('.')

    if kwargs['dir_BC'] is not None:
      self.dir['BC'] = os.path.abspath(kwargs['dir_BC'])
    else:
      self.dir['BC'] = self.dir['CD']  # Default that may be
      # overwritten by stored directory

    # Load previously stored file names and arguments
    FNs = OrderedDict()
    args = OrderedDict()
    for p in ['CD', 'BC']:
      params = self._load_pkl_gz(p, kwargs['pose'])
      if params is not None:
        (fn_dict, param_dict) = params
        FNs[p] = dictionary_tools.convert_dictionary_relpath(
          fn_dict, relpath_o=self.dir[p], relpath_n=None)
        args[p] = param_dict
        if (p=='CD') and (kwargs['dir_BC'] is None) and \
           ('dir_BC' in FNs[p].keys()) and \
           (FNs[p]['dir_BC'] is not None):
          self.dir['BC'] = FNs[p]['dir_BC']
      else:
        FNs[p] = OrderedDict()
        args[p] = OrderedDict()

    print '\n*** Directories ***'
    print dictionary_tools.dict_view(self.dir)

    # Identify tarballs
    tarFNs = [kwargs[prefix + '_tarball'] \
      for prefix in ['ligand','receptor','complex'] \
      if (prefix + '_tarball') in kwargs.keys() and
      kwargs[(prefix + '_tarball')] is not None]
    for p in ['BC', 'CD']:
      if (p in FNs.keys()) and ('tarball' in FNs[p].keys()):
        tarFNs += [tarFN for tarFN in FNs[p]['tarball'].values() \
          if tarFN is not None]
    tarFNs = set([FN for FN in tarFNs if os.path.isfile(FN)])

    # Identify files to look for in the tarballs
    seekFNs = []
    if len(tarFNs) > 0:
      # From the keyword arguments
      for prefix in ['ligand', 'receptor', 'complex']:
        for postfix in ('database', 'prmtop', 'inpcrd', 'mol2', 'fixed_atoms'):
          key = '%s_%s' % (prefix, postfix)
          if (key in kwargs.keys()) and (kwargs[key] is not None):
            FN = os.path.abspath(kwargs[key])
            if not os.path.isfile(FN):
              seekFNs.append(os.path.basename(FN))
      if kwargs['score'] != 'default':
        seekFNs.append(kwargs['score'])
      # From files in a previous instance
      for p in ['BC', 'CD']:
        if p in FNs.keys():
          for level1 in ['ligand_database','receptor_database', \
              'prmtop','inpcrd','fixed_atoms']:
            if level1 in FNs[p].keys():
              if isinstance(FNs[p][level1], dict):
                for level2 in ['L', 'R', 'RL']:
                  if level2 in FNs[p][level1].keys():
                    seekFNs.append(os.path.basename(FNs[p][level1][level2]))
              else:
                seekFNs.append(os.path.basename(FNs[p][level1]))
      seekFNs = set(seekFNs)
    seek_frcmod = (kwargs['frcmodList'] is None) or \
      (not os.path.isfile(kwargs['frcmodList'][0]))

    if kwargs['keep_tar']:
      print 'Files extracted from tarballs will be kept\n'

    # Decompress tarballs into self.dir['CD']
    self.toClear = []

    if len(seekFNs) > 0:
      import tarfile

      print ">>> Decompressing tarballs"
      print 'looking for:\n  ' + '\n  '.join(seekFNs)
      if seek_frcmod:
        print '  and frcmod files'

      for tarFN in tarFNs:
        print 'reading ' + tarFN
        tarF = tarfile.open(tarFN, 'r')
        for member in tarF.getmembers():
          for seekFN in seekFNs:
            if member.name.endswith(seekFN):
              tarF.extract(member, path=self.dir['CD'])
              if not kwargs['keep_tar']:
                self.toClear.append(os.path.join(self.dir['CD'], seekFN))
              print '  extracted ' + seekFN
          if seek_frcmod and member.name.endswith('frcmod'):
            FN = os.path.abspath(os.path.join(self.dir['CD'], member.name))
            if not os.path.isfile(FN):
              tarF.extract(member, path=self.dir['CD'])
              kwargs['frcmodList'] = [FN]
              if not kwargs['keep_tar']:
                self.toClear.append(FN)
              print '  extracted ' + FN
      print

    # Set up file name dictionary
    print '*** Files ***'

    for p in ['BC', 'CD']:
      if p in FNs.keys():
        if FNs[p] != {}:
          print 'previously stored in %s directory:' % p
          print dictionary_tools.dict_view(FNs[p], relpath=self.dir['start'])

    if not (FNs['BC'] == {} and FNs['CD'] == {}):
      print 'from arguments and defaults:'

    def cdir_or_dir_CD(FN):
      if FN is not None:
        return path_tools.findPath([FN, os.path.join(self.dir['CD'], FN)])
      else:
        return None

    if kwargs['frcmodList'] is not None:
      if isinstance(kwargs['frcmodList'], str):
        kwargs['frcmodList'] = [kwargs['frcmodList']]
      kwargs['frcmodList'] = [cdir_or_dir_CD(FN) \
        for FN in kwargs['frcmodList']]

    if 'score' in kwargs.keys() and \
        (kwargs['score'] is not None) and \
        (kwargs['score'] != 'default'):
      kwargs['score'] = path_tools.findPath([kwargs['score'], \
        os.path.join(self.dir['CD'],kwargs['score'])])

    FFpath = path_tools.search_paths['gaff'] \
      if 'gaff' in path_tools.search_paths.keys() else []
    FNs['new'] = OrderedDict([
      ('ligand_database',cdir_or_dir_CD(kwargs['ligand_database'])),
      ('receptor_database',cdir_or_dir_CD(kwargs['receptor_database'])),
      ('forcefield',path_tools.findPath(\
        [kwargs['forcefield'],'../Data/gaff2.dat'] + FFpath)),
      ('frcmodList',kwargs['frcmodList']),
      ('tarball',OrderedDict([
        ('L',path_tools.findPath([kwargs['ligand_tarball']])),
        ('R',path_tools.findPath([kwargs['receptor_tarball']])),
        ('RL',path_tools.findPath([kwargs['complex_tarball']]))])),
      ('prmtop',OrderedDict([
        ('L',cdir_or_dir_CD(kwargs['ligand_prmtop'])),
        ('R',cdir_or_dir_CD(kwargs['receptor_prmtop'])),
        ('RL',cdir_or_dir_CD(kwargs['complex_prmtop']))])),
      ('inpcrd',OrderedDict([
        ('L',cdir_or_dir_CD(kwargs['ligand_inpcrd'])),
        ('R',cdir_or_dir_CD(kwargs['receptor_inpcrd'])),
        ('RL',cdir_or_dir_CD(kwargs['complex_inpcrd']))])),
      ('mol2',OrderedDict([
        ('L',cdir_or_dir_CD(kwargs['ligand_mol2']))])),
      ('fixed_atoms',OrderedDict([
        ('R',cdir_or_dir_CD(kwargs['receptor_fixed_atoms'])),
        ('RL',cdir_or_dir_CD(kwargs['complex_fixed_atoms']))])),
      ('grids',OrderedDict([
        ('LJr',path_tools.findPath([kwargs['grid_LJr'],
          os.path.join(kwargs['dir_grid'],'LJr.nc'),
          os.path.join(kwargs['dir_grid'],'LJr.dx'),
          os.path.join(kwargs['dir_grid'],'LJr.dx.gz')])),
        ('LJa',path_tools.findPath([kwargs['grid_LJa'],
          os.path.join(kwargs['dir_grid'],'LJa.nc'),
          os.path.join(kwargs['dir_grid'],'LJa.dx'),
          os.path.join(kwargs['dir_grid'],'LJa.dx.gz')])),
        ('sELE',path_tools.findPath([kwargs['grid_sELE'],
          kwargs['grid_ELE'],
          os.path.join(kwargs['dir_grid'],'pb.nc'),
          os.path.join(kwargs['dir_grid'],'pbsa.nc'),
          os.path.join(kwargs['dir_grid'],'direct_ELE.nc')])),
        ('ELE',path_tools.findPath([kwargs['grid_ELE'],
          os.path.join(kwargs['dir_grid'],'direct_ELE.nc'),
          os.path.join(kwargs['dir_grid'],'pb.nc'),
          os.path.join(kwargs['dir_grid'],'pbsa.nc')])),
        ('desolv',path_tools.findPath([kwargs['grid_desolv'],
          os.path.join(kwargs['dir_grid'],'desolv.nc'),
          os.path.join(kwargs['dir_grid'],'desolv.dx'),
          os.path.join(kwargs['dir_grid'],'desolv.dx.gz')]))])),
      ('score',kwargs['score']),
      ('dir_BC',self.dir['BC'])])

    if not (FNs['BC'] == {} and FNs['CD'] == {}):
      print dictionary_tools.dict_view(FNs['new'], relpath=self.dir['start'])
      print 'to be used:'

    self.FNs = dictionary_tools.merge_dictionaries(
      [FNs[src] for src in ['new', 'BC', 'CD']])

    # Default: a force field modification is in the same directory as the ligand
    if (self.FNs['frcmodList'] is None):
      if self.FNs['prmtop']['L'] is not None:
        dir_lig = os.path.dirname(self.FNs['prmtop']['L'])
        frcmodpaths = [os.path.abspath(os.path.join(dir_lig, \
          os.path.basename(self.FNs['prmtop']['L'])[:-7]+'.frcmod'))]
      else:
        dir_lig = '.'
        frcmodpaths = []
      if kwargs['frcmodList'] is None:
        frcmodpaths.extend([\
          os.path.abspath(os.path.join(dir_lig,'lig.frcmod')),\
          os.path.abspath(os.path.join(dir_lig,'ligand.frcmod'))])
        frcmod = path_tools.findPath(frcmodpaths)
        self.FNs['frcmodList'] = [frcmod]
    elif not isinstance(self.FNs['frcmodList'], list):
      self.FNs['frcmodList'] = [self.FNs['frcmodList']]

    # Check for existence of required files
    do_CD = (hasattr(args,'run_type') and \
              (args.run_type not in ['store_params', 'BC']))

    for key in ['ligand_database', 'forcefield']:
      if (self.FNs[key] is None) or (not os.path.isfile(self.FNs[key])):
        raise Exception('File for %s is missing!' % key)

    for (key1, key2) in [('prmtop', 'L'), ('inpcrd', 'L')]:
      FN = self.FNs[key1][key2]
      if (FN is None) or (not os.path.isfile(FN)):
        raise Exception('File for %s %s is missing' % (key1, key2))

    for (key1,key2) in [\
        ('prmtop','RL'), ('inpcrd','RL'), \
        ('grids','LJr'), ('grids','LJa'), ('grids','ELE')]:
      FN = self.FNs[key1][key2]
      errstring = 'Missing file %s %s required for CD!' % (key1, key2)
      if (FN is None) or (not os.path.isfile(FN)):
        if do_CD:
          raise Exception(errstring)
        else:
          print errstring

    if ((self.FNs['inpcrd']['RL'] is None) and \
        (self.FNs['inpcrd']['R'] is None)):
      if do_CD:
        raise Exception('Receptor coordinates needed for CD!')
      else:
        print 'Receptor coordinates needed for CD!'

    print dictionary_tools.dict_view(self.FNs,
                                     relpath=self.dir['start'],
                                     show_None=True)

    args['default_BC'] = OrderedDict([
      ('protocol', 'Adaptive'), ('therm_speed', 30.0), ('T_HIGH', 600.),
      ('T_SIMMIN', 300.), ('T_TARGET', 300.),
      ('H_mass', 4.0), ('delta_t', 4.0), ('sampler', 'NUTS'),
      ('steps_per_seed', 1000), ('seeds_per_state', 50), ('darts_per_seed', 0),
      ('repX_cycles', 20), ('min_repX_acc', 0.4), ('sweeps_per_cycle', 1000),
      ('snaps_per_cycle', 50), ('attempts_per_sweep', 25),
      ('steps_per_sweep', 50), ('darts_per_sweep', 0),
      ('phases', ['NAMD_Gas', 'NAMD_OBC']),
      ('sampling_importance_resampling', False), ('solvation', 'Desolvated'),
      ('keep_intermediate', False), ('GMC_attempts', 0),
      ('GMC_tors_threshold', 0.0)
    ])

    args['default_CD'] = OrderedDict(args['default_BC'].items() + [
      ('temperature_scaling','Linear'),
      ('site',None),
      ('site_center',None),
      ('site_direction',None),
      ('site_max_Z',None),
      ('site_max_R',None),
      ('site_density',50.),
      ('site_measured',None),
      ('pose',-1),
      ('k_pose', 1000.0), #  * MMTK.Units.kJ / MMTK.Units.mol / MMTK.Units.K
      ('MCMC_moves',1),
      ('rmsd',False)] + \
      [('receptor_'+phase,None) for phase in arguments.allowed_phases])
    args['default_CD']['snaps_per_cycle'] = 50

    # Store passed arguments in dictionary
    for p in ['BC', 'CD']:
      args['new_' + p] = OrderedDict()
      for key in args['default_' + p].keys():
        specific_key = p + '_' + key
        if (specific_key in kwargs.keys()) and \
           (kwargs[specific_key] is not None):
          # Use the specific key if it exists
          args['new_' + p][key] = kwargs[specific_key]
        elif (key in ['site_center', 'site_direction'] +
                     ['receptor_'+phase for phase in arguments.allowed_phases]) and \
             (kwargs[key] is not None):
          # Convert these to arrays of floats
          args['new_' + p][key] = np.array(kwargs[key], dtype=float)
        elif key in kwargs.keys():
          # Use the general key
          args['new_' + p][key] = kwargs[key]

    self.params = OrderedDict()
    for p in ['BC', 'CD']:
      self.params[p] = dictionary_tools.merge_dictionaries(
        [args[src] for src in ['new_' + p, p, 'default_' + p]])

    # Check that phases are permitted
    for phase in (self.params['BC']['phases'] + self.params['CD']['phases']):
      if phase not in arguments.allowed_phases:
        raise Exception(phase + ' phase is not supported!')

    # Make sure prerequistite phases are included:
    #   sander_Gas is necessary for any sander or gbnsr6 phase
    #   NAMD_Gas is necessary for APBS_PBSA
    for process in ['BC', 'CD']:
      phase_list = self.params[process]['phases']
      if (not 'sander_Gas' in phase_list) and \
          len([p for p in phase_list \
            if p.startswith('sander') or p.startswith('gbnsr6')])>0:
        phase_list.append('sander_Gas')
      if (not 'NAMD_Gas' in phase_list) and ('APBS_PBSA' in phase_list):
        phase_list.append('NAMD_Gas')

    # Variables dependent on the parameters
    self.original_Es = [[{}]]
    for phase in arguments.allowed_phases:
      if self.params['CD']['receptor_' + phase] is not None:
        self.original_Es[0][0]['R'+phase] = \
          np.atleast_2d(self.params['CD']['receptor_'+phase])
      else:
        self.original_Es[0][0]['R' + phase] = None

    print '\n*** Simulation parameters and constants ***'
    for p in ['BC', 'CD']:
      print '\nfor %s:' % p
      print dictionary_tools.dict_view(self.params[p])[:-1]

  def _load_pkl_gz(self, p, pose):
    """Loads parameters from a pickled gzip file

    """
    if p == 'CD' and pose > -1:
      progress_FN = os.path.join(self.dir[p],
                                 '%s_progress_pose%03d.pkl.gz' % (p, pose))
    else:
      progress_FN = os.path.join(self.dir[p], '%s_progress.pkl.gz' % (p))

    saved = load_pkl_gz(progress_FN)
    if (saved is None):
      if os.path.isfile(progress_FN):
        os.remove(progress_FN)
      if p == 'CD' and pose > -1:
        progress_FN = os.path.join(
          self.dir[p], '%s_progress_pose%03d.pkl.gz.BAK' % (p, pose))
      else:
        progress_FN = os.path.join(self.dir[p], '%s_progress.pkl.gz.BAK' % (p))

      saved = load_pkl_gz(progress_FN)
      if (saved is None):
        print '  no progress information for %s' % p
      else:
        print '  using stored progress and data in %s' % p
    # self._clear(p)

    params = None
    if saved is not None:
      params = saved[0]
    return params

  def _save_pkl_gz(self, p, data):
    """
    Saves simulation parameters, simulation files, and the protocol.

    Parameters
    ----------
    p : str
      Process, either 'BC' or 'CD'
    data : simulation_data.SimulationData object
      Simulation data, to store the protocol
    """
    random_orient = None
    if p == 'CD' and hasattr(self, '_n_trans'):
      random_orient = (self._n_trans, self._max_n_trans, self._random_trans, \
         self._n_rot, self._max_n_rot, self._random_rotT)

    param_dict = dict([tp for tp in self.params[p].items() \
                      if not tp[0] in ['repX_cycles']])
    if p == 'BC':
      fn_dict = dictionary_tools.convert_dictionary_relpath(
        {
          'ligand_database': self.FNs['ligand_database'],
          'forcefield': self.FNs['forcefield'],
          'frcmodList': self.FNs['frcmodList'],
          'tarball': {
            'L': self.FNs['tarball']['L']
          },
          'prmtop': {
            'L': self.FNs['prmtop']['L']
          },
          'inpcrd': {
            'L': self.FNs['inpcrd']['L']
          }
        },
        relpath_o=None,
        relpath_n=self.dir['BC'])
    elif p == 'CD':
      fn_dict = dictionary_tools.convert_dictionary_relpath(
        dict(self.FNs.items()), relpath_o=None, relpath_n=self.dir['CD'])
    params = (fn_dict, param_dict)

    saved = {'progress': (params, data.protocol, data.cycle)}

    key = 'progress'
    if p == 'CD' and self.params['CD']['pose'] > -1:
      saved_FN = os.path.join(self.dir[p],'%s_%s_pose%03d.pkl.gz'%(\
        p, key, self.params['CD']['pose']))
    else:
      saved_FN = os.path.join(self.dir[p], '%s_%s.pkl.gz' % (p, key))
    if not os.path.isdir(self.dir[p]):
      os.system('mkdir -p ' + self.dir[p])
    if os.path.isfile(saved_FN):
      os.rename(saved_FN, saved_FN + '.BAK')
    return write_pkl_gz(saved_FN, saved[key]) + '\n  saved %s progress' % p
