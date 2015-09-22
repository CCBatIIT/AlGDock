allowed_phases = ['Gas','GBSA','PBSA','NAMD_Gas','NAMD_GBSA',\
  'OpenMM_Gas','OpenMM_GBn','OpenMM_GBn2','OpenMM_OBC1','OpenMM_OBC2','OpenMM_HCT',\
  'APBS']

arguments = {
  'dir_dock':{'help':'Directory where docking results are stored'},
  'dir_cool':{'help':'Directory where cooling results are stored'},
  'namd':{'help':'Location of Not Another Molecular Dynamics (NAMD)'},
  'vmd':{'help':'Location of Visual Molecular Dynamics (VMD)'},
  'sander':{'help':'Location of sander (from AMBER)'},
  'convert':{'help':'Location of convert (from ImageMagick)'},
  'font':{'help':'Location of font file (readable by PIL)'},
  #   Stored in both dir_cool and dir_dock
  'ligand_tarball':{'help':'tar file that contains database, prmtop, and inpcrd for ligand'},
  'ligand_database':{'help':'MMTK molecule definition for ligand'},
  'forcefield':{'help':'AMBER force field file'},
  'frcmodList':{'help':'AMBER force field modifications file(s)', 'nargs':'+'},
  'ligand_prmtop':{'help':'AMBER prmtop for the ligand'},
  'ligand_inpcrd':{'help':'AMBER coordinates for the ligand'},
  #   Stored in dir_dock
  'receptor_tarball':{'help':'tar file that contains prmtop, inpcrd, and fixed_atoms files for the receptor'},
  'receptor_prmtop':{'help':'AMBER prmtop for the receptor'},
  'receptor_inpcrd':{'help':'AMBER coordinates for the receptor'},
  'receptor_fixed_atoms':{'help':'PDB file with fixed atoms labeled by 1 in the occupancy column'},
  'complex_tarball':{'help':'tar file that contains prmtop, inpcrd, and fixed_atoms files for the complex'},
  'complex_prmtop':{'help':'AMBER prmtop file for the complex'},
  'complex_inpcrd':{'help':'AMBER coordinates for the complex'},
  'complex_fixed_atoms':{'help':'PDB file with fixed atoms labeled by 1 in the occupancy column'},
  'dir_grid':{'help':'Directory containing potential energy grids'},
  'grid_LJr':{'help':'DX file for Lennard-Jones repulsive grid'},
  'grid_LJa':{'help':'DX file for Lennard-Jones attractive grid'},
  'grid_ELE':{'help':'DX file for electrostatic grid'},
  'score':{'help':"Starting configuration(s) for replica exchange. Can be a mol2 file, .pkl.gz file, or 'default'"},
  # Simulation settings and constants
  #   Run-dependent 
  'cool_repX_cycles':{'type':int,
    'help':'Number of replica exchange cycles for cooling'},
  'dock_repX_cycles':{'type':int,
    'help':'Number of replica exchange cycles for docking'},
  'run_type':{'choices':['pose_energies','minimized_pose_energies',
              'store_params', 'cool', \
              'dock','timed','postprocess',\
              'redo_postprocess','free_energies','redo_free_energies', 'all', \
              'render_docked', 'render_intermediates', \
              'clear_intermediates', None],
    'help':'Type of calculation to run'},
  'max_time':{'type':int, 'default':180, \
    'help':'For timed calculations, the maximum amount of wall clock time, in minutes'},
  'cores':{'type':int, \
    'help':'Number of CPU cores to use'},
  'rotate_matrix':{'help':'Rotation matrix for viewing'},
  #   Defaults
  'protocol':{'choices':['Adaptive','Set'],
    'help':'Approach to determining series of thermodynamic states'},
  'therm_speed':{'type':float,
    'help':'Thermodynamic speed during adaptive simulation'},
  'sampler':{
    'choices':['HMC','NUTS','VV','TDHMC'],
    'help':'Sampling method'},
  'MCMC_moves':{'type':int,
    'help':'Types of MCMC moves to use'},
  'T_HIGH':{'type':float, 'default':600.0,
    'help':'High temperature'},
  'T_TARGET':{'type':float, 'default':300.0,
    'help':'Target temperature'},
  # For initialization
  'seeds_per_state':{'type':int,
    'help':'Number of starting configurations in each state during initialization'},
  'steps_per_seed':{'type':int,
    'help':'Number of MD steps per state during initialization'},
  # For replica exchange
  'repX_cycles':{'type':int,
    'help':'Number of replica exchange cycles for docking and cooling'},
  'sweeps_per_cycle':{'type':int,
    'help':'Number of replica exchange sweeps per cycle'},
  'attempts_per_sweep':{'type':int,
    'help':'Number of replica exchange attempts per sweep'},
  'steps_per_sweep':{'type':int,
    'help':'Number of MD steps per replica exchange sweep'},
  'snaps_per_independent':{'type':int,
    'help':'Number of snapshots per independent sample'},
  'keep_intermediate':{'action':'store_true',
    'help':'Keep configurations for intermediate states?'},
  'min_repX_acc':{'type':float,
    'help':'Minimum value for replica exchange acceptance rate'},
  # for GMC
  'GMC_attempts':{'type':int, 'default': 0,
    'help': 'Number of attempts is K * GMC_attempts. Zero means not to do GMC' },
  'GMC_tors_threshold': {'type':float, 'default': 0.0,
    'help': 'The torsion threshold (in radian) below which no crossover will be attempted' },
  # For postprocessing
  'phases':{'nargs':'+', 'help':'Phases to use in postprocessing'},
  'rmsd':{'nargs':'?', 'const':True, 'default':False,
    'help':'Calculate rmsd between snapshots and a configuration or set of configurations, which may be passed as an argument. The default configuration is from the ligand_inpcrd argument.'},
  # Binding site
  'site':{'choices':['Sphere','Cylinder','Measure'], 'default':'Sphere', \
    'help':'Type of binding site. "Measure" means that parameters' + \
           ' for a sphere will be measured from docked poses.'},
  'site_center':{'nargs':3, 'type':float,
    'help':'Position of binding site center'},
  'site_direction':{'nargs':3, 'type':float,
    'help':'Principal axis of a cylindrical binding site'},
  'site_max_X':{'type':float,
    'help':'Maximum position along principal axis in a cylindrical binding site'},
  'site_max_R':{'type':float,
    'help':'Maximum radial position for a spherical or cylindrical binding site'},
  'site_density':{'type':float,
    'help':'Density of center-of-mass points in the first docking state'}}

import copy
for process in ['cool','dock']:
  for key in ['protocol', 'therm_speed', 'sampler',
      'seeds_per_state', 'steps_per_seed',
      'sweeps_per_cycle', 'attempts_per_sweep', 'steps_per_sweep',
      'snaps_per_independent', 'keep_intermediate']:
    arguments[process+'_'+key] = copy.deepcopy(arguments[key])

for phase in allowed_phases:
  arguments['receptor_'+phase] = {'type':float, 'nargs':'+',
    'help':'Receptor potential energy in %s (kJ/mol)'%phase}
