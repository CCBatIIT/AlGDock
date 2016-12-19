# NAMD_OBC corresponds to igb=5 in AMBER. It is the fastest GBSA implicit solvent.
#
# For sander and OpenMM, GBSA implicit solvents
# are labeled according to OpenMM nomenclature
#    HCT  - Hawkins-Cramer-Truhlar GBSA model (igb=1 in AMBER)
#    OBC1 - Onufriev-Bashford-Case GBSA model
#           using the GBOBCI parameters (igb=2 in AMBER)
#    OBC2 - Onufriev-Bashford-Case GBSA model
#           using the GBOBCII parameters (igb=5 in AMBER).
#    GBn  - GBn solvation model (igb=7 in AMBER).
#    GBn2 - GBn2 solvation model (igb=8 in AMBER).
# Several sander phases also have an ALPB option, which is applied for a receptor.
# ALPB is based on a spherical solute, which doesn't make sense for a ligand.
#
# GBNSR6 has two implicit phases.
#    Still - The original equation.
#    CHA   - Uses the charge hydration asymmetry correction.
# ALPB is used for anything with a receptor, as in sander.
#
# APBS implements the Poisson-Boltzmann implicit solvent

allowed_phases = ['NAMD_Gas','NAMD_OBC'] + \
  ['sander_'+p for p in ['Gas','HCT','OBC1','OBC2','GBn','GBn2','PBSA']] + \
  ['sander_ALPB_'+p for p in ['HCT','OBC1','OBC2','GBn']] + \
  ['gbnsr6_'+p for p in ['Still','CHA']] + \
  ['APBS_PBSA']

try:
  import simtk.openmm
  allowed_phases += ['OpenMM_'+p for p in ['Gas','HCT','OBC1','OBC2','GBn','GBn2']]
except ImportError:
  pass

arguments = {
  'dir_dock':{'help':'Directory where docking results are stored'},
  'dir_cool':{'help':'Directory where cooling results are stored'},
  'vmd':{'help':'Location of Visual Molecular Dynamics (VMD)'},
  'convert':{'help':'Location of convert (from ImageMagick)'},
  'font':{'help':'Location of font file (readable by PIL)'},
  #   Stored in both dir_cool and dir_dock
  'ligand_tarball':{'help':'tar file that contains database, prmtop, and inpcrd for ligand'},
  'ligand_database':{'help':'MMTK molecule definition for ligand'},
  'forcefield':{'help':'AMBER force field file'},
  'frcmodList':{'help':'AMBER force field modifications file(s)', 'nargs':'+'},
  'ligand_prmtop':{'help':'AMBER prmtop for the ligand'},
  'ligand_inpcrd':{'help':'AMBER coordinates for the ligand'},
  'ligand_mol2':{'help':'mol2 file for the ligand'},
  'ligand_rb':{'help':'rigid body specification for the ligand'},
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
  'run_type':{'choices':['configuration_energies','minimized_configuration_energies',
              'store_params', 'initial_cool', 'cool', \
              'dock','timed','postprocess',\
              'redo_postprocess','free_energies','redo_free_energies', 'all', \
              'render_docked', 'render_intermediates', \
              'clear_intermediates', None],
    'help':'Type of calculation to run'},
  'max_time':{'type':int, 'default':180, \
    'help':'For timed calculations, the maximum amount of wall clock time, in minutes'},
  'cores':{'type':int, 'help':'Number of CPU cores to use'},
  'rotate_matrix':{'help':'Rotation matrix for viewing'},
  'random_seed':{'type':int, 'help':'Random number seed'},
  #   Defaults
  'protocol':{'choices':['Adaptive','Geometric'],
    'help':'Approach to determining series of thermodynamic states'},
  'therm_speed':{'type':float,
    'help':'Thermodynamic speed during adaptive simulation'},
  'sampler':{
    'choices':['MixedHMC','HMC','NUTS','VV'],
    'help':'Sampling method'},
  'MCMC_moves':{'type':int,
    'help':'Types of MCMC moves to use'},
  'fraction_TD':{'type':float, 'default':0.5, \
    'help':'For the MixedHMC integrator, the fraction of MD steps spent with torsional MD'},
  'TD_steps_per_trial':{'type':int, 'default':5, \
    'help':'The number of steps for torsional dynamics HMC trials'},
  'delta_t_TD':{'type':float, 'default':4.0, \
    'help':'The time step for torsional dynamics HMC trials'},
  'delta_t':{'type':float, 'default':3.0, \
    'help':'The default time step, in fs'},
  'T_HIGH':{'type':float, 'default':600.0,
    'help':'High temperature'},
  'T_SIMMIN':{'type':float, 'default':325.0,
    'help':'Minimum simulation temperature'},
  'T_TARGET':{'type':float, 'default':300.0,
    'help':'Target temperature'},
  'H_mass':{'type':float, 'default':4.0, \
    'help':'The repartitioned mass of hydrogen. Set negative to turn off HMR'},
  # For initialization
  'seeds_per_state':{'type':int,
    'help':'Number of starting configurations in each state during initialization'},
  'steps_per_seed':{'type':int,
    'help':'Number of MD steps per state during initialization'},
  'darts_per_seed':{'type':int,
    'help':'Number of smart darting attempts for each seed during initialization'},
  # For replica exchange
  'repX_cycles':{'type':int,
    'help':'Number of replica exchange cycles for docking and cooling'},
  'sweeps_per_cycle':{'type':int,
    'help':'Number of replica exchange sweeps per cycle'},
  'attempts_per_sweep':{'type':int,
    'help':'Number of replica exchange attempts per sweep'},
  'steps_per_sweep':{'type':int,
    'help':'Number of MD steps per replica exchange sweep'},
  'darts_per_sweep':{'type':int,
    'help':'Number of smart darting attempts per replica exchange sweep'},
  'snaps_per_independent':{'type':int,
    'help':'Number of snapshots per independent sample'},
  'sampling_importance_resampling':{'action':'store_true',
    'help':'perfom sampling importance resampling'}, # TODO: Set new default after testing
  'keep_intermediate':{'action':'store_true',
    'help':'Keep configurations for intermediate states'},
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
    'help':'Density of center-of-mass points in the first docking state'},
  # Binding pose
  'pose':{'type':int, 'default':-1,
    'help':'Confine the ligand to a specific pose from the "score" file. If the argument is negative there is no pose restriction. Otherwise the argument is the index of the pose.'},
  'k_pose':{'type':float, \
    'help':'Spring constant for restraint to the binding pose.'}}

import copy
for process in ['cool','dock']:
  for key in ['protocol', 'therm_speed', 'sampler',
      'seeds_per_state', 'steps_per_seed', 'darts_per_seed',
      'sweeps_per_cycle', 'attempts_per_sweep',
      'steps_per_sweep', 'darts_per_sweep',
      'snaps_per_independent', 'keep_intermediate']:
    arguments[process+'_'+key] = copy.deepcopy(arguments[key])

for phase in allowed_phases:
  arguments['receptor_'+phase] = {'type':float, 'nargs':'+',
    'help':'Receptor potential energy in %s (kJ/mol)'%phase}
