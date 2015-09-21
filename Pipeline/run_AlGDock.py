#!/usr/bin/env python

import argparse
parser = argparse.ArgumentParser()

parser.add_argument('--saved_arguments', default='saved_arguments.py',
  help='File containing default values of parameters' + \
       '(can be overwritten by flags)')
parser.add_argument('--include_ligand', nargs='+', default=None, \
  help='Only runs AlGDock for these ligands')
parser.add_argument('--exclude_ligand', nargs='+', default=None, \
  help='Does not run AlGDock for these ligands')
parser.add_argument('--new_instance_per_ligand', action='store_true', \
  default=False, \
  help='Runs run_AlGDock.py for each ligand')
parser.add_argument('--include_receptor', nargs='+', default=None, \
  help='Only runs AlGDock for these receptors')
parser.add_argument('--exclude_receptor', nargs='+', default=None, \
  help='Does not run AlGDock for these receptor')
parser.add_argument('--older_than', type=int, default=None, \
  help='Only runs redo_free_energies jobs if f_RL.pkl.gz ' + \
       'is older than OLDER_THAN hours')
parser.add_argument('--check_complete', action='store_true', default=False, \
  help='Checks whether f_RL.pkl.gz is complete')

# Arguments related to file locations
parser.add_argument('--ligand', default='../ligand/AlGDock_in/', \
  help='The directory/file tree to look for ligand files')
parser.add_argument('--library_requirement', default=None, \
  help='The ligand file name must contain the string LIBRARY_REQUIREMENT')
parser.add_argument('--receptor', default='../receptor/amber_in', \
  help='The directory/file to look for receptor files (prmtop and inpcrd)')
parser.add_argument('--receptor_grids', default='../receptor/AlGDock_in', \
  help='The directory to look for receptor grids (nc)')
parser.add_argument('--complex', default='../complex/AlGDock_in', \
  help='The directory tree to look for complex files (prmtop and inpcrd)')
parser.add_argument('--site_info', \
  default='../receptor/2-binding_site/measured_binding_site.py', \
  help='Python script with binding site parameters ' + \
       '(site_R, half_edge_length)')
parser.add_argument('--dock6', default='../dock6/', \
  help='The directory to look for dock6 results (*.mol2.gz)')
parser.add_argument('--tree_cool', default='cool/',
  help='Directory tree to store cooling results')
parser.add_argument('--tree_dock', default='dock/',
  help='Directory tree to store docking results')
# Arguments related to job management
parser.add_argument('--reps', default=None, nargs=2, type=int, \
  help='Range of repetitions')
parser.add_argument('--max_jobs', default=None, type=int)
parser.add_argument('--first_ligand', default=None, type=int)
parser.add_argument('--max_ligands', default=None, type=int)
parser.add_argument('--dry', action='store_true', default=False, \
  help='Does not actually submit the job to the queue')
parser.add_argument('--no_release', action='store_true', default=False, \
  help='Does not release held jobs')
parser.add_argument('--interactive', action='store_true', default=False, \
  help='Output command for running in an interactive python environment')
parser.add_argument('--check_tarballs', action='store_true', default=False, \
  help='Check inside tarballs for files')
parser.add_argument('--skip_onq', action='store_true', default=False, \
  help='Skips looking for the job on the queue')
# Arguments related to scoring and assessment
parser.add_argument('--score', choices=['xtal','dock',None], default=None,
  help='Does a scoring rather than ab initio docking. ' + \
  'xtal means the crystal structure and dock means from another docking program.')
parser.add_argument('--rmsd', choices=['xtal',None], default=None,
  help='Calulates the rmsd between snapshots a configuration. ' + \
       'xtal means the crystal structure.')

# Simulation settings and constants
#   Run-dependent 
parser.add_argument('--cool_repX_cycles', type=int,
  help='Number of replica exchange cycles for cooling')
parser.add_argument('--dock_repX_cycles', type=int,
  help='Number of replica exchange cycles for docking')
parser.add_argument('--run_type',
  choices=['pose_energies','minimized_pose_energies', \
           'store_params', 'cool', \
           'dock','timed','postprocess',\
           'redo_postprocess','free_energies','redo_free_energies', 'all', \
           'render_docked', 'render_intermediates', \
           'clear_intermediates', None],
  help='Type of calculation to run')
parser.add_argument('--max_time', type=int, default = 180, \
  help='For timed calculations, the maximum amount of wall clock time, in minutes')
parser.add_argument('--cores', type=int, \
  help='Number of CPU cores to use')
#   Defaults
parser.add_argument('--protocol', choices=['Adaptive','Set'],
  help='Approach to determining series of thermodynamic states')
parser.add_argument('--therm_speed', type=float,
  help='Thermodynamic speed during adaptive simulation')
parser.add_argument('--sampler',
  choices=['HMC','NUTS','VV'],
  help='Sampling method')
parser.add_argument('--MCMC_moves', type=int,
  help='Types of MCMC moves to use')
parser.add_argument('--T_HIGH', type=float,
  help='High temperature')
parser.add_argument('--T_TARGET', type=float,
  help='Target temperature')
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
parser.add_argument('--attempts_per_sweep', type=int,
  help='Number of replica exchange attempts per sweeps')
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
parser.add_argument('--cool_attempts_per_sweep', type=int,
  help='Number of replica exchange attempts per sweeps')
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
parser.add_argument('--dock_attempts_per_sweep', type=int,
  help='Number of replica exchange attempts per sweeps')
parser.add_argument('--dock_steps_per_seed', type=int,
  help='Number of MD steps per state during initialization')
# For replica exchange
parser.add_argument('--dock_sweeps_per_cycle', type=int,
  help='Number of replica exchange sweeps per cycle')
parser.add_argument('--dock_steps_per_sweep', type=int,
  help='Number of MD steps per replica exchange sweep')
parser.add_argument('--dock_keep_intermediate', action='store_true',
  default=None, help='Keep configurations for intermediate states?')

parser.add_argument('--site',
  choices=['Sphere','Cylinder','Measure'], \
  help='Type of binding site. "Measure" means that parameters' + \
         ' for a sphere will be measured from docked poses.')
parser.add_argument('--site_center', nargs=3, type=float,
  help='Position of binding site center')
#  parser.add_argument('--site_direction', nargs=3, type=float,
#    help='Principal axis of a cylindrical binding site')
#  parser.add_argument('--site_max_X', type=float,
#    help='Maximum position along principal axis in a cylindrical binding site')
parser.add_argument('--site_max_R', type=float,
  help='Maximum radial position for a spherical or cylindrical binding site')
parser.add_argument('--site_density', type=float,
  help='Density of center-of-mass points in the first docking stage')
# Additional calculations
parser.add_argument('--receptor_Gas', type=float, nargs='+',
  help='Receptor potential energies in AMBER Gas implicit solvent (in units of kJ/mol)')
parser.add_argument('--receptor_GBSA', type=float, nargs='+',
  help='Receptor potential energies in AMBER GBSA implicit solvent (in units of kJ/mol)')
parser.add_argument('--receptor_PBSA', type=float, nargs='+',
  help='Receptor potential energies in AMBER PBSA implicit solvent (in units of kJ/mol)')
parser.add_argument('--receptor_NAMD_Gas', type=float, nargs='+',
  help='Receptor potential energies in gas phase (in units of kJ/mol)')
parser.add_argument('--receptor_NAMD_GBSA', type=float, nargs='+',
  help='Receptor potential energies in NAMD GBSA implicit solvent (in units of kJ/mol)')
args_in = parser.parse_args()

# Check for the existence of input files
import os, glob

# Look for gaff.dat, AlGDock, and qsub_command.py
import inspect
dirs = {}
dirs['current'] = os.getcwd()
dirs['script'] = os.path.dirname(os.path.abspath(\
  inspect.getfile(inspect.currentframe())))
execfile(os.path.join(dirs['script'],'_external_paths.py'))

for path in ['ligand','receptor','receptor_grids','complex',\
             'tree_dock','tree_cool','dock6']:
  setattr(args_in,path,os.path.abspath(getattr(args_in,path)))

if not args_in.skip_onq:
  execfile(os.path.join(dirs['script'],'_jobs_on_queue.py'))
  onq = jobs_on_queue()
else:
  onq = []
command_paths = findPaths(['qsub_command','gaff.dat'])
algdock_path = findPath(search_paths['algdock'])

import numpy as np

def nonzero(path):
  return (os.path.isfile(path) and os.path.getsize(path)>0) or \
         (os.path.isdir(path))

# Get other parameters
general_sim_arg_keys = ['protocol', 'therm_speed',\
  'sampler', 'seeds_per_state', 'steps_per_seed', 'repX_cycles', \
  'sweeps_per_cycle', 'steps_per_sweep', 'keep_intermediate', 'phases']
sim_arg_keys = general_sim_arg_keys + \
  ['cool_'+a for a in general_sim_arg_keys] + \
  ['dock_'+a for a in general_sim_arg_keys] + \
  ['MCMC_moves', 'score', 'rmsd', 'run_type', 'cores'] + \
  ['site', 'site_center', 'site_direction', \
   'site_max_X', 'site_max_R', 'site_density','reps']

if (args_in.saved_arguments is not None) and nonzero(args_in.saved_arguments):
  print 'Passed arguments:'
  print 'Loading default simulation arguments from '+args_in.saved_arguments
  execfile(args_in.saved_arguments)
  for key in sim_arg_keys:
    if hasattr(args_saved,key) and \
      ((not hasattr(args_in,key)) or (getattr(args_in,key) is None)):
      setattr(args_in,key,getattr(args_saved,key))
else:
  class args_saved:
    pass

# Determine the binding site radius and half edge length
site = args_in.site
if hasattr(args_in,'site_center') and hasattr(args_in,'site_max_R'):
  site_center = args_in.site_center
  site_max_R = args_in.site_max_R
else:
  # If no argument, load in a file
  if os.path.isfile(args_in.site_info):
    execfile(args_in.site_info)
  else:
    raise Exception('No binding site information')
  # These should be in nanometers, not Angstroms
  site_center = [half_edge_length*0.1, half_edge_length*0.1, half_edge_length*0.1]
  site_max_R = site_R/10.

# Look for ligand files
if os.path.isfile(args_in.ligand):
  ligand_FNs = [os.path.abspath(args_in.ligand)]
elif os.path.isdir(args_in.ligand):
  ligand_FNs = glob.glob(os.path.join(args_in.ligand,'*/*.tar.gz'))
  ligand_FNs = sorted([os.path.abspath(FN) for FN in ligand_FNs])
else:
  raise Exception('Ligand input %s is not a file or directory!'%args_in.ligand)
# Check that ligand tarball has nonzero size
ligand_FNs = [FN for FN in ligand_FNs if os.path.getsize(FN)>0]
if len(ligand_FNs)>1:
  # Filter ligands
  if args_in.include_ligand is not None:
    ligand_FNs = [FN for FN in ligand_FNs \
      if np.array([FN.find(ligN)!=-1 \
      for ligN in args_in.include_ligand]).any()]
  if args_in.exclude_ligand is not None:
    ligand_FNs = [FN for FN in ligand_FNs \
      if not np.array([FN.find(ligN)!=-1 \
      for ligN in args_in.exclude_ligand]).any()]

  if args_in.new_instance_per_ligand:
    import sys
    for ligand_FN in ligand_FNs:
      command = ' '.join(sys.argv + ['--ligand',ligand_FN]) + ' &'
      print command
      os.system(command)
    sys.exit()

# Look for receptor files
if os.path.isfile(args_in.receptor):
  receptor_FNs = [args_in.receptor]
elif os.path.isdir(args_in.receptor):
  receptor_FNs = glob.glob(os.path.join(args_in.receptor,'*.prmtop'))
else:
  raise Exception('Receptor input %s is not a file or directory!'%args_in.receptor)
receptor_FNs = [os.path.abspath(FN) for FN in receptor_FNs]
# Filter receptors
if args_in.include_receptor is not None:
  receptor_FNs = [FN for FN in receptor_FNs \
    if np.array([FN.find(recN)!=-1 \
    for recN in args_in.include_receptor]).any()]
if args_in.exclude_receptor is not None:
  receptor_FNs = [FN for FN in receptor_FNs \
    if not np.array([FN.find(recN)!=-1 \
    for recN in args_in.exclude_receptor]).any()]

# Require inpcrd as well as prmtop files
receptor_FNs = [FN for FN in receptor_FNs if
  np.array([nonzero(FN[:-6]+key) \
    for key in ['prmtop','inpcrd']]).all()]
# Require grid files
receptor_FNs = [FN for FN in receptor_FNs if np.array(\
  [nonzero(os.path.join(args_in.receptor_grids,os.path.basename(FN)[:-6]+key)) \
    for key in ['LJa.25.nc','LJr.25.nc','PB.nc']]).all()]

# Look for complex files
if os.path.isfile(args_in.complex):
  complex_FNs = [args_in.complex]
elif os.path.isdir(args_in.complex):
  complex_FNs = glob.glob(os.path.join(args_in.complex,'*/*/*.tar.gz'))
else:
  raise Exception('Complex input %s is not a file or directory!'%args_in.complex)
# Require that complex tarball has nonzero size
complex_FNs = [FN for FN in complex_FNs if nonzero(FN)]

print 'Found %d ligands, %d receptors, and %d complexes ready for AlGDock'%(\
  len(ligand_FNs),len(receptor_FNs),len(complex_FNs))

if (args_in.library_requirement is not None):
  ligand_FNs = [FN for FN in ligand_FNs \
    if FN.find(args_in.library_requirement)>-1]
  print '%d ligand(s) meet the library requirement'%(len(ligand_FNs))

if args_in.reps is None:
  args_in.reps = [0,1]
if args_in.first_ligand is None:
  args_in.first_ligand = 0
if args_in.max_ligands is None:
  args_in.max_ligands = len(ligand_FNs)

print 'Arguments:'
print args_in

namespace = locals()

job_status = {'submitted':0, 'skipped':0, 'no_complex':0, 'no_poses':0, \
  'no_dock6':0, 'missing_file':0, 'onq':0, 'recently_redone':0, 'complete':0}

import tarfile
checked = []

for rep in range(args_in.reps[0],args_in.reps[1]):
  for ligand_FN in ligand_FNs[args_in.first_ligand:args_in.first_ligand+args_in.max_ligands]:
    labels = {}
    if os.path.isfile(args_in.ligand):
      lib_and_key_prefix = ligand_FN.split('/')[-2]
    else:
      lib_and_key_prefix = os.path.dirname(ligand_FN[len(args_in.ligand)+1:])
    labels['library'] = '.'.join(lib_and_key_prefix.split('.')[:-1])
    labels['key'] = os.path.basename(ligand_FN[:-7])    
    labels['lib_subdir'] = labels['library']+'.'+labels['key'][:-2]+'__'
    labels['ligand'] = labels['library']+'.'+labels['key']
    paths = {'dir_cool':os.path.join(args_in.tree_cool, \
                labels['lib_subdir'], '%s-%d'%(labels['key'],rep)),
             'forcefield':command_paths['gaff.dat'],
             'ligand_tarball':ligand_FN}

    # Define and check files within the ligand tarball
    paths_in_tar = {'ligand_database':labels['ligand'].lower()+'.db'}
    for key in [('ligand_prmtop','prmtop'),('ligand_inpcrd','inpcrd'), \
                ('frcmodList','frcmod')]:
      paths_in_tar[key[0]] = labels['ligand']+'.'+key[1]
    if args_in.check_tarballs and (ligand_FN not in checked):
      tarF = tarfile.open(ligand_FN)
      names = [m.name for m in tarF.getmembers()]
      not_found = [paths_in_tar[key] for key in paths_in_tar.keys() \
        if not paths_in_tar[key] in names]
      if len(not_found)>0:
        print 'The following files were missing in '+ligand_FN+':'
        print ' '.join(not_found)
        continue
      else:
        checked.append(ligand_FN)
        
    if not os.path.isdir(paths['dir_cool']):
      os.system('mkdir -p '+paths['dir_cool'])
    if (args_in.run_type in ['initial_cool','cool']) and \
        nonzero(os.path.join(paths['dir_cool'],'f_L.pkl.gz')):
      job_status['complete'] += 1
      continue # Cooling is already done
    for receptor_FN in receptor_FNs:
      labels['receptor'] = os.path.basename(receptor_FN)[:-7]
      labels['complex'] = labels['library']+'.'+labels['key']+'-'+labels['receptor']
      # Identify receptor files
      for key in ['prmtop','inpcrd']:
        paths['receptor_'+key] = os.path.abspath(receptor_FN[:-6]+key)
      if os.path.isfile(os.path.abspath(receptor_FN[:-6]+'pdb')):
        paths['receptor_fixed_atoms'] = os.path.abspath(receptor_FN[:-6]+'pdb')
      for key in [('grid_LJa','LJa.25.nc'), ('grid_LJr','LJr.25.nc'), \
                  ('grid_ELE','PB.nc')]:
        paths[key[0]] = os.path.join(args_in.receptor_grids, \
          '%s.%s'%(labels['receptor'],key[1]))

      # Identify complex tarball
      complex_tar_FN = os.path.join(args_in.complex, \
        labels['lib_subdir'], labels['key'], labels['receptor']+'.tar.gz')
      if not (os.path.isfile(complex_tar_FN)):
        print 'No complex tarfile ' + complex_tar_FN
        job_status['no_complex'] += 1
        continue # Complex files are missing
      paths['complex_tarball'] = complex_tar_FN

      # Define and check files within the complex tarball
      for key in ['prmtop','inpcrd']:
        paths_in_tar['complex_'+key] = labels['complex']+'.'+key
      if 'receptor_fixed_atoms' in paths.keys():
        paths_in_tar['complex_fixed_atoms'] = labels['complex']+'.pdb'
      if args_in.check_tarballs and (complex_tar_FN not in checked):
        tarF = tarfile.open(complex_tar_FN)
        names = [m.name for m in tarF.getmembers()]
        not_found = [paths_in_tar[key] for key in paths_in_tar.keys() \
          if key.startswith('complex') and not paths_in_tar[key] in names]
        if len(not_found)>0:
          print 'The following files were missing in '+complex_tar_FN+':'
          print ' '.join(not_found)
          continue
        else:
          checked.append(complex_tar_FN)

      input_FNs = paths.values()
      input_FNs_missing = np.array([not nonzero(FN) for FN in input_FNs])
      if input_FNs_missing.any():
        print 'Necessary files:'
        print paths
        print 'Files are missing: ' + ', '.join(np.array(input_FNs)[input_FNs_missing])
        job_status['missing_file'] += 1
        continue # Files are missing

      # Convert relative path to absolute paths
      for key in paths.keys():
        paths[key] = os.path.abspath(paths[key])

      labels['job'] = '%s-%d'%(labels['complex'],rep)
      paths['dir_dock'] = os.path.join(args_in.tree_dock, \
        labels['lib_subdir'], labels['key'], '%s-%d'%(labels['receptor'],rep))
      if not os.path.isdir(paths['dir_dock']):
        os.system('mkdir -p '+paths['dir_dock'])
        
      f_RL_FN = os.path.join(paths['dir_dock'],'f_RL.pkl.gz')
      if nonzero(f_RL_FN):
        import time
        # Check if the calculation was recently redone
        if (args_in.run_type=='redo_free_energies') and \
            (args_in.older_than is not None) and \
            (time.time()-os.path.getmtime(f_RL_FN))/60./60.<args_in.older_than:
          job_status['recently_redone'] += 1
          continue

        if (args_in.run_type in \
            ['random_dock','initial_dock', 'dock','all','timed']):
          if args_in.check_complete:
            import gzip, pickle
            F = gzip.open(f_RL_FN,'r')
            dat = pickle.load(F)
            F.close()
            try:
              completed_cycles = np.min([len(dat[-1][p+'_MBAR']) for p in args_in.phases])
              complete = (completed_cycles >= int(args_in.dock_repX_cycles))
            except:
              print 'Error in '+f_RL_FN
              completed_cycles = 0
              complete = False
            if complete:
              # Docking is done
              job_status['complete'] += 1
              continue
            else:
              print '%d/%d cycles in %s'%(\
                completed_cycles, args_in.dock_repX_cycles, paths['dir_dock'])
          else:
              # Docking is done
              job_status['complete'] += 1
              continue
      elif (args_in.run_type=='redo_free_energies'):
        job_status['missing_file'] += 1
        continue
        
      jobname = '-'.join(dirs['current'].split('/')[-2:]+[labels['job']])
      if jobname in onq:
        job_status['onq'] += 1
        print jobname + ' is on the queue'
        continue # Job is on the queue

      # Consolidate paths to pass
      paths_to_pass = {}
      for key in paths.keys():
        paths_to_pass[key] = paths[key]
      for key in paths_in_tar.keys():
        paths_to_pass[key] = paths_in_tar[key]

      interactive_to_pass = []
      terminal_to_pass = []
      skip_job = False
      for key in (paths_to_pass.keys() + sim_arg_keys):
        # Priority is passed arguments (which may include saved arguments),
        #   local variables,
        #   and then the path dictionary
        if hasattr(args_in,key) and (getattr(args_in,key) is not None):
          val = getattr(args_in,key)
        elif (key in namespace.keys()) and (namespace[key] is not None):
          val = namespace[key]
        elif key in paths_to_pass.keys():
          val = paths_to_pass[key]
        else:
          continue
        # Special cases
        if key=='rmsd':
          val = {None: False, 'xtal':True}[val]
        elif key=='score':
          if val is None:
            val = False
          elif val=='xtal':
            val = 'default'
          elif val=='dock':
            val = os.path.abspath(os.path.join(args_in.dock6, \
              labels['lib_subdir'], labels['key'], labels['receptor'] + '.nc'))
            if not nonzero(val):
              if os.path.isfile(val[:-3]+'.mol2.gz'):
                job_status['no_poses'] += 1
                skip_job = True
                break # No poses in dock6
              else:
                print 'No dock6 output in '+val
                job_status['no_dock6'] += 1
                skip_job = True
                break # Dock6 files are missing
        elif key=='frcmodList':
          val = [val]
        # Actual strings to pass
        if isinstance(val,str):
          interactive_to_pass.append("%s = '%s'"%(key,val))
          terminal_to_pass.append("--%s %s"%(key,val))
        elif isinstance(val,bool):
          if val:
            interactive_to_pass.append("%s = %s"%(key,val))
            terminal_to_pass.append("--%s"%(key))
        elif isinstance(val,int):
          interactive_to_pass.append("%s = %d"%(key,val))
          terminal_to_pass.append("--%s %d"%(key,val))
        elif isinstance(val,float):
          interactive_to_pass.append("%s = %.5f"%(key,val))
          terminal_to_pass.append("--%s %.5f"%(key,val))
        elif isinstance(val,list):
          if isinstance(val[0],str):
            interactive_to_pass.append("%s = ['%s']"%(key,
              "', '".join([a for a in val])))
            terminal_to_pass.append("--%s %s"%(key,
              " ".join([a for a in val])))
          elif isinstance(val[0],float):
            interactive_to_pass.append("%s = [%s]"%(key,
              ", ".join(['%.5f'%a for a in val])))
            terminal_to_pass.append("--%s %s"%(key,
              " ".join(['%.5f'%a for a in val])))
        else:
          print 'Value:', val
          raise Exception('Type not known!')

      if skip_job:
        continue

      outputFNs = {}
      for FN in ['cool_log.txt',
          'cool_progress.pkl.gz','cool_progress.pkl.gz.BAK',
          'cool_data.pkl.gz','cool_data.pkl.gz.BAK',
          'f_L.pkl.gz']:
        outputFNs[FN] = os.path.join(paths_to_pass['dir_cool'],FN)
      for FN in ['dock_log.txt',
          'dock_progress.pkl.gz', 'dock_progress.pkl.gz.BAK',
          'dock_data.pkl.gz', 'dock_data.pkl.gz.BAK',
          'f_RL.pkl.gz']:
        outputFNs[FN] = os.path.join(paths_to_pass['dir_dock'],FN)
      for k in outputFNs.keys():
        if not os.path.isfile(outputFNs[k]):
          open(outputFNs[k], 'a').close()
      transfer_output_remaps = [f(key) \
        for key in outputFNs.keys() \
          for f in (lambda x:key,lambda x:outputFNs[key])]

      interactive_command = "self = AlGDock.BindingPMF.BPMF(" + \
        ', \\\n  '.join(interactive_to_pass) + \
        ")"
      if algdock_path is None:
        terminal_command = '$ALGDOCK '
      else:
        if algdock_path.find('.py')!=-1:
          terminal_command = 'python ' + algdock_path + ' '
        else:
          terminal_command = algdock_path + ' '
      terminal_command += ' \\\n  '.join(terminal_to_pass)

      if args_in.interactive:
        print interactive_command
      else:
        if args_in.run_type=='random_dock':
          mem=16
        else:
          mem=2
        if args_in.run_type in ['initial_cool','cool']:
          os.chdir(paths['dir_cool'])
        else:
          os.chdir(paths['dir_dock'])
        import subprocess
        subprocess.call(['python', command_paths['qsub_command'], \
          jobname, terminal_command, '--mem', '%d'%mem] + \
          ['--input_files'] + outputFNs.values() + \
          ['--output_files'] + list(outputFNs) + \
          ['--output_remaps'] + transfer_output_remaps + \
          ['--comment', interactive_command.replace(' \\\n','')] + \
          {True:['--dry'],False:[]}[args_in.dry] + \
          {True:['--no_release'],False:[]}[args_in.no_release])
        os.chdir(dirs['current'])
      
      job_status['submitted'] += 1

      if (args_in.max_jobs is not None) and \
         (job_status['submitted']>=args_in.max_jobs):
        break
    if (args_in.max_jobs is not None) and \
       (job_status['submitted']>=args_in.max_jobs):
      break
  if (args_in.max_jobs is not None) and \
     (job_status['submitted']>=args_in.max_jobs):
    break

format_str = "Jobs: {submitted} submitted, {skipped} skipped, " + \
  "{no_complex} without complex files, {no_dock6} without dock6 files, " + \
  "{no_poses} have no docked poses, " + \
  "{missing_file} missing other files, {onq} on the queue, " + \
  "{recently_redone} recently redone, {complete} complete"
print format_str.format(**job_status)
