#!/usr/bin/env python

# AlGDock Pipeline script to actually run BPMF calculations


import tarfile
from os.path import abspath, basename, dirname, exists, getsize, \
    isdir, isfile, join
import os
import numpy as np
import glob
import inspect
import argparse
import sys
sys.path.insert(0, '/anvil/projects/x-mcb150144/ellanguyen/Targets/Installers/AlGDock-0.0.1/AlGDock')
import IO

def nonzero(path):
    return (isfile(path) and getsize(path) > 0) or \
           (isdir(path))


# Look for AlGDock and qsub_command.py
dirs = {}
dirs['current'] = os.getcwd()
dirs['script'] = '/anvil/projects/x-mcb150144/ellanguyen/Targets/Installers/AlGDock-0.0.1/Pipeline'

execfile(join(dirs['script'], '_external_paths.py'))
command_paths = findPaths(['qsub_command'])
algdock_path = findPath(search_paths['algdock'])
ff_path = findPath(search_paths['gaff'])

# Parse arguments
parser = argparse.ArgumentParser()

parser.add_argument('--script', default=algdock_path, help='AlGDock script')
parser.add_argument('--saved_arguments', default='saved_arguments.py',
                    help='File containing default values of parameters' +
                    '(can be overwritten by flags)')
parser.add_argument('--include_ligand', nargs='+', default=None,
                    help='Only runs AlGDock for these ligands')
parser.add_argument('--exclude_ligand', nargs='+', default=None,
                    help='Does not run AlGDock for these ligands')
parser.add_argument('--new_instance_per_ligand', action='store_true',
                    default=False,
                    help='Runs run_AlGDock.py for each ligand')
parser.add_argument('--include_receptor', nargs='+', default=None,
                    help='Only runs AlGDock for these receptors')
parser.add_argument('--exclude_receptor', nargs='+', default=None,
                    help='Does not run AlGDock for these receptor')
parser.add_argument('--older_than', type=int, default=None,
                    help='Only runs redo_free_energies jobs if f_RL.pkl.gz ' +
                    'is older than OLDER_THAN hours')
parser.add_argument('--check_complete', action='store_true', default=False,
                    help='Checks whether f_RL{0}.pkl.gz, where {0} is pose_string, is complete')
parser.add_argument('--clear_locks', action='store_true', default=False,
                    help='clears locks')
parser.add_argument('--create_empty_files', action='store_true', default=False)
parser.add_argument('--jobname_prefix', default=None,
                    help='Prefix for job names')

# Arguments related to file locations
parser.add_argument(
    '--forcefield', default=ff_path, help='File for Generalized AMBER Force Field')
parser.add_argument('--ligand', default='../ligand/AlGDock_in/',
                    help='The directory/file tree to look for ligand files')
parser.add_argument('--library_requirement', default=None,
                    help='The ligand file name must contain the string LIBRARY_REQUIREMENT')
parser.add_argument('--receptor', default='../receptor/amber_in',
                    help='The directory/file to look for receptor files (prmtop and inpcrd)')
parser.add_argument('--receptor_grids', default='../receptor/AlGDock_in',
                    help='The directory to look for receptor grids (nc)')
parser.add_argument('--complex', default='../complex/AlGDock_in',
                    help='The directory tree to look for complex files (prmtop and inpcrd)')
parser.add_argument('--site_info',
                    default='../receptor/2-binding_site/measured_binding_site.py',
                    help='Python script with binding site parameters ' +
                    '(site_R, half_edge_length)')
parser.add_argument('--dock6', default='../dock6/',
                    help='The directory to look for dock6 results (*.mol2.gz)')
parser.add_argument('--tree_cool', default='BC/',
                    help='Directory tree to store cooling results')
parser.add_argument('--tree_dock', default='CD/',
                    help='Directory tree to store docking results')
# Arguments related to job management
parser.add_argument('--reps', default=None, nargs=2, type=int,
                    help='Range of repetitions')
parser.add_argument('--max_jobs', default=None, type=int)
parser.add_argument('--calcs_per_job', default=1, type=int,
                    help='Number of jobs to include in a single submission script.')
parser.add_argument('--first_ligand', default=None, type=int)
parser.add_argument('--max_ligands', default=None, type=int)
parser.add_argument('--dry', action='store_true', default=False,
                    help='Does not actually submit the job to the queue')
parser.add_argument('--no_release', action='store_true', default=False,
                    help='Does not release held jobs')
parser.add_argument('--interactive', action='store_true', default=False,
                    help='Output command for running in an interactive python environment')
parser.add_argument('--check_tarballs', action='store_true', default=False,
                    help='Check inside tarballs for files')
parser.add_argument('--skip_onq', action='store_true', default=False,
                    help='Skips looking for the job on the queue')
# Arguments related to scoring and assessment
parser.add_argument('--score', default=None,
                    help='File for starting pose or keyword: xtal (crystal structure), ' +
                    'dock (another docking program), xtal_plus_dock')
parser.add_argument('--poses', default=[-1, 0], nargs=2, type=int,
                    help='Confine the ligand to specific poses from the "score" file. ' +
                    'If the argument is "-1 0" there is no pose restriction. ' +
                    'Otherwise the argument is the range of poses.')
parser.add_argument('--rmsd', choices=['xtal', None], default=None,
                    help='Calulates the rmsd between snapshots a configuration. ' +
                    'xtal means the crystal structure.')

# Simulation settings and constants
#   Run-dependent
parser.add_argument('--BC_repX_cycles', type=int,
                    help='Number of replica exchange cycles for cooling')
parser.add_argument('--CD_repX_cycles', type=int,
                    help='Number of replica exchange cycles for docking')
parser.add_argument('--run_type',
                    choices=['configuration_energies',
                             'minimized_configuration_energies',
                             'store_params', 'initial_BC', 'BC',
                             'initial_CD', 'CD', 'postprocess',
                             'redo_postprocess', 'free_energies', 'redo_free_energies', 'all',
                             'timed', 'timed_BC', 'timed_CD',
                             'render_docked', 'render_intermediates',
                             'clear_intermediates', None],
                    help='Type of calculation to run')
parser.add_argument('--max_time', type=int, default=600,
                    help='For timed calculations, the maximum amount of wall clock time, ' +
                    'in minutes')
parser.add_argument('--keep_tar', action='store_true',
                    help='Keep files extracted from tar input')
parser.add_argument('--cores', type=int,
                    help='Number of CPU cores to use')
#   Defaults
parser.add_argument('--protocol', choices=['Adaptive', 'Geometric'],
                    help='Approach to determining series of thermodynamic states')
parser.add_argument('--therm_speed', type=float,
                    help='Thermodynamic speed during adaptive simulation')
parser.add_argument('--sampler',
                    choices=['HMC', 'MixedHMC', 'NUTS', 'VV'],
                    help='Sampling method')
parser.add_argument('--MCMC_moves', type=int,
                    help='Types of MCMC moves to use')
parser.add_argument('--T_HIGH', type=float,
                    help='High temperature')
parser.add_argument('--T_SIMMIN', type=float,
                    help='Lowest simulated temperature')
parser.add_argument('--T_TARGET', type=float,
                    help='Target temperature')
parser.add_argument('--temperature_scaling',
                    choices=['Linear', 'Quadratic'],
                    help='Determines whether the temperature changes linearly with the docking progress variable or quadratically with the grid scaling progress variable. (only for docking)')
parser.add_argument('--H_mass', type=float,
                    help='The repartitioned mass of hydrogen. Set negative to turn off HMR')
parser.add_argument('--delta_t', type=float,
                    help='The default time step, in fs')
# For initialization
parser.add_argument('--seeds_per_state', type=int,
                    help='Number of starting configurations in each state during initialization')
parser.add_argument('--steps_per_seed', type=int,
                    help='Number of MD steps per state during initialization')
parser.add_argument('--darts_per_seed', type=int,
                    help='Number of smart darting attempts for each seed during initialization')
# For replica exchange
parser.add_argument('--repX_cycles', type=int,
                    help='Number of replica exchange cycles for docking and cooling')
parser.add_argument('--sweeps_per_cycle', type=int,
                    help='Number of replica exchange sweeps per cycle')
parser.add_argument('--snaps_per_cycle', type=int,
                    help='Number of snapshots to save per cycle')
parser.add_argument('--attempts_per_sweep', type=int,
                    help='Number of replica exchange attempts per sweeps')
parser.add_argument('--steps_per_sweep', type=int,
                    help='Number of MD steps per replica exchange sweep')
parser.add_argument('--darts_per_sweep', type=int,
                    help='Number of smart darting attempts per replica exchange sweep')
parser.add_argument('--sampling_importance_resampling', action='store_true', default=None,
                    help='Use sampling importance resamping')
parser.add_argument('--solvation', default='Desolvated',
                    choices=['Desolvated', 'Reduced', 'Full', 'Fractional'],
                    help='How to use OBC implicit solvent during sampling.')
parser.add_argument('--keep_intermediate', action='store_true', default=None,
                    help='Keep configurations for intermediate states?')
parser.add_argument('--min_repX_acc', type=float,
                    help='Minimum value for replica exchange acceptance rate')
# For postprocessing
parser.add_argument('--phases', nargs='+',
                    help='Phases to use in postprocessing')
#   Stored in dir_cool
parser.add_argument('--BC_protocol', choices=['Adaptive', 'Set'],
                    help='Approach to determining series of thermodynamic states')
parser.add_argument('--BC_therm_speed', type=float,
                    help='Thermodynamic speed during adaptive simulation')
parser.add_argument('--BC_sampler',
                    choices=['HMC', 'NUTS', 'VV'],
                    help='Sampling method')
# For initialization
parser.add_argument('--BC_seeds_per_state', type=int,
                    help='Number of starting configurations in each state during initialization')
parser.add_argument('--BC_steps_per_seed', type=int,
                    help='Number of MD steps per state during initialization')
# For replica exchange
parser.add_argument('--BC_sweeps_per_cycle', type=int,
                    help='Number of replica exchange sweeps per cycle')
parser.add_argument('--BC_attempts_per_sweep', type=int,
                    help='Number of replica exchange attempts per sweeps')
parser.add_argument('--BC_steps_per_sweep', type=int,
                    help='Number of MD steps per replica exchange sweep')
parser.add_argument('--BC_keep_intermediate', action='store_true',
                    default=None, help='Keep configurations for intermediate states?')
#   Stored in dir_dock
parser.add_argument('--CD_protocol', choices=['Adaptive', 'Set'],
                    help='Approach to determining series of thermodynamic states')
parser.add_argument('--CD_therm_speed', type=float,
                    help='Thermodynamic speed during adaptive simulation')
parser.add_argument('--CD_sampler',
                    choices=['HMC', 'NUTS', 'VV'],
                    help='Sampling method')
# For initialization
parser.add_argument('--CD_seeds_per_state', type=int,
                    help='Number of starting configurations in each state during initialization')
parser.add_argument('--CD_attempts_per_sweep', type=int,
                    help='Number of replica exchange attempts per sweeps')
parser.add_argument('--CD_steps_per_seed', type=int,
                    help='Number of MD steps per state during initialization')
# For replica exchange
parser.add_argument('--CD_sweeps_per_cycle', type=int,
                    help='Number of replica exchange sweeps per cycle')
parser.add_argument('--CD_steps_per_sweep', type=int,
                    help='Number of MD steps per replica exchange sweep')
parser.add_argument('--CD_keep_intermediate', action='store_true',
                    default=None, help='Keep configurations for intermediate states?')

parser.add_argument('--site',
                    choices=['Sphere', 'Cylinder', 'Measure'],
                    help='Type of binding site. "Measure" means that parameters' +
                    ' for a sphere will be measured from docked configurations.')
parser.add_argument('--site_center', nargs='+',
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
parser.add_argument('--receptor_NAMD_OBC', type=float, nargs='+',
                    help='Receptor potential energies in NAMD GBSA implicit solvent (in units of kJ/mol)')
args_in = parser.parse_args()

# Convert argument paths to absolute paths
for path in ['ligand', 'receptor', 'receptor_grids', 'complex',
             'tree_dock', 'tree_cool', 'dock6']:
    setattr(args_in, path, os.path.abspath(getattr(args_in, path)))


# Determine which jobs are on the queue
if not args_in.skip_onq:
    execfile(join(dirs['script'], '_jobs_on_queue.py'))
    onq = jobs_on_queue()
else:
    onq = []

# Get other parameters
general_sim_arg_keys = ['protocol', 'therm_speed',
                        'sampler', 'seeds_per_state', 'steps_per_seed', 'repX_cycles',
                        'sweeps_per_cycle', 'snaps_per_cycle',
                        'attempts_per_sweep', 'steps_per_sweep', 'darts_per_sweep',
                        'sampling_importance_resampling', 'solvation',
                        'keep_intermediate', 'phases']
sim_arg_keys = general_sim_arg_keys + \
    ['BC_'+a for a in general_sim_arg_keys] + \
    ['CD_'+a for a in general_sim_arg_keys] + \
    ['MCMC_moves', 'score', 'rmsd', 'pose', 'run_type',
        'max_time', 'keep_tar', 'cores',
        'site', 'site_center', 'site_direction',
        'site_max_X', 'site_max_R', 'site_density',
        'T_HIGH', 'T_TARGET', 'T_SIMMIN', 'reps']

# Load saved arguments
if (args_in.saved_arguments is not None) and nonzero(args_in.saved_arguments):
    print 'Passed arguments:'
    print 'Loading default simulation arguments from '+args_in.saved_arguments
    execfile(args_in.saved_arguments)
    for key in sim_arg_keys:
        if hasattr(args_saved, key) and \
                ((not hasattr(args_in, key)) or (getattr(args_in, key) is None)):
            setattr(args_in, key, getattr(args_saved, key))
else:
    class args_saved:
        pass

# Report arguments in args_in that are not listed in sim_arg_keys
print "\nThese arguments will not be passed to AlGDock:"
for key in args_in.__dict__.keys():
    if key not in sim_arg_keys and getattr(args_in, key) is not None:
        print '{0}={1}'.format(key, getattr(args_in, key))
print

# Determine the binding site radius and half edge length
site = args_in.site
if hasattr(args_in, 'site_center') and hasattr(args_in, 'site_max_R'):
    site_center = args_in.site_center
    site_max_R = args_in.site_max_R
else:
    # If no argument, load in a file
    if isfile(args_in.site_info):
        execfile(args_in.site_info)
        # These should be in nanometers, not Angstroms
        site_center = [half_edge_length*0.1,
                       half_edge_length*0.1, half_edge_length*0.1]
        site_max_R = site_R/10.
    else:
        print 'No binding site information'

# Look for ligand files
if isfile(args_in.ligand):
    ligand_FNs = [abspath(args_in.ligand)]
elif isdir(args_in.ligand):
    ligand_FNs = glob.glob(join(args_in.ligand, '*/*.tar.gz'))
    ligand_FNs = sorted([abspath(FN) for FN in ligand_FNs])
else:
    raise Exception('Ligand input %s is not a file or directory!' %
                    args_in.ligand)
# Check that ligand tarball has nonzero size
ligand_FNs = [FN for FN in ligand_FNs if getsize(FN) > 0]
if len(ligand_FNs) > 1:
    # Filter ligands
    if args_in.include_ligand is not None:
        ligand_FNs = [FN for FN in ligand_FNs
                      if np.array([FN.find(ligN) != -1
                                   for ligN in args_in.include_ligand]).any()]
    if args_in.exclude_ligand is not None:
        ligand_FNs = [FN for FN in ligand_FNs
                      if not np.array([FN.find(ligN) != -1
                                       for ligN in args_in.exclude_ligand]).any()]

    if args_in.new_instance_per_ligand:
        import sys
        for ligand_FN in ligand_FNs:
            command = ' '.join(sys.argv + ['--ligand', ligand_FN]) + ' &'
            print command
            os.system(command)
        sys.exit()
# Sort ligands in reverse order by size. This way,
# slower jobs will be queued sooner and
# jobs that start around the same time will likely finish around the same time
ligand_FNs.sort(key=lambda x: getsize(x), reverse=True)

# Look for receptor files
if isfile(args_in.receptor):
    receptor_FNs = [args_in.receptor]
elif isdir(args_in.receptor):
    receptor_FNs = glob.glob(join(args_in.receptor, '*.prmtop'))
else:
    raise Exception('Receptor input %s is not a file or directory!' %
                    args_in.receptor)
receptor_FNs = [abspath(FN) for FN in receptor_FNs]
# Filter receptors
if args_in.include_receptor is not None:
    receptor_FNs = [FN for FN in receptor_FNs
                    if np.array([FN.find(recN) != -1
                                 for recN in args_in.include_receptor]).any()]
if args_in.exclude_receptor is not None:
    receptor_FNs = [FN for FN in receptor_FNs
                    if not np.array([FN.find(recN) != -1
                                     for recN in args_in.exclude_receptor]).any()]

# Require inpcrd as well as prmtop files
receptor_FNs = [FN for FN in receptor_FNs if
                np.array([nonzero(FN[:-6]+key)
                          for key in ['prmtop', 'inpcrd']]).all()]
# Require grid files #Testing use ELE from AlGDock not APBS
receptor_FNs = [FN for FN in receptor_FNs if np.array(
    [nonzero(join(args_in.receptor_grids, basename(FN)[:-6]+key))
     for key in ['LJa.25.nc', 'LJr.25.nc', 'PB.nc']]).all()] #Change from PB to ele

# Look for complex files
if isfile(args_in.complex):
    complex_FNs = [args_in.complex]
elif isdir(args_in.complex):
    complex_FNs = glob.glob(join(args_in.complex, '*/*/*.tar.gz'))
    if len(complex_FNs) == 0:
        complex_FNs = glob.glob(join(args_in.complex, '*/*.tar.gz'))
else:
    raise Exception('Complex input %s is not a file or directory!' %
                    args_in.complex)
# Require that complex tarball has nonzero size
complex_FNs = [FN for FN in complex_FNs if nonzero(FN)]

print 'Found %d ligands, %d receptors, and %d complexes ready for AlGDock' % (
    len(ligand_FNs), len(receptor_FNs), len(complex_FNs))

if (args_in.library_requirement is not None):
    ligand_FNs = [FN for FN in ligand_FNs
                  if FN.find(args_in.library_requirement) > -1]
    print '%d ligand(s) meet the library requirement' % (len(ligand_FNs))

if args_in.reps is None:
    args_in.reps = [0, 1]
if args_in.first_ligand is None:
    args_in.first_ligand = 0
if args_in.max_ligands is None:
    args_in.max_ligands = len(ligand_FNs)

print 'Arguments:'
print args_in

namespace = locals()

status = {'jobs': 0, 'submitted': 0, 'skipped': 0,
          'no_complex': 0, 'no_configurations': 0, 'no_BC': 0, 'no_dock6': 0,
          'missing_file': 0, 'onq': 0, 'recently_redone': 0, 'locked': 0, 'complete': 0}

checked = []
terminal_commands = []

# Loops are over ligands, repetitions, receptors, and poses

# Loop over ligands
for ligand_FN in ligand_FNs[args_in.first_ligand:args_in.first_ligand+args_in.max_ligands]:
    labels = {}
    if isfile(args_in.ligand):
        lib_and_key_prefix = ligand_FN.split('/')[-2]
    else:
        lib_and_key_prefix = dirname(ligand_FN[len(args_in.ligand)+1:])
    labels['library'] = '.'.join(lib_and_key_prefix.split('.')[:-1])
    labels['key'] = basename(ligand_FN[:-7])
    labels['lib_subdir'] = labels['library']+'.'+labels['key'][:-2]+'__'
    labels['ligand'] = labels['library']+'.'+labels['key']

    # Define and check files within the ligand tarball
    paths_in_tar = {'ligand_database': labels['ligand'].lower()+'.db'}
    for key in [('ligand_prmtop', 'prmtop'), ('ligand_inpcrd', 'inpcrd'),
                ('ligand_mol2', 'mol2'),
                ('frcmodList', 'frcmod')]:
        paths_in_tar[key[0]] = labels['ligand']+'.'+key[1]
    if args_in.check_tarballs and (ligand_FN not in checked):
        tarF = tarfile.open(ligand_FN)
        names = [m.name for m in tarF.getmembers()]
        not_found = [paths_in_tar[key] for key in paths_in_tar.keys()
                     if not paths_in_tar[key] in names]
        if len(not_found) > 0:
            print 'The following files were missing in '+ligand_FN+':'
            print ' '.join(not_found)
            continue
        else:
            checked.append(ligand_FN)

    #print('# Loop over repetitions')
    
    for rep in range(args_in.reps[0], args_in.reps[1]):
        paths = {'dir_BC': join(args_in.tree_cool,
                                labels['lib_subdir'], '%s-%d' % (labels['key'], rep)),
                 'forcefield': args_in.forcefield,
                 'ligand_tarball': ligand_FN}

        if not isdir(paths['dir_BC']):
            os.system('mkdir -p '+paths['dir_BC'])
        if nonzero(join(paths['dir_BC'], 'f_L.pkl.gz')):
            if (args_in.run_type in ['initial_BC', 'BC', 'timed_BC']):
                status['complete'] += 1
                continue  # Cooling is already done
        else:
            if (args_in.run_type in ['CD', 'timed_CD']):
                status['no_BC'] += 1
                continue  # Not ready for docking because cooling isn't complete

        
        # Loop over receptor files
        for receptor_FN in receptor_FNs:
            labels['receptor'] = basename(receptor_FN)[:-7]
            labels['complex'] = labels['library'] + \
                '.'+labels['key']+'-'+labels['receptor']
            # Identify receptor files
            for key in ['prmtop', 'inpcrd']:
                paths['receptor_'+key] = abspath(receptor_FN[:-6]+key)
            if isfile(abspath(receptor_FN[:-6]+'pdb')):
                paths['receptor_fixed_atoms'] = abspath(receptor_FN[:-6]+'pdb')
            for key in [('grid_LJa', 'LJa.25.nc'), ('grid_LJr', 'LJr.25.nc'),
                        ('grid_ELE', 'PB.nc')]: 
                paths[key[0]] = join(args_in.receptor_grids,
                                     '%s.%s' % (labels['receptor'], key[1]))

            # Identify complex tarball
            complex_tar_FN = join(args_in.complex,
                                  labels['lib_subdir'], labels['key'], labels['receptor']+'.tar.gz')
            if not (isfile(complex_tar_FN)):
                complex_tar_FN = join(args_in.complex, labels['lib_subdir'],
                                      '{0}.{1}-{2}.tar.gz'.format(
                    labels['library'], labels['key'], labels['receptor']))
                if not (isfile(complex_tar_FN)):
                    print 'No complex tarfile ' + complex_tar_FN
                    status['no_complex'] += 1
                    continue  # Complex files are missing
            paths['complex_tarball'] = complex_tar_FN

           # Identify docking directory
            paths['dir_CD'] = join(args_in.tree_dock,
                                   labels['lib_subdir'], labels['key'], '%s-%d' % (labels['receptor'], rep))
            if not isdir(paths['dir_CD']):
                os.system('mkdir -p '+paths['dir_CD'])

            # Convert relative path to absolute paths
            for key in paths.keys():
                #print(key, paths[key])
                paths[key] = abspath(paths[key])

            # Define and check files within the complex tarball
            for key in ['prmtop', 'inpcrd']:
                paths_in_tar['complex_'+key] = labels['complex']+'.'+key
            if 'receptor_fixed_atoms' in paths.keys():
                paths_in_tar['complex_fixed_atoms'] = labels['complex']+'.pdb'
            if args_in.check_tarballs and (complex_tar_FN not in checked):
                tarF = tarfile.open(complex_tar_FN)
                names = [m.name for m in tarF.getmembers()]
                not_found = [paths_in_tar[key] for key in paths_in_tar.keys()
                             if key.startswith('complex') and not paths_in_tar[key] in names]
                if len(not_found) > 0:
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
                status['missing_file'] += 1
                continue  # Files are missing

            # Determine the score argument
            if args_in.score is None:
                score = False
            elif args_in.score == 'xtal':
                score = 'default'
            elif args_in.score == 'dock':
                score = abspath(join(args_in.dock6,
                                          labels['lib_subdir'], labels['key'],
                                          labels['receptor'] + '.nc'))
                if not nonzero(score):
                    if isfile(score[:-3]+'.mol2.gz'):
                        status['no_configurations'] += 1
                        skip_job = True
                        break  # No configurations in dock6
                    else:
                        print 'No dock6 output in '+args_in.score
                        status['no_dock6'] += 1
                        skip_job = True
                        break  # Dock6 files are missing
            elif args_in.score == 'xtal_plus_dock':
                score = abspath(join(dir_dock6,
                                     'xtal_plus_dock6_scored.mol2'))
            else:
                score = args_in.score

            # Check for completed configuration energies
            if args_in.run_type.endswith('configuration_energies') and \
                    (not args_in.interactive):
                prefix = 'xtal' if score == 'default' else \
                    os.path.basename(score).split('.')[0]
                if args_in.run_type == 'minimized_configuration_energies':
                    prefix = 'min_' + prefix
                FN = join(paths['dir_CD'], prefix + '.pkl.gz')
                complete = False
                if os.path.isfile(FN):
                    import gzip
                    import pickle
                    F = gzip.open(FN, 'r')
                    (confs, Es) = pickle.load(F)
                    F.close()
                    complete = True
                    for phase in args_in.phases:
                        for moeity in ['R', 'L', 'RL']:
                            if moeity+phase not in Es.keys():
                                print '{0} is missing {1}'.format(score, moeity+phase)
                                complete = False
                # Pose energy is done
                if complete:
                    status['complete'] += 1
                    continue

            # Find out the number of poses
            pose_upper_range = args_in.poses[1]
            if isfile(complex_tar_FN):
                tarF = tarfile.open(complex_tar_FN)
                try:
                    poses_in_mdcrd = 'poses.mdcrd' in tarF.getnames()
                except IOError:
                    print 'Error with '+complex_tar_FN
                    import sys
                    sys.exit()
                if poses_in_mdcrd:
                    tarF.extract('poses.mdcrd')
                    IO_crd = IO.crd()
                    (crd, title) = IO_crd.read('poses.mdcrd', return_title=True)
                    os.remove('poses.mdcrd')
                    pose_upper_range = min(
                        int(title.split()[0]), args_in.poses[1])

            # Loop over poses
            for pose in range(args_in.poses[0], pose_upper_range):
                if pose == -1:
                    pose_string = ''
                else:
                    pose_string = '_pose%03d' % pose

                f_RL_FN = join(paths['dir_CD'],
                               'f_RL{0}.pkl.gz'.format(pose_string))
                if nonzero(f_RL_FN):
                    import time
                    # Check if the calculation was recently redone
                    if (args_in.run_type == 'redo_free_energies') and \
                        (args_in.older_than is not None) and \
                            (time.time()-os.path.getmtime(f_RL_FN))/60./60. < args_in.older_than:
                        status['recently_redone'] += 1
                        continue

                    if (args_in.run_type in
                            ['random_CD', 'initial_CD', 'CD', 'all', 'timed']):
                        print('GO with CD')
                        if args_in.check_complete:
                            import gzip
                            import pickle
                            F = gzip.open(f_RL_FN, 'r')
                            dat = pickle.load(F)
                            F.close()
                            try:
                                completed_cycles = np.min(
                                    [len(dat[-1][p+'_MBAR']) for p in args_in.phases])
                                complete = (completed_cycles >= int(
                                    args_in.CD_repX_cycles))
                            except:
                                print 'Error in '+f_RL_FN
                                completed_cycles = 0
                                complete = False
                            if complete:
                                status['complete'] += 1  # Docking is done
                                continue
                            else:
                                print '%d/%d cycles in %s' % (
                                    completed_cycles, args_in.CD_repX_cycles, paths['dir_CD'])
                        else:
                            status['complete'] += 1  # Docking is done
                            continue

                elif (args_in.run_type == 'redo_free_energies'):
                    status['missing_file'] += 1
                    continue
                
                jobname = '' if args_in.jobname_prefix is None \
                    else (args_in.jobname_prefix + '.')
                jobname += '{0}-{1}{2}.{3}'.format(labels['complex'],
                                                   rep, pose_string, args_in.run_type)
                if jobname in onq:
                    status['onq'] += 1
                    print jobname + ' is on the queue'
                    continue  # Job is on the queue

                if args_in.clear_locks:
                    if not (jobname in onq):
                        if exists(join(paths['dir_BC'], '.lock')):
                            # There is a lock in the cooling directory
                            print '# Removing lock in cooling directory %s' % (paths['dir_BC'])
                            os.remove(join(paths['dir_BC'], '.lock'))
                        if exists(join(paths['dir_CD'], '.lock'+pose_string)):
                            # There is a lock in the docking directory
                            print '# Removing lock in docking directory %s' % (paths['dir_CD'])
                            os.remove(
                                join(paths['dir_CD'], '.lock'+pose_string))

                # Consolidate paths to pass
                paths_to_pass = {}
                for key in paths.keys():
                    paths_to_pass[key] = paths[key]
                for key in paths_in_tar.keys():
                    paths_to_pass[key] = paths_in_tar[key]

                interactive_to_pass = []
                terminal_to_pass = []
                skip_job = False
                for key in sorted(paths_to_pass.keys() + sim_arg_keys):
                    # Priority is passed arguments (which may include saved arguments),
                    #   local variables,
                    #   and then the path dictionary
                    if hasattr(args_in, key) and (getattr(args_in, key) is not None):
                        val = getattr(args_in, key)
                    elif (key in namespace.keys()) and (namespace[key] is not None):
                        val = namespace[key]
                    elif key in paths_to_pass.keys():
                        val = paths_to_pass[key]
                    else:
                        continue
                    # Special cases
                    if key == 'rmsd':
                        val = {None: False, 'xtal': True}[val]
                    elif key == 'score':
                        val = score
                    elif key == 'frcmodList':
                        val = [val]
                    # Actual strings to pass
                    if isinstance(val, str):
                        interactive_to_pass.append("%s = '%s'" % (key, val))
                        terminal_to_pass.append("--%s %s" % (key, val))
                    elif isinstance(val, bool):
                        if val:
                            interactive_to_pass.append("%s = %s" % (key, val))
                            terminal_to_pass.append("--%s" % (key))
                    elif isinstance(val, int):
                        interactive_to_pass.append("%s = %d" % (key, val))
                        terminal_to_pass.append("--%s %d" % (key, val))
                    elif isinstance(val, float):
                        interactive_to_pass.append("%s = %.5f" % (key, val))
                        terminal_to_pass.append("--%s %.5f" % (key, val))
                    elif isinstance(val, list):
                        if isinstance(val[0], str):
                            interactive_to_pass.append("%s = ['%s']" % (key,
                                                                        "', '".join([a for a in val])))
                            terminal_to_pass.append("--%s %s" % (key,
                                                                 " ".join([a for a in val])))
                        elif isinstance(val[0], float):
                            interactive_to_pass.append("%s = [%s]" % (key,
                                                                      ", ".join(['%.5f' % a for a in val])))
                            terminal_to_pass.append("--%s %s" % (key,
                                                                 " ".join(['%.5f' % a for a in val])))
                    else:
                        print 'Value:', val
                        raise Exception('Type not known!')

                if skip_job:
                    continue

                outputFNs = {}
                for FN in ['cool_log.txt',
                           'cool_progress.pkl.gz', 'cool_progress.pkl.gz.BAK',
                           'cool_data.pkl.gz', 'cool_data.pkl.gz.BAK',
                           'f_L.pkl.gz']:
                    outputFNs[FN] = join(paths_to_pass['dir_BC'], FN)
                for FN in ['dock%s_log.txt' % pose_string,
                           'dock_progress%s.pkl.gz' % pose_string,
                           'dock_progress%s.pkl.gz.BAK' % pose_string,
                           'dock_data%s.pkl.gz' % pose_string,
                           'dock_data%s.pkl.gz.BAK' % pose_string,
                           'f_RL%s.pkl.gz' % pose_string]:
                    outputFNs[FN] = join(paths_to_pass['dir_CD'], FN)
                for k in outputFNs.keys():
                    if args_in.create_empty_files and (not isfile(outputFNs[k])):
                        open(outputFNs[k], 'a').close()
                transfer_output_remaps = [f(key)
                                          for key in outputFNs.keys()
                                          for f in (lambda x:key, lambda x:outputFNs[key])]

                interactive_command = "self = AlGDock.BindingPMF.BPMF(" + \
                    ', \\\n  '.join(sorted(interactive_to_pass)) + \
                    ")"
                if args_in.script is None:
                    terminal_command = '$ALGDOCK '
                else:
                    if args_in.script.find('.py') != -1:
                        terminal_command = 'python ' + args_in.script + ' '
                    else:
                        terminal_command = args_in.script + ' '
                terminal_command += ' \\\n  '.join(terminal_to_pass)

                if args_in.interactive:
                    print interactive_command
                else:
                    import subprocess
                    if args_in.calcs_per_job > 1:
                        terminal_commands.append(terminal_command)
                        if len(terminal_commands) == args_in.calcs_per_job:
                            command_string = ' &; '.join([
                                command.replace(' \\\n', '').replace('  ', ' ')
                                for command in terminal_commands]) + ' &;'
                            subprocess.call(['python', command_paths['qsub_command'],
                                             jobname, command_string] +
                                            {True: ['--dry'], False: []}[args_in.dry])
                            status['jobs'] += 1
                            terminal_commands = []
                    else:
                        subprocess.call(['python', command_paths['qsub_command'],
                                         jobname, terminal_command, '--mem', '20'] +
                                        ['--input_files'] + outputFNs.values() +
                                        ['--output_files'] + list(outputFNs) +
                                        ['--output_remaps'] + transfer_output_remaps +
                                        ['--comment', interactive_command.replace(' \\\n', '')] +
                                        {True: ['--dry'], False: []}[args_in.dry] +
                                        {True: ['--no_release'], False: []}[args_in.no_release])
                        status['jobs'] += 1

                status['submitted'] += 1

                if (args_in.max_jobs is not None) and \
                   (status['jobs'] >= args_in.max_jobs):
                    break
            if (args_in.max_jobs is not None) and \
               (status['jobs'] >= args_in.max_jobs):
                break
        if (args_in.max_jobs is not None) and \
           (status['jobs'] >= args_in.max_jobs):
            break
    if (args_in.max_jobs is not None) and \
       (status['jobs'] >= args_in.max_jobs):
        break

# Submit the remaining calculations in a final job
if (args_in.calcs_per_job > 1) and len(terminal_commands) > 0 \
        and (not args_in.interactive):
    import subprocess
    command_string = ' &; '.join([
        command.replace(' \\\n', '').replace('  ', ' ')
        for command in terminal_commands]) + ' &;'
    subprocess.call(['python', command_paths['qsub_command'],
                     jobname, command_string] +
                    {True: ['--dry'], False: []}[args_in.dry])
    status['jobs'] += 1
    terminal_commands = []

format_str = "Calculation status: {jobs} jobs, {submitted} submitted, " + \
    "{skipped} skipped, " + \
    "{no_complex} without complex files, {no_BC} without cooling files" + \
    "{no_dock6} without dock6 files, " + \
    "{no_configurations} have no docked configurations, " + \
    "{missing_file} missing other files, {onq} on the queue, " + \
    "{recently_redone} recently redone, {complete} complete"
print format_str.format(**status)
