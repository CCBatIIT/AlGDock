#!/usr/bin/env python

import argparse
parser = argparse.ArgumentParser()

parser.add_argument('--saved_arguments', default='saved_arguments.py',
  help='File containing default values of parameters' + \
       '(can be overwritten by flags)')

# Arguments related to file locations
parser.add_argument('--ligand', default='../ligand/AlGDock_in/', \
  help='The directory/file to look for ligand files' + \
       ' (prmtop, inpcrd, frcmod, and db)')
parser.add_argument('--receptor', default='../receptor/amber_in', \
  help='The directory/file to look for receptor files (prmtop and inpcrd)')
parser.add_argument('--receptor_grids', default='../receptor/AlGDock_in', \
  help='The directory to look for receptor grids (nc)')
parser.add_argument('--complex', default='../complex/AlGDock_in', \
  help='The directory to look for complex files (prmtop and inpcrd)')
parser.add_argument('--site_info', \
  default='../receptor/2-binding_site/measured_binding_site.py', \
  help='Python script with binding site parameters' + \
       ' (site_R, half_edge_length)')
parser.add_argument('--dock6', default='../dock6/', \
  help='The directory to look for dock6 results (*.mol2.gz)')
parser.add_argument('--tree_cool', default='cool/',
  help='Directory tree to store cooling results')
parser.add_argument('--tree_dock', default='dock/',
  help='Directory tree to store docking results')
# Arguments related to job management
parser.add_argument('--reps', default=[0,1], nargs=2, type=int, \
  help='Range of repetitions')
parser.add_argument('--max_jobs', default=None, type=int)
parser.add_argument('--max_ligands', default=None, type=int)
parser.add_argument('--dry', action='store_true', default=False, \
  help='Does not actually submit the job to the queue')
parser.add_argument('--no_release', action='store_true', default=False, \
  help='Does not release held jobs')
parser.add_argument('--interactive', action='store_true', default=False, help='Output command for running in an interactive python environment')
# Arguments related to scoring and assessment
parser.add_argument('--score', choices=['xtal','dock',None], default=None,
  help='Does a scoring rather than ab initio docking. xtal means the crystal structure and dock means from another docking program.')
parser.add_argument('--rmsd', choices=['xtal',None], default=None,
  help='Calulates the rmsd between snapshots a configuration. xtal means the crystal structure.')

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
parser.add_argument('--cores', type=int, \
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

#  parser.add_argument('--site',
#    choices=['Sphere','Cylinder'], help='Type of binding site')
#  parser.add_argument('--site_center', nargs=3, type=float,
#    help='Position of binding site center')
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
execfile(os.path.join(dirs['script'],'_jobs_on_queue.py'))
onq = jobs_on_queue()
command_paths = findPaths(['qsub_command','gaff.dat'])
algdock_path = findPath(search_paths['algdock'])

import numpy as np
def nonzero(path):
  return (os.path.isfile(path) and os.path.getsize(path)>0) or \
         (os.path.isdir(path))

# Look for ligand files
if os.path.isfile(args_in.ligand):
  ligand_FNs = [args_in.ligand]
elif os.path.isdir(args_in.ligand):
  ligand_FNs = glob.glob(os.path.join(args_in.ligand,'*.prmtop'))
else:
  raise Exception('Ligand input %s is not a file or directory!'%args_in.ligand)
# Require inpcrd, frcmod, and db files
ligand_FNs = [FN for FN in ligand_FNs if
  (np.array([nonzero(FN[:-6]+key) \
    for key in ['prmtop','inpcrd','frcmod']]).all() and
   nonzero(os.path.join(os.path.dirname(FN),
                 os.path.basename(FN)[:-6].lower()+'db')))]

# Look for receptor files
if os.path.isfile(args_in.receptor):
  receptor_FNs = [args_in.receptor]
elif os.path.isdir(args_in.receptor):
  receptor_FNs = glob.glob(os.path.join(args_in.receptor,'*.prmtop'))
else:
  raise Exception('Receptor input %s is not a file or directory!'%args_in.receptor)
# Require inpcrd as well as prmtop files
receptor_FNs = [FN for FN in receptor_FNs if
  np.array([nonzero(FN[:-6]+key) \
    for key in ['prmtop','inpcrd']]).all()]
# Require inpcrd as well as prmtop files
receptor_FNs = [FN for FN in receptor_FNs if np.array(\
  [nonzero(os.path.join(args_in.receptor_grids,os.path.basename(FN)[:-6]+key)) \
    for key in ['LJa.25.nc','LJr.25.nc','PB.nc']]).all()]

if os.path.isfile(args_in.complex):
  complex_FNs = [args_in.complex]
elif os.path.isdir(args_in.complex):
  complex_FNs = glob.glob(os.path.join(args_in.complex,'*.prmtop.gz'))
else:
  raise Exception('Complex input %s is not a file or directory!'%args_in.complex)
# Require gzipped inpcrd and frcmod files
complex_FNs = [FN for FN in complex_FNs if
  np.array([nonzero(FN[:-9]+key) \
    for key in ['prmtop.gz','inpcrd.gz']]).all()]
complex_labels = [os.path.basename(FN)[:-10] for FN in complex_FNs]

print 'Found %d ligands, %d receptors, and %d complexes ready for AlGDock'%(\
  len(ligand_FNs),len(receptor_FNs),len(complex_FNs))

# Load the binding site radius and half edge length
site = 'Sphere'
if os.path.isfile(args_in.site_info):
  execfile(args_in.site_info)
else:
  raise Exception('No binding site information')
# These should be in nanometers, not Angstroms
site_center = [half_edge_length*0.1, half_edge_length*0.1, half_edge_length*0.1]
site_max_R = site_R/10.

# Get other parameters
general_sim_arg_keys = ['protocol','no_protocol_refinement','therm_speed',\
  'sampler', 'seeds_per_state', 'steps_per_seed', 'repX_cycles', \
  'sweeps_per_cycle', 'steps_per_sweep', 'keep_intermediate', 'phases']
sim_arg_keys = general_sim_arg_keys + \
  ['cool_'+a for a in general_sim_arg_keys] + \
  ['dock_'+a for a in general_sim_arg_keys] + \
  ['MCMC_moves', 'score', 'rmsd', 'run_type', 'cores'] + \
  ['site', 'site_center', 'site_direction', \
   'site_max_X', 'site_max_R', 'site_density']

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

namespace = locals()

job_status = {'submitted':0, 'no_complex':0, 'no_dock6':0, 'missing_file':0, \
              'onq':0, 'complete':0}
for ligand_FN in ligand_FNs:
  labels = {'ligand':os.path.basename(ligand_FN)[:-7]}
  dir_cool = os.path.join(args_in.tree_cool,labels['ligand'])
  paths = [('forcefield',command_paths['gaff.dat'])]
  paths += [(key[0], ligand_FN[:-6]+key[1]) \
    for key in [('ligand_prmtop','prmtop'),('ligand_inpcrd','inpcrd'), \
                ('frcmodList','frcmod')]]
  paths += [('ligand_database', \
    os.path.join(os.path.dirname(ligand_FN),
                 os.path.basename(ligand_FN)[:-6].lower()+'db'))]
  if not os.path.isdir(dir_cool):
    os.system('mkdir -p '+dir_cool)
  paths += [('dir_cool', os.path.join(args_in.tree_cool,labels['ligand']))]
  if (args_in.run_type in ['initial_cool','cool']) and \
      nonzero(join(dir_cool,'f_L.pkl.gz')):
    job_status['complete'] += 1
    continue # Cooling is already done
  for receptor_FN in receptor_FNs:
    labels['receptor'] = os.path.basename(receptor_FN)[:-7]
    labels['complex'] = labels['ligand']+'-'+labels['receptor']
    if not (labels['complex'] in complex_labels):
      print 'No prmtop and inpcrd for '+labels['complex']
      job_status['no_complex'] += 1
      continue # Complex files are missing
    paths += [('receptor_'+key, os.path.abspath(receptor_FN[:-6]+key)) \
      for key in ['prmtop','inpcrd']]
    paths += \
      [('complex_'+key[0],
       os.path.join(args_in.complex,'%s.%s'%(labels['complex'],key[1]))) \
        for key in [('prmtop','prmtop.gz'),('inpcrd','inpcrd.gz')]]
    paths += \
      [(key[0], \
       os.path.join(args_in.receptor_grids, '%s.%s'%(labels['receptor'],key[1]))) \
        for key in [('grid_LJa','LJa.25.nc'),('grid_LJr','LJr.25.nc'),\
                    ('grid_ELE','PB.nc')]]

    if not np.array([nonzero(FN[1]) for FN in paths]).all():
      print 'Files are missing'
      job_status['missing_file']
      continue # Files are missing

    for rep in range(args_in.reps[0],args_in.reps[1]):
      labels['job'] = '%s-%d'%(labels['complex'],rep)
      dir_dock = os.path.join(args_in.tree_dock,labels['job'])
      paths += [('dir_dock',dir_dock)]
      if not os.path.isdir(dir_dock):
        os.system('mkdir -p '+dir_dock)
      if (args_in.run_type in ['random_dock','initial_dock', \
                               'dock','all','one_step']) and \
          nonzero(os.path.join(dir_dock,'f_RL.pkl.gz')):
        job_status['complete'] += 1
        continue # Docking is done
      jobname = '-'.join(dirs['current'].split('/')[-2:]+[labels['job']])
      if jobname in onq:
        job_status['onq'] += 1
        continue # Job is on the queue

      path_keys = [p[0] for p in paths]
      paths = dict(paths)
      for key in path_keys:
        paths[key] = os.path.abspath(paths[key])

      interactive_to_pass = []
      terminal_to_pass = []
      passError = False
      for key in (path_keys + sim_arg_keys):
        # Priority is passed arguments (which may include saved arguments),
        #   local variables,
        #   and then the path dictionary
        if hasattr(args_in,key) and (getattr(args_in,key) is not None):
          val = getattr(args_in,key)
        elif (key in namespace.keys()) and (namespace[key] is not None):
          val = namespace[key]
        elif key in path_keys:
          val = paths[key]
        else:
          continue
        # Special cases
        if key=='rmsd':
          val = {None: False, 'xtal':True}[val]
        elif key=='score':
          if val is None:
            val = False
          elif val=='xtal':
            val = True
          elif val=='dock':
            val = os.path.abspath(os.path.join(args_in.dock6, \
                               'ancg-'+labels['complex']+'.mol2.gz'))
            if (not nonzero(val)):
              print 'No dock6 output in '+val
              job_status['no_dock6'] += 1
              passError = True
              continue # Dock6 files are missing
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
              "' '".join([a for a in val])))
          elif isinstance(val[0],float):
            interactive_to_pass.append("%s = [%s]"%(key,
              ", ".join(['%.5f'%a for a in val])))
            terminal_to_pass.append("--%s %s"%(key,
              " ".join(['%.5f'%a for a in val])))
        else:
          raise Exception('Type not known!')

      if passError:
        continue
      
      outputFNs = {}
      for FN in ['cool_log.txt',
          'cool_params.pkl.gz','cool_params.pkl.gz',
          'cool_progress.pkl.gz','cool_progress.pkl.gz.BAK',
          'cool_data.pkl.gz','cool_data.pkl.gz.BAK',
          'f_L.pkl.gz']:
        outputFNs[FN] = os.path.join(dir_cool,FN)
      for FN in ['dock_log.txt',
          'dock_params.pkl.gz', 'dock_params.pkl.gz.BAK',
          'dock_progress.pkl.gz', 'dock_progress.pkl.gz.BAK',
          'dock_data.pkl.gz', 'dock_data.pkl.gz.BAK',
          'f_RL.pkl.gz']:
        outputFNs[FN] = os.path.join(dir_dock,FN)
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
        import subprocess
        subprocess.call(['python', command_paths['qsub_command'], \
          jobname, terminal_command, '--mem', '%d'%mem] + \
          ['--input_files'] + outputFNs.values() + \
          ['--output_files'] + list(outputFNs) + \
          ['--output_remaps'] + transfer_output_remaps + \
          ['--comment', interactive_command.replace(' \\\n','')] + \
          {True:['--dry'],False:[]}[args_in.dry] + \
          {True:['--no_release'],False:[]}[args_in.no_release])
      
      job_status['submitted'] += 1

      num_ligands = sum(job_status.values())
      if (args_in.max_jobs is not None) and \
         (job_status['submitted']>=args_in.max_jobs):
        break
      if (args_in.max_ligands is not None) and \
         (num_ligands>=args_in.max_ligands):
        break

    num_ligands = sum(job_status.values())
    if (args_in.max_jobs is not None) and \
       (job_status['submitted']>=args_in.max_jobs):
      break
    if (args_in.max_ligands is not None) and \
       (num_ligands>=args_in.max_ligands):
      break
      
  num_ligands = sum(job_status.values())
  if (args_in.max_jobs is not None) and \
     (job_status['submitted']>=args_in.max_jobs):
    break
  if (args_in.max_ligands is not None) and \
     (num_ligands>=args_in.max_ligands):
    break

print "Jobs: {submitted} submitted, {no_complex} without complex files, {no_dock6} without dock6 files, {missing_file} missing other files, {onq} on the queue, {complete} complete".format(**job_status)