import os

from AlGDock import path_tools
from AlGDock.BindingPMF import DEBUG
from AlGDock.BindingPMF import HMStime

import time
import numpy as np

import MMTK
import MMTK.Units
from MMTK.ParticleProperties import Configuration

import multiprocessing
from multiprocessing import Process

# In APBS, minimum ratio of PB grid length to maximum dimension of solute
LFILLRATIO = 4.0  # For the ligand
RFILLRATIO = 2.0  # For the receptor/complex

class Postprocessing:
  """Postprocesses configurations

  Attributes
  ----------
  args : AlGDock.simulation_arguments.SimulationArguments
    Simulation arguments
  log : AlGDock.logger.Logger
    Simulation log
  top : AlGDock.topology.Topology
    Topology of the ligand
  top_RL : AlGDock.topology.Topology
    Topology of the complex
  system : AlGDock.system.System
    Simulation system
  data : AlGDock.simulation_data.SimulationData
    Stores results from the simulation
  """
  def __init__(self, args, log, top, top_RL, system, data, save):
    """Initializes the class

    Parameters
    ----------
    args : AlGDock.simulation_arguments.SimulationArguments
      Simulation arguments
    log : AlGDock.logger.Logger
      Simulation log
    top : AlGDock.topology.Topology
      Topology of the ligand
    top_RL : AlGDock.topology.Topology
      Topology of the complex
    system : AlGDock.system.System
      Simulation system
    data : AlGDock.simulation_data.SimulationData
      Stores results from the simulation
    save : AlGDock.BindingPMF.save
      Saves the data
    """
    self.args = args
    self.log = log
    self.top = top
    self.top_RL = top_RL
    self.system = system
    self.data = data
    self.save = save

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

  def run(self,
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
    if hasattr(self, 'elsize'):
      return

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

    (focus_dims, focus_center, focus_spacing) = self.system.getGridParams()
    focus_dims = roundUpDime(focus_dims)

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
