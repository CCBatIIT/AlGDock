# This module implements a Smart Darting "integrator"
import numpy as np
# from Cython.Debugger.libpython import get_inferior_unicode_postfix

try:
  from MMTK import Configuration, Dynamics, Environment, Features, Trajectory, Units
  import MMTK_dynamics
except ImportError:
  MMTK = None

try:
  import openmm
  import openmm.unit as unit
  from openmm.app import AmberPrmtopFile, AmberInpcrdFile, Simulation, NoCutoff
  from openmm import *
except ImportError:
  OpenMM = None

R = 8.3144621*Units.J/Units.mol/Units.K
twoPi = 2*np.pi

#
# Smart Darting integrator
#
class SmartDartingIntegrator(Dynamics.Integrator):
  def __init__(self, universe, molecule, extended, confs=None, **options):
    """
    confs - configurations to dart to
    extended - whether or not to use external coordinates
    """
    Dynamics.Integrator.__init__(self, universe, options)
    # Supported features: none for the moment, to keep it simple
    self.features = []

    self.molecule = molecule
    self.extended = extended
    
    # Converter between Cartesian and BAT coordinates
    #TODO: MODIFY RIDIGBODY
    import AlGDock.rigid_bodies
    id = AlGDock.rigid_bodies.identifier(self.universe, self.molecule)
    #e.g.  id, <AlGDock.RigidBodies.identifier instance at 0x7faed3cce550>,
    #e.g.  id.initial_atom, Atom ligand.db.H4i3
    from AlGDock.BAT import converter
    self._BAT_util = converter(self.universe, self.molecule, \
      initial_atom = id.initial_atom)
    # BAT coordinates to perturb: external coordinates and primary torsions
    dof = self.universe.configuration().array.shape[0]*3 if extended \
      else self.universe.configuration().array.shape[0]*3 - 6
    self._BAT_to_perturb = range(6) if extended else []
    self._BAT_to_perturb += self._BAT_util.getFirstTorsionInds(extended)
    print('DEBUGGING !!!! _BAT_to_perturb', self._BAT_to_perturb)

    if confs is None:
      self.confs = None
    else:
      self.set_confs(confs)

  def set_confs(self, confs, rmsd_threshold=0.05, period_frac_threshold=0.35, \
      append=False):

    import time
    start_time = time.time()

    nconfs_attempted = len(confs)
    if append and (self.confs is not None):
      nconfs_o = len(self.confs)
      confs = confs + self.confs
    else:
      nconfs_o = 0

    # Minimize configurations
    if MMTK:
      from MMTK.Minimization import SteepestDescentMinimizer # @UnresolvedImport
      minimizer = SteepestDescentMinimizer(self.universe)

      minimized_confs = []
      minimized_energies = []
      for conf in confs:
        self.universe.setConfiguration(Configuration(self.universe, conf))
        x_o = np.copy(self.universe.configuration().array)
        e_o = self.universe.energy()
        for rep in range(50):
          minimizer(steps = 25)
          x_n = np.copy(self.universe.configuration().array)
          e_n = self.universe.energy()
          diff = abs(e_o-e_n)
          if np.isnan(e_n) or diff<0.05 or diff>1000.:
            self.universe.setConfiguration(Configuration(self.universe, x_o))
            break
          else:
            x_o = x_n
            e_o = e_n
        if not np.isnan(e_o):
          minimized_confs.append(x_o)
          minimized_energies.append(e_o)
    else:
      minimized_confs = []
      minimized_energies = []
      for conf in confs:
        #universe = top.OMM_simulation
        self.universe.context.setPositions(conf)
        state = self.universe.context.getState(getEnergy=True, getPositions=True)
        current_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        max_iterations = 1000
        tolerance = 0.01 * unit.kilojoule/unit.mole
        for _ in range(max_iterations):
          openmm.LocalEnergyMinimizer.minimize(self.universe.context, tolerance, max_iterations)
          state = self.universe.context.getState(getEnergy=True, getPositions=True)
          new_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
          energy_diff = abs(current_energy - new_energy)

          if np.isnan(new_energy) or energy_diff < 0.05 or energy_diff > 1000.0:
            break
          else:
            current_energy = new_energy

        if not np.isnan(current_energy):
          minimized_confs.append(state.getPositions(asNumpy=True))
          minimized_energies.append(current_energy)

    confs = minimized_confs
    energies = minimized_energies
    
    # Sort by increasing energy
    energies, confs = (list(l) \
      for l in zip(*sorted(zip(energies, confs), key=lambda p:p[0])))

    # Only keep configurations with energy with 12 kJ/mol of the lowest energy
    confs = [confs[i] for i in range(len(confs)) \
      if (energies[i]-energies[0])<12.]
    energies = energies[:len(confs)]

    if self.extended:
      # Keep only unique configurations, using rmsd as a threshold
      inds_to_keep = [0]
      for j in range(1,len(confs)):
        min_rmsd = np.min([np.sqrt(np.sum(np.square(\
          confs[j][self.molecule.heavy_atoms,:] - \
          confs[k][self.molecule.heavy_atoms,:]))/self.molecule.nhatoms) \
          for k in inds_to_keep])
        if min_rmsd>rmsd_threshold:
          inds_to_keep.append(j)
      confs = [confs[i] for i in inds_to_keep]
      energies = [energies[i] for i in inds_to_keep]

    confs_BAT = [self._BAT_util.BAT(confs[n], extended=self.extended) \
      for n in range(len(confs))]
    confs_BAT_tp = [confs_BAT[c][self._BAT_to_perturb] \
      for c in range(len(confs_BAT))]

    if not self.extended:
      # Keep only unique configurations based on period fraction threshold
      inds_to_keep = [0]
      for j in range(1,len(confs)):
        period_fracs = np.array([np.abs(confs_BAT_tp[j] - confs_BAT_tp[k]) \
          for k in inds_to_keep])/twoPi
        period_fracs = np.min([period_fracs,1-period_fracs],0)
        if (np.max(period_fracs,1)>period_frac_threshold).all():
          inds_to_keep.append(j)
      confs = [confs[i] for i in inds_to_keep]
      confs_BAT = [confs_BAT[i] for i in inds_to_keep]
      confs_BAT_tp = [confs_BAT_tp[i] for i in inds_to_keep]
      energies = [energies[i] for i in inds_to_keep]
    if MMTK:
      confs_ha = [confs[c][self.molecule.heavy_atoms,:] \
      for c in range(len(confs))]
    else:
      heavy_atoms = []
      for atom in self.molecule.atoms():
        if atom.element is not None and atom.element.atomic_number > 1:
          heavy_atoms.append(atom.index)
      confs_ha = [confs[c][heavy_atoms,:] \
      for c in range(len(confs))]

    if len(confs)>1:
      # Probabilty of jumping to a conformation k
      # is proportional to exp(-E/(R*600.)).
      logweight = np.array(energies)/(R*600.)
      weights = np.exp(-logweight+min(logweight))
      self.weights = weights/sum(weights)

      # self.darts[j][k] will jump from conformation j to conformation k
      self.darts = [[confs_BAT[j][self._BAT_to_perturb] - \
        confs_BAT[k][self._BAT_to_perturb] \
        for j in range(len(confs))] for k in range(len(confs))]

      # Finds the minimum distance between target conformations.
      # This is the maximum allowed distance to permit a dart.
      if self.extended:
        ssd = np.concatenate([[np.sum(np.square(confs_ha[j] - confs_ha[k]))
                for j in range(k+1,len(confs))] \
                  for k in range(len(confs))])
        # Uses the minimum distance or rmsd of 0.25 A
        self.epsilon = min(np.min(ssd)*3/4., confs_ha[0].shape[0]*0.025*0.025)
      else:
        period_fracs = [[np.abs(self.darts[j][k])/twoPi \
          for j in range(len(confs))] for k in range(len(confs))]
        period_fracs = [[np.min([period_fracs[j][k], 1-period_fracs[j][k]],0) \
          for j in range(len(confs))] for k in range(len(confs))]
        spf = np.concatenate([[\
          np.sum(period_fracs[j][k]) \
          for j in range(k+1,len(confs))] for k in range(len(confs))])
        self.epsilon = np.min(spf)*3/4.
    else:
      self.epsilon = 0.

    # Set the universe to the lowest-energy configuration
    if MMTK:
      self.universe.setConfiguration(Configuration(self.universe,np.copy(confs[0])))
    else:
      self.universe.context.setPositions(np.copy(confs[0]))
    self.confs = confs
    self.confs_ha = confs_ha
    self.confs_BAT = confs_BAT
    self.confs_BAT_tp = confs_BAT_tp
    self.period_frac_threshold = period_frac_threshold

    from AlGDock.BindingPMF import HMStime
    report = "  attempted %d and finished with" + \
      " %d smart darting targets (eps=%.4f) in " + \
      HMStime(time.time()-start_time)
    report = report%(nconfs_attempted, len(self.confs), self.epsilon)
    if len(self.confs)>1:
      report += "\n  the lowest %d energies are: "%(min(len(confs),10)) + \
        ', '.join(['%.2f'%e for e in energies[:10]])
    return report

  def show_confs(self, confs=None):
    if confs==None:
      if self.extended:
        confs = self.confs
      else:
        confs = [self._BAT_util.Cartesian(bat) for bat in self.confs_BAT]
    import AlGDock.IO
    IO_dcd = AlGDock.IO.dcd(self.molecule)
    IO_dcd.write('confs.dcd', confs)
    self._BAT_util.showMolecule(dcdFN='confs.dcd')
    import os
    os.remove('confs.dcd')

  def _closest_pose_Cartesian(self, conf_ha):
    # Closest pose has smallest sum of square distances between heavy atom coordinates
    ssd = np.array([(np.square(self.confs_ha[c] - conf_ha)).sum() \
      for c in range(len(self.confs_ha))])
    closest_pose_index = np.argmin(ssd)
    return (closest_pose_index, ssd[closest_pose_index])

  def _closest_pose_BAT(self, conf_BAT_tp):
    # Closest pose has smallest sum of period fractions between torsion angles
    # For only torsion angles, differences in units of periods (between 0 and 1)
    period_fracs = (np.abs(np.array([self.confs_BAT_tp[c] - conf_BAT_tp \
      for c in range(len(self.confs_BAT_tp))]))%twoPi)/twoPi
    # Wraps around the period
    period_fracs = np.min([period_fracs,1-period_fracs],0)
    spf = np.sum(period_fracs,1)
    closest_pose_index = np.argmin(spf)
    return (closest_pose_index, spf[closest_pose_index])

  def _p_attempt(self, dart_from, dart_to):
    return self.weights[dart_to]/(1.-self.weights[dart_from])

  def __call__(self, **options):
    if (self.confs is None) or len(self.confs)<2:
      return ([self.universe.configuration().array], [self.universe.energy()], \
        0, 0, 0.0)
    
    # Process the keyword arguments
    self.setCallOptions(options)
    # Check if the universe has features not supported by the integrator
    Features.checkFeatures(self, self.universe)
  
    RT = R*self.getOption('T')
    ntrials = min(self.getOption('ntrials'), len(self.confs))

    # Seed the random number generator
    if 'random_seed' in self.call_options.keys():
      np.random.seed(self.getOption('random_seed'))

    acc = 0.
    energies = []
    closest_poses = []
    if MMTK:
      xo_Cartesian = np.copy(self.universe.configuration().array)
    else:
      state = self.universe.context.getState(getPositions=True, getEnergy=True)
      xo_Cartesian = state.getPositions(asNumpy=True)  # Returns positions as a NumPy array
      energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

    if self.extended:
      (closest_pose_o, distance_o) = self._closest_pose_Cartesian(\
        xo_Cartesian[self.molecule.heavy_atoms,:])
    else:
      xo_BAT = self._BAT_util.BAT(xo_Cartesian, extended=self.extended)
      (closest_pose_o, distance_o) = self._closest_pose_BAT(xo_BAT[self._BAT_to_perturb])

    # Only attempt smart darting within
    # a sphere of radius epsilon of the minimum
    if distance_o > self.epsilon:
      if MMTK:
        return ([self.universe.configuration().array], [self.universe.energy()], \
        0, 0, 0.0)
      else:
        return ([xo_Cartesian], [energy], \
        0, 0, 0.0)

    if self.extended:
      xo_BAT = self._BAT_util.BAT(xo_Cartesian, extended=self.extended)
    if MMTK:
      eo = self.universe.energy()
    else:
      eo = energy

#    report = ''
    for t in range(ntrials):
      # Choose a pose to dart towards
      dart_towards = closest_pose_o
      while dart_towards==closest_pose_o:
        dart_towards = np.random.choice(len(self.weights), p=self.weights)
      # Generate a trial move
      xn_BAT = np.copy(xo_BAT)
      xn_BAT[self._BAT_to_perturb] = xo_BAT[self._BAT_to_perturb] + self.darts[closest_pose_o][dart_towards]

      # Check that the trial move is closest to dart_towards
      if self.extended:
        xn_Cartesian = self._BAT_util.Cartesian(xn_BAT)
        (closest_pose_n, distance_n) = self._closest_pose_Cartesian(\
          xn_Cartesian[self.molecule.heavy_atoms,:])
        if (closest_pose_n!=dart_towards):
          print '    attempted dart from pose %d (%f) to pose %d'%(\
            closest_pose_o,distance_o,dart_towards) + \
            ' landed near pose %d (%f)!'%(closest_pose_n, distance_n)
          continue
      else:
        (closest_pose_n, distance_n) = self._closest_pose_BAT(\
          xn_BAT[self._BAT_to_perturb])
        if (closest_pose_n!=dart_towards):
          print '    attempted dart from pose %d (%f) to pose %d'%(\
            closest_pose_o,distance_o,dart_towards) + \
            ' landed near pose %d (%f)!'%(closest_pose_n, distance_n)
          continue
        xn_Cartesian = self._BAT_util.Cartesian(xn_BAT)

      # Determine energy of new state
      if MMTK:
        self.universe.setConfiguration(Configuration(self.universe, xn_Cartesian))
        en = self.universe.energy()
      else:
        self.universe.context.setPositions(xn_Cartesian)
        state = self.universe.context.getState(getPositions=True, getEnergy=True)
        en = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
      #      report += 'Attempting move from near pose %d with energy %f to pose %d with energy %f. '%(closest_pose_o,eo,closest_pose_n,en)

      # Accept or reject the trial move
      if (abs(en-eo)<1000) and \
          ((en<eo) or (np.random.random()<np.exp(-(en-eo)/RT)*\
            self._p_attempt(closest_pose_n,closest_pose_o)/\
            self._p_attempt(closest_pose_o,closest_pose_n))):
        xo_Cartesian = xn_Cartesian
        xo_BAT = xn_BAT
        eo = 1.*en
        closest_pose_o = closest_pose_n
        distance_o = distance_n
        acc += 1
#        report += 'Accepted.\n'
#      else:
#        report += 'Rejected.\n'
#
#    print report
    if MMTK:
      self.universe.setConfiguration(Configuration(self.universe, xo_Cartesian))
    else:
      self.universe.context.setPositions(xo_Cartesian)
    return ([np.copy(xo_Cartesian)], [eo], acc, ntrials, 0.0)