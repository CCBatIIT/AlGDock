# This module implements a Smart Darting "integrator"

from MMTK import Configuration, Dynamics, Environment, Features, Trajectory, Units
import MMTK_dynamics
import numpy as np

R = 8.3144621*Units.J/Units.mol/Units.K

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
    import AlGDock.RigidBodies
    id = AlGDock.RigidBodies.identifier(self.universe, self.molecule)
    from BAT import converter
    self._BAT_util = converter(self.universe, self.molecule, \
      initial_atom = id.initial_atom)
    # BAT coordinates to perturb: external coordinates and primary torsions
    dof = self.universe.configuration().array.shape[0]*3 if extended \
      else self.universe.configuration().array.shape[0]*3 - 6
    self._BAT_to_perturb = range(6) if extended else []
    self._BAT_to_perturb += self._BAT_util.getFirstTorsionInds(extended)
#    self._BAT_to_perturb += list(np.array(\
#      self._BAT_util.getFirstTorsionInds(extended))-1)
#    self._BAT_to_perturb += list(np.array(\
#      self._BAT_util.getFirstTorsionInds(extended))-2)
#    self._BAT_to_perturb.sort()

    if confs is None:
      self.confs = None
    else:
      self.set_confs(confs)

  def set_confs(self, confs, rmsd_threshold=0.05, period_threshold=0.3, \
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
    from MMTK.Minimization import SteepestDescentMinimizer # @UnresolvedImport
    minimizer = SteepestDescentMinimizer(self.universe)

    minimized_confs = []
    minimized_energies = []
    for conf in confs:
      self.universe.setConfiguration(Configuration(self.universe, conf))
      x_o = np.copy(self.universe.configuration().array)
      e_o = self.universe.energy()
      for rep in range(20):
        minimizer(steps = 50)
        x_n = np.copy(self.universe.configuration().array)
        e_n = self.universe.energy()
        diff = abs(e_o-e_n)
        if np.isnan(e_n) or diff<1. or diff>1000.:
          self.universe.setConfiguration(Configuration(self.universe, x_o))
          break
        else:
          x_o = x_n
          e_o = e_n
      if not np.isnan(e_o):
        minimized_confs.append(x_o)
        minimized_energies.append(e_o)
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
      for j in range(len(confs)):
        min_rmsd = np.min([np.sqrt((confs[j][self.molecule.heavy_atoms,:] - \
          confs[k][self.molecule.heavy_atoms,:])**2).sum()/self.molecule.nhatoms \
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
      # Keep only unique configurations, using period as a threshold
      inds_to_keep = [0]
      for j in range(len(confs)):
        diffs = np.array([(confs_BAT_tp[j] - confs_BAT_tp[k]) \
          for k in inds_to_keep])/(2*np.pi)
        diffs = np.min([diffs,1-diffs],0)
        if (np.max(np.abs(diffs),1)>period_threshold).all():
          inds_to_keep.append(j)
      confs = [confs[i] for i in inds_to_keep]
      confs_BAT = [confs_BAT[i] for i in inds_to_keep]
      confs_BAT_tp = [confs_BAT_tp[i] for i in inds_to_keep]
      energies = [energies[i] for i in inds_to_keep]

    self.confs = confs
    self.confs_ha = [self.confs[c][self.molecule.heavy_atoms,:] \
      for c in range(len(self.confs))]

    if len(self.confs)>1:
      # Probabilty of jumping to a conformation k
      # is proportional to exp(-E/(R*600.)).
      logweight = np.array(energies)/(R*600.)
      weights = np.exp(-logweight+min(logweight))
      self.weights = weights/sum(weights)

      self.confs_BAT = confs_BAT
      self.confs_BAT_tp = confs_BAT_tp
      # self.darts[j][k] will jump from conformation j to conformation k
      self.darts = [[self.confs_BAT[j][self._BAT_to_perturb] - \
        self.confs_BAT[k][self._BAT_to_perturb] \
        for j in range(len(confs))] for k in range(len(confs))]

    # Set the universe to the lowest-energy configuration
    self.universe.setConfiguration(Configuration(self.universe,np.copy(confs[0])))

    from AlGDock.BindingPMF import HMStime
    report = "  started with %d, attempted %d, and ended with" + \
      " %d smart darting targets in " + HMStime(time.time()-start_time)
    report = report%(nconfs_o, nconfs_attempted, len(self.confs))
    report += "\n  the lowest %d energies are: "%min(len(confs),10) + \
        ', '.join(['%.2f'%e for e in energies[:10]])
    return report

  def show_confs(self):
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
    return np.argmin(np.array([((self.confs_ha[c] - conf_ha)**2).sum() \
      for c in range(len(self.confs))]))

  def _closest_pose_BAT(self, conf_BAT_tp):
    # Closest pose has smallest sum of square distances between torsion angles
    # For only torsion angles, differences in units of periods (between 0 and 1)
    diffs = np.abs(np.array([self.confs_BAT_tp[c] - conf_BAT_tp \
      for c in range(len(self.confs_BAT_tp))]))/(2*np.pi)
    # Wraps around the period
    diffs = np.min([diffs,1-diffs],0)
    return np.argmin(np.sum(diffs**2,1))

  def _p_attempt(self, dart_from, dart_to):
    return self.weights[dart_to]/(1.-self.weights[dart_from])

  def __call__(self, **options):
    if (self.confs is None) or len(self.confs)<2:
      return ([self.universe.configuration().array], [self.universe.energy()], 0.0, 0.0)
    
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

    xo_Cartesian = np.copy(self.universe.configuration().array)
    xo_BAT = self._BAT_util.BAT(xo_Cartesian, extended=self.extended)
    eo = self.universe.energy()
    if self.extended:
      closest_pose_o = self._closest_pose_Cartesian(\
        xo_Cartesian[self.molecule.heavy_atoms,:])
    else:
      closest_pose_o = self._closest_pose_BAT(xo_BAT[self._BAT_to_perturb])
  
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
        closest_pose_n = self._closest_pose_Cartesian(\
          xn_Cartesian[self.molecule.heavy_atoms,:])
        if (closest_pose_n!=dart_towards):
          continue
      else:
        closest_pose_n = self._closest_pose_BAT(xn_BAT[self._BAT_to_perturb])
        if (closest_pose_n!=dart_towards):
          continue
        xn_Cartesian = self._BAT_util.Cartesian(xn_BAT)

      # Determine energy of new state
      self.universe.setConfiguration(Configuration(self.universe, xn_Cartesian))
      en = self.universe.energy()
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
        acc += 1
#        report += 'Accepted.\n'
#      else:
#        report += 'Rejected.\n'
#
#    print report
    self.universe.setConfiguration(Configuration(self.universe, xo_Cartesian))
    return ([np.copy(xo_Cartesian)], [eo], float(acc)/float(ntrials), 0.0)