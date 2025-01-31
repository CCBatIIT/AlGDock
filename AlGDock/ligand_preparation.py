from AlGDock.logger import NullDevice
from AlGDock.BindingPMF import HMStime

import os
import sys
import time
import numpy as np

try:
  import MMTK
  import MMTK.Units
  from MMTK.ParticleProperties import Configuration
except ImportError:
  MMTK = None

try:
  import openmm
  import openmm.unit as unit
  from openmm.app import AmberPrmtopFile, AmberInpcrdFile, Simulation, NoCutoff
  from openmm import *
except ImportError:
  OpenMM = None

class LigandPreparation():
    """Minimizes and equilibrates the ligand

    Attributes
    ----------
    args : AlGDock.simulation_arguments.SimulationArguments
      Simulation arguments
    log : AlGDock.logger.Logger
      Simulation log
    top : AlGDock.topology.TopologyMMTK
      Topology of the ligand
    system : AlGDock.system.System
      Simulation system
    _get_confs_to_rescore : function
      Returns the configurations to rescore
    iterator : AlGDock.simulation_iterator.SimulationIterator
      Performs an iteration on one thermodynamic state
    data : AlGDock.simulation_data.SimulationData
      Stores results from the simulation
    """

    def __init__(self, args, log, top, system, _get_confs_to_rescore,
                 iterator, data):
        """Initializes the class

        Parameters
        ----------
        args : AlGDock.simulation_arguments.SimulationArguments
          Simulation arguments
        log : AlGDock.logger.Logger
          Simulation log
        top : AlGDock.topology.TopologyMMTK
          Topology of the ligand
        system : AlGDock.system.System
          Simulation system
        _get_confs_to_rescore : function
          Returns the configurations to rescore
        iterator : AlGDock.simulation_iterator.SimulationIterator
          Performs an iteration on one thermodynamic state
        data : AlGDock.simulation_data.SimulationData
          Stores results from the simulation
        """
        self.args = args
        self.log = log
        self.top = top
        self.system = system
        self._get_confs_to_rescore = _get_confs_to_rescore
        self.iterator = iterator
        self.data = data

    def run(self, process):
        """Performs the ligand preparation

        Also stores minimized ligand in milestone B and keeps track of timing

        Parameters
        ----------
        process : str
          Process, either 'BC' or 'CD'

        Returns
        -------
        seeds : list of np.array
          Minimized ligand configurations for the appropriate thermodynamic state
        """
        self.log.set_lock(process)
        self.log.recordStart(process+' ligand prep')
        self.log.tee("\n>>> Ligand preparation for %s, starting at " % process +
                     time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()))
        if process == 'BC':
            seeds = self._prepare_ligand_BC()
        elif process == 'CD':
            seeds = self._prepare_ligand_CD()
        self.log.tee("Elapsed time for ligand preparation for %s: " % process +
                     HMStime(self.log.timeSince(process+' ligand prep')))
        self.log.clear_lock(process)
        return seeds

    def _prepare_ligand_BC(self):
        """Prepares the ligand for BC simulations

        Returns
        -------
        seeds : list of np.array
          Ligand configurations minimized in milestone B
        """
        if self.data['BC'].protocol == []:

            # Set up the force field
            params_o = self.system.paramsFromAlpha(1.0, 'BC', site=False)
            self.system.setParams(params_o)

            # Get starting configurations
            basename = os.path.basename(self.args.FNs['score'])
            basename = basename[:basename.find('.')]
            dirname = os.path.dirname(self.args.FNs['score'])
            minimizedB_FN = os.path.join(dirname, basename + '_minB.nc')
            if os.path.isfile(minimizedB_FN):
                from netCDF4 import Dataset
                dock6_nc = Dataset(minimizedB_FN, 'r')
                if MMTK:
                    minimizedConfigurations = [
                    dock6_nc.variables['confs'][n][self.top.inv_prmtop_atom_order_L, :]
                    for n in range(dock6_nc.variables['confs'].shape[0])
                    ]
                else:
                    minimizedConfigurations = [
                    dock6_nc.variables['confs'][n] for n in range(dock6_nc.variables['confs'].shape[0])]

                Es = dict([(key, dock6_nc.variables[key][:])
                           for key in dock6_nc.variables.keys() if key != 'confs'])
                dock6_nc.close()
            else:
                (minimizedConfigurations, Es) = self._get_confs_to_rescore(
                    site=False, minimize=True)

                from netCDF4 import Dataset
                dock6_nc = Dataset(minimizedB_FN, 'w')
                dock6_nc.createDimension(
                    'n_confs', len(minimizedConfigurations))
                dock6_nc.createDimension(
                    'n_atoms', minimizedConfigurations[0].shape[0])
                dock6_nc.createDimension('n_cartesian', 3)
                dock6_nc.createDimension('one', 1)
                dock6_nc.createVariable(
                    'confs', 'f8', ('n_confs', 'n_atoms', 'n_cartesian'))
                for n in range(len(minimizedConfigurations)):
                    dock6_nc.variables['confs'][n] = minimizedConfigurations[n][self.top.prmtop_atom_order_L, :]
                for key in Es.keys():
                    dock6_nc.createVariable(key, 'f8', ('one', 'n_confs'))
                    dock6_nc.variables[key][:] = Es[key]
                dock6_nc.close()

            # initializes smart darting for BC
            # and sets the universe to the lowest energy configuration
            if MMTK:
                self.iterator.initializeSmartDartingConfigurations(
                minimizedConfigurations, 'BC', self.log, self.data)
                if len(minimizedConfigurations) > 0:
                    self.top.universe.setConfiguration(
                    Configuration(self.top.universe, minimizedConfigurations[-1]))
            else:
                #TODO: NEED TO MINIMIZE the data
                self.iterator.initializeSmartDartingConfigurations(
                minimizedConfigurations, 'BC', self.log, self.data) # minimize the date
                if len(minimizedConfigurations) > 0:
                    self.top.OMM_simulation.context.setPositions(minimizedConfigurations[-1])

            self.data['BC'].confs['starting_poses'] = minimizedConfigurations

            # Ramp the temperature from 0 to the desired starting temperature using HMC
            self._ramp_T(params_o['T'], normalize=True)

            # Run at starting temperature
            seeds = [np.copy(self.top.universe.configuration().array)
                     for n in range(self.args.params['BC']['seeds_per_state'])]
        else:
            seeds = None
        return seeds

    def _prepare_ligand_CD(self):
        """Prepares the ligand for CD simulations

        Returns
        -------
        seeds : list of np.array
          Ligand configurations minimized in milestone D
        """
        if self.data['CD'].protocol == []:
            params_o = self.system.paramsFromAlpha(1.0, 'CD')
            self.system.setParams(params_o)

            if (self.args.params['CD']['pose'] == -1):
                seeds = self._get_confs_to_rescore(site=True, minimize=True)[0] 
                self.data['CD'].confs['starting_poses'] = seeds
                
            else:
                # For pose BPMF, starting_poses is defined in _set_universe_evaluator
                seeds = self.data['CD'].confs['starting_poses']

            if seeds == []:
                seeds = None
            else:
                # initializes smart darting for CD and sets the universe
                # to the lowest energy configuration
                if MMTK:
                    self.iterator.initializeSmartDartingConfigurations(
                    seeds, 'CD', self.log, self.data)
                    if len(seeds) > 0:
                        self.top.universe.setConfiguration(
                        Configuration(self.top.universe, np.copy(seeds[-1])))
                else:
                    # TODO: NEED TO MINIMIZE the data
                    self.iterator.initializeSmartDartingConfigurations(
                    seeds, 'CD', self.log, self.data)
                    if len(seeds) > 0:
                        self.top.OMM_simulation.context.setPositions(
                        np.copy(seeds[-1]))

                print ('RAMPING UP')
                # Ramp up the temperature using HMC
                self._ramp_T(self.args.params['BC']
                             ['T_TARGET'], normalize=False)

                seeds = [np.copy(self.top.universe.configuration().array)
                         for n in range(self.args.params['CD']['seeds_per_state'])]
        return seeds

    def _ramp_T(self, T_START, T_LOW=20., normalize=False):
        """Ramp the temperature from T_LOW to T_START

        Parameters
        ----------
        T_START : float
          The final temperature
        T_LOW : float
          The lowest temperature in the ramp
        normalize : bool
          If True, the ligand center of mass will be normalized
        """
        self.log.recordStart('T_ramp')

        # First minimize the energy
        from MMTK.Minimization import SteepestDescentMinimizer  # @UnresolvedImport
        minimizer = SteepestDescentMinimizer(self.top.universe)

        original_stderr = sys.stderr
        sys.stderr = NullDevice()  # Suppresses warnings for minimization

        x_o = np.copy(self.top.universe.configuration().array)
        e_o = self.top.universe.energy()
        for rep in range(5000):
            minimizer(steps=10)
            x_n = np.copy(self.top.universe.configuration().array)
            e_n = self.top.universe.energy()
            diff = abs(e_o - e_n)
            if np.isnan(e_n) or diff < 0.05 or diff > 1000.:
                self.top.universe.setConfiguration(
                    Configuration(self.top.universe, x_o))
                break
            else:
                x_o = x_n
                e_o = e_n

        sys.stderr = original_stderr
        self.log.tee("  minimized to %.3g kcal/mol over %d steps" % (e_o, 10 *
                                                                     (rep + 1)))

        # Then ramp the energy to the starting temperature
        from AlGDock.Integrators.HamiltonianMonteCarlo.HamiltonianMonteCarlo \
            import HamiltonianMonteCarloIntegrator
        sampler = HamiltonianMonteCarloIntegrator(self.top.universe)

        e_o = self.top.universe.energy()
        T_LOW = 20.
        T_SERIES = T_LOW * (T_START / T_LOW)**(np.arange(30) / 29.)
        for T in T_SERIES:
            delta_t = 2.0 * MMTK.Units.fs
            steps_per_trial = 10
            attempts_left = 10
            while attempts_left > 0:
                random_seed = int(T*10000) + attempts_left + \
                    int(self.top.universe.configuration().array[0][0]*10000)
                if self.args.random_seed == 0:
                    random_seed += int(time.time() * 1000)
                random_seed = random_seed % 32767
                (xs, energies, acc, ntrials, delta_t) = \
                    sampler(steps=2500, steps_per_trial=10, T=T,
                            delta_t=delta_t, random_seed=random_seed)
                attempts_left -= 1
                acc_rate = float(acc) / ntrials
                if acc_rate < 0.4:
                    delta_t -= 0.25 * MMTK.Units.fs
                else:
                    attempts_left = 0
                if delta_t < 0.1 * MMTK.Units.fs:
                    delta_t = 0.1 * MMTK.Units.fs
                    steps_per_trial = max(int(steps_per_trial / 2), 1)
            fmt = "  T = %d, delta_t = %.3f fs, steps_per_trial = %d, acc_rate = %.3f"
            if acc_rate < 0.01:
                print self.top.universe.energyTerms()
            self.log.tee(fmt % (T, delta_t * 1000, steps_per_trial, acc_rate))
        if normalize:
            self.top.universe.normalizePosition()
        e_f = self.top.universe.energy()

        self.log.tee("  ramped temperature from %d to %d K in %s, " % (
            T_LOW, T_START, HMStime(self.log.timeSince('T_ramp'))) +
            "changing energy to %.3g kcal/mol\n" % (e_f))
