# This module implements a Hamiltonian Monte Carlo "integrator"
# It is a stripped-down version of the velocity verlet integrator
# with a Metropolis acceptance criterion.
# It requires the option 'T' in addition to velocity verlet options.

from MMTK import Dynamics, Environment, Features, Trajectory, Units, Random
import bTD_dynamics
from Scientific import N

import os
import time
from time import sleep
import signal
import subprocess
import random
import numpy as np

R = 8.3144621*Units.J/Units.mol/Units.K

#
# Torsional Dynamics Hamiltonian Monte Carlo integrator
#
class TDHamiltonianMonteCarloIntegrator(Dynamics.Integrator):

    def __init__(self, universe, server="a", ligand_dir="b", dyn_type="c", seed=int(time.time() + os.getpid()), **options):
        Dynamics.Integrator.__init__(self, universe, options)
        # Supported features: none for the moment, to keep it simple
        self.features = []

    def __call__(self, **options):
        # Process the keyword arguments
        self.setCallOptions(options)

        random.seed(self.getOption('seed'))
        #print "seed_0=   ", self.getOption('seed')
        np.random.seed(self.getOption('seed'))
        memid = np.random.randint(2000,9999)
        #print "seed_np1= ", self.getOption('seed')
        #print " RAND memid= " + str(memid)
        np.random.seed(int(time.time() + os.getpid()))
        #print "seed_np2= ", int(time.time() + os.getpid())

        Random.initializeRandomNumbersFromTime()

        # Check if the universe has features not supported by the integrator
        Features.checkFeatures(self, self.universe)
      
        RT = R*self.getOption('T')
        delta_t = self.getOption('delta_t')

        server = self.getOption('server')
        ligand_dir = self.getOption('ligand_dir')
        dyn_type = self.getOption('dyn_type')
       
        server_cmd = ''
        server_cmd = server + ' -mol2 ' + ligand_dir + 'ligand.mol2' +  ' -rb ' + ligand_dir + 'ligand.rb' + ' -gaff gaff.dat -frcmod ' + ligand_dir + 'ligand.frcmod -ictd ' + dyn_type + ' -memid ' + str(memid) + ' > tdout'
        #print server_cmd
        tdserver = subprocess.Popen(server_cmd, shell=True, preexec_fn=os.setsid)
        sleep(0.1)
 
        if 'steps_per_trial' in self.call_options.keys():
          steps_per_trial = self.getOption('steps_per_trial')
          ntrials = self.getOption('steps')/steps_per_trial
        else:
          steps_per_trial = self.getOption('steps')
          ntrials = 1
  
        if 'normalize' in self.call_options.keys():
          normalize = self.getOption('normalize')
        else:
          normalize = False          
      
        # Get the universe variables needed by the integrator
        masses = self.universe.masses()
        fixed = self.universe.getAtomBooleanArray('fixed')
        nt = self.getOption('threads')
        comm = self.getOption('mpi_communicator')
        evaluator = self.universe.energyEvaluator(threads=nt,
                                                  mpi_communicator=comm)
        evaluator = evaluator.CEvaluator()

 
        # Deal with reordering
        #print "objectList = ", self.universe.objectList()
        #print "objectList prmtop_order = ", self.universe.objectList()[0].prmtop_order
        #print "atomList = ", self.universe.objectList()[0].atomList()
        names = []
        order = []
        for a in self.universe.objectList()[0].atomList():
          names.append(a.name)
        for i in self.universe.objectList()[0].prmtop_order:
          order.append(i.number)
        #print "names = ", names
        #print "order = ", order
        descorder = list('')
        for o in order:
          descorder = descorder + list(str(o))
          descorder = descorder + list('|')
        descorder = descorder + list('1')
        descorder = descorder + list('|')
        descorder = descorder + list(str(memid))
        descorder = descorder + list('|')
        #print descorder

      
        #print 'Py T ', self.getOption('T'), ' dt ', delta_t, ' trls ', ntrials, 'stps/trl ', steps_per_trial 

        xs = []
        energies = []
        k = 0
        acc = 0
        for t in range(ntrials):
          # Initialize the velocity
          self.universe.initializeVelocitiesToTemperature(self.getOption('T'))
          #print "Python vels= ",
          #for t in self.universe.velocities(): print t,
          #print

          # Store previous configuration and initial energy
          xo = self.universe.copyConfiguration()
          pe_o = self.universe.energy()
          eo = pe_o + self.universe.kineticEnergy()
          ke_o = self.universe.kineticEnergy()

          # Run the velocity verlet integrator
          late_args = (
                  masses.array, fixed.array, evaluator,
                  N.zeros((0, 2), N.Int), N.zeros((0, ), N.Float),
                  N.zeros((1,), N.Int),
                  N.zeros((0,), N.Float), N.zeros((2,), N.Float),
                  N.zeros((0,), N.Float), N.zeros((1,), N.Float),
                  delta_t, self.getOption('first_step'),
                  steps_per_trial, self.getActions(),
                  ''.join(descorder))
          #self.run(bTD_dynamics.integrateRKM_SPLIT,
          #tdserver = subprocess.Popen(server_cmd, shell=True, preexec_fn=os.setsid)
          #sleep(0.1)
          self.run(bTD_dynamics.integrateRKM_INTER3,
            (self.universe,
             self.universe.configuration().array,
             self.universe.velocities().array) + late_args)

          # Decide whether to accept the move
          pe_n = self.universe.energy()
          ke_n = self.universe.kineticEnergy()
          en = pe_n + ke_n
          random_number = np.random.random()
          expression = np.exp(-(en-eo)/RT)
          #expression = np.exp(-(pe_n-pe_o)/RT)
          #print ' stps Py PEo ', pe_o, ' KEo ', ke_o, ' Eo ', eo, ' PEn ', pe_n, ' KEn ', ke_n, ' En ', en, ' rand ', random_number, ' exp ', expression
          if (((en<eo) or (random_number<expression)) and (not np.isnan(en))):
          #if ((pe_n<pe_o) or (random_number<expression)): #or (t==0): # Always acc 1st step for now
            acc += 1
            #print 'Py MC acc'
            if normalize:
              self.universe.normalizeConfiguration()
            descorder[-7] = '1'
          else:
            #print 'Py MC nacc'
            self.universe.setConfiguration(xo)
            pe_n = pe_o
            descorder[-7] = '0'
          xs.append(self.universe.copyConfiguration().array)
          energies.append(pe_n)

        # Kill the server
        tdserver.terminate()
        os.killpg(tdserver.pid, signal.SIGTERM)  # Send the signal to all the process groups
        os.system('ipcrm -M ' + str(memid))      # Free shared memory
        
        #Return
        return (xs, energies, float(acc)/float(ntrials), delta_t)
