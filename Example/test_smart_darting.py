import AlGDock.BindingPMF_plots
import os, shutil, glob
import MMTK
import numpy as np

self = AlGDock.BindingPMF_plots.BPMF_plots(\
    dir_dock='dock', dir_cool='cool',\
    ligand_database='prmtopcrd/ligand.db', \
    forcefield='prmtopcrd/gaff.dat', \
    ligand_prmtop='prmtopcrd/ligand.prmtop', \
    ligand_inpcrd='prmtopcrd/ligand.trans.inpcrd', \
    receptor_prmtop='prmtopcrd/receptor.prmtop', \
    receptor_inpcrd='prmtopcrd/receptor.trans.inpcrd', \
    receptor_fixed_atoms='prmtopcrd/receptor.pdb', \
    complex_prmtop='prmtopcrd/complex.prmtop', \
    complex_inpcrd='prmtopcrd/complex.trans.inpcrd', \
    complex_fixed_atoms='prmtopcrd/complex.pdb', \
    score = 'prmtopcrd/anchor_and_grow_scored.mol2', \
    dir_grid='grids', \
    protocol='Adaptive', cool_therm_speed=1.5, dock_therm_speed=1.5,\
    sampler='NUTS', \
    MCMC_moves=1, \
    seeds_per_state=10, steps_per_seed=200,
    sweeps_per_cycle=10, attempts_per_sweep=100, steps_per_sweep=50,
    cool_repX_cycles=2, dock_repX_cycles=3, \
    site='Sphere', site_center=[1.91650, 1.91650, 1.91650], site_max_R=0.01, \
    site_density=10., \
    phases=['NAMD_Gas','NAMD_GBSA'], \
    cores=-1, \
    rmsd=True)

# These are cartesian coordinates of docked poses
(confs, Es) = self._get_confs_to_rescore(minimize=False)

#  len(confs)
#  for key in Es: #Es has a lot of keys, including 'total'. The 'total'key shows the energies associated with each docked pose
#      print key
#  Es['total'] # Use this to obtain the individual values of energy for each docked pose
#  len(Es['total']) # 180 values of energy

# Convert to bond-angle-torsion (BAT)
# Figure out how to move from one docked pose to another. This is your MCMC move set.

import AlGDock.BAT
c = AlGDock.BAT.converter(self.universe, self.molecule)
#  x = c.BAT() # x stores the universe BAT coordinates
#  x = c.BAT(Cartesian=confs[0]) # x stores the universe BAT coordinates

converted_confs_BAT = []
for i in range(len(confs)):
  converted_confs_BAT.append(c.BAT(Cartesian=confs[i]))

# Basically need to loope over docked poses to get "darts."
# _get_confs_to_rescore will get you docked poses in AlGDock

# What pose are you in?
# Your current configuration is: self.universe.configuration().array
# Get root mean square deviation (rmsd) values between your current pose and pose of interest. Here a quick one-liner that will do it:
# np.array([np.sqrt(((conf_c - confs[c][self.molecule.heavy_atoms,:])**2).sum()/self.molecule.nhatoms) for c in range(len(confs))])
""" The line above shows an example of "comprehensive for" structure"""

def calculate_RMSD_confs(confs):
  current_conf = self.universe.configuration().array
  RMSD_of_confs = []
  for docked_pose in range(0,len(confs)):
    RMSD_value = 0
    for i in range(0,len(confs[docked_pose])):
      RMSD_value = RMSD_value + (current_conf[i][0]-confs[docked_pose][i][0])**2 + (current_conf[i][1]-confs[docked_pose][i][1])**2 + (current_conf[i][2]-confs[docked_pose][i][2])**2
    RMSD_value = np.sqrt(RMSD_value/len(confs[docked_pose]))
    RMSD_of_confs.append(RMSD_value)
  return RMSD_of_confs

def calculate_closest_pose(confs):
    return np.argmin(calculate_RMSD_confs(confs))

# All of the following are part of the MCMC move:

dart_from = calculate_closest_pose(confs)

# Chooses one of the configurations to jump to with probability proportional
# to exp(-E/RT).
weights = np.exp(-np.array(Es['total'])/(8.314*300))
weights[dart_from] = 0
weights = weights/sum(weights)
dart_to = np.random.choice(len(weights), p=weights)

# Here is what you will need to add to the current BAT to do the move:
differece_array = np.array(converted_confs_BAT[dart_to])-np.array(converted_confs_BAT[dart_from])

# 1. Select elements of the difference array that correspond to external
# coordinates and torsions. This makes your "dart"
# 2. Add these darts to the current BAT
# 3. Convert the BAT of the new configuration to Cartesian.
# 4. Evaluate the energy of the new configuration.
# 5. Evaluate the probability of jumping from the new configuration
#    to the original configuration.
#    (a) What is the closest pose to the new configuration. dart_from_n
#    (b) What is the probability of going from dart_from_n to dart_from



# To satisfy the Metropolis criterion, you will need: weights[dart_to]
# Here is an example of an MCMC move
#    acc = 0
#    xo = self.universe.copyConfiguration()
#    eo = self.universe.energy()
#    com = self.universe.centerOfMass().array
#    
#    for c in range(trials):
#      step = np.random.randn(3)*step_size
#      xn = np.dot((xo.array - com), random_rotate()) + com + step
#      self.universe.setConfiguration(Configuration(self.universe,xn))
#      en = self.universe.energy()
#      # The following line is the Metropolis acceptance criterion:
#      if ((en<eo) or (np.random.random()<np.exp(-(en-eo)/self.RT))):
#        acc += 1
#        xo = self.universe.copyConfiguration()
#        eo = en
#        com += step
#      else:
#        self.universe.setConfiguration(xo)


# Select a new pose to dart to with a probability that satisfies detailed balance. Look up Metropolis acceptance criterion. Make sure it includes the proposal probability.
# See _MC_translate_rotate in AlGDock.BindingPMF.BPMF for an example for how to do an MCMC move.

""" THIS IS HOW YOU CREATE A PBD FILE TO LATER RUN ON VMD:
    from MMTK.PDB import PDBOutputFile
    pdb_file = PDBOutputFile('test.pdb')
    pdb_file.write(self.universe, self.universe.configuration() )
    pdb_file.close() """

""" This is how to run a short MD simulation
from NUTS import NUTSIntegrator  # @UnresolvedImport
NUTS_sampler = NUTSIntegrator(self.universe)
NUTS_sampler(steps=500,T=300) """

