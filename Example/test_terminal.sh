#!/bin/bash

# time python
# kernprof -l -v
# python -m memory_profiler

$ALGDOCK --dir_dock dock --dir_cool cool \
  --ligand_database prmtopcrd/ligand.db \
  --forcefield prmtopcrd/gaff.dat \
  --ligand_prmtop prmtopcrd/ligand.prmtop \
  --ligand_inpcrd prmtopcrd/ligand.trans.inpcrd \
  --receptor_prmtop prmtopcrd/receptor.prmtop \
  --receptor_inpcrd prmtopcrd/receptor.trans.inpcrd \
  --receptor_fixed_atoms prmtopcrd/receptor.pdb \
  --complex_prmtop prmtopcrd/complex.prmtop \
  --complex_inpcrd prmtopcrd/complex.trans.inpcrd \
  --complex_fixed_atoms prmtopcrd/complex.pdb \
  --rmsd \
  --dir_grid grids \
  --protocol Adaptive --cool_therm_speed 0.5 --dock_therm_speed 0.5 \
  --sampler NUTS \
  --MCMC_moves 1 \
  --seeds_per_state 10 \
  --steps_per_seed 200 \
  --sweeps_per_cycle 10 --attempts_per_sweep 100 --steps_per_sweep 50 \
  --cool_repX_cycles 2 --dock_repX_cycles 3 \
  --site Sphere --site_center 1.91650 1.91650 1.91650 \
  --site_max_R 0.01 \
  --site_density 10. \
  --phases NAMD_Gas NAMD_GBSA \
  --cores -1 \
  --score prmtopcrd/anchor_and_grow_scored.mol2 \
  --rmsd \
  --run_type timed \
  --max_time 10

#  --score_multiple \
#  --phases NAMD_Gas NAMD_GBSA Gas GBSA APBS PBSA \
#  --site Sphere --site_center 1.91650 1.91650 1.91650 \
#  --site_max_R 0.01 \

# For cooling, there are ~9 states with a thermodynamic speed of 0.7
# For undocking, there are ~27 states with a thermodynamic speed of 0.7

# Here are example cooling free energies:
#  calculated NAMD_Gas solvation free energy of 0.669438 RT using cycles 0 to 0
#  calculated NAMD_Gas solvation free energy of 0.694241 RT using cycles 1 to 1
#  calculated NAMD_GBSA solvation free energy of -13.435825 RT using cycles 0 to 0
#  calculated NAMD_GBSA solvation free energy of -13.662444 RT using cycles 1 to 1
#  calculated cooling free energy of 60.797888 RT using MBAR for cycles 0 to 0
#  calculated cooling free energy of 60.563113 RT using MBAR for cycles 1 to 1

# For warming, there are ~10 states with a thermodynamic speed of 0.7
# Free energies are about the same.

# For undocking with a temperature protocol that follows the normal grid strength,
# there are ~22 states with a thermodynamic speed of 0.7
