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
  --protocol Adaptive --cool_therm_speed 1.25 --dock_therm_speed 1.25 \
  --no_protocol_refinement \
  --sampler NUTS \
  --MCMC_moves 1 \
  --seeds_per_state 10 \
  --steps_per_seed 250 \
  --sweeps_per_cycle 5 --attempts_per_sweep 100 --steps_per_sweep 1000 \
  --cool_repX_cycles 2 --dock_repX_cycles 2 \
  --site Sphere --site_center 1.91650 1.91650 1.91650 \
  --site_max_R 0.01 --site_density 10. \
  --phases NAMD_Gas NAMD_GBSA \
  --cores -1 \
  --score prmtopcrd/anchor_and_grow_scored.mol2 \
  --rmsd \
  --run_type timed \
  --max_time 10

#  --no_stored_evaluators \
#  --phases NAMD_Gas NAMD_GBSA Gas GBSA APBS PBSA \
