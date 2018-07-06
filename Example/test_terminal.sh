#!/bin/bash

# time python
# kernprof -l -v
# python -m memory_profiler

# The example is 1of6

$ALGDOCK --dir_dock dock --dir_cool cool \
  --ligand_tarball prmtopcrd/ligand.tar.gz \
  --ligand_database prmtopcrd/ligand.db \
  --forcefield prmtopcrd/gaff2.dat \
  --ligand_prmtop prmtopcrd/ligand.prmtop \
  --ligand_inpcrd prmtopcrd/ligand.trans.inpcrd \
  --ligand_mol2 prmtopcrd/ligand.mol2 \
  --ligand_rb prmtopcrd/ligand.rb \
  --receptor_tarball prmtopcrd/receptor.tar.gz \
  --receptor_prmtop prmtopcrd/receptor.prmtop \
  --receptor_inpcrd prmtopcrd/receptor.trans.inpcrd \
  --receptor_fixed_atoms prmtopcrd/receptor.pdb \
  --complex_tarball prmtopcrd/complex.tar.gz \
  --complex_prmtop prmtopcrd/complex.prmtop \
  --complex_inpcrd prmtopcrd/complex.trans.inpcrd \
  --complex_fixed_atoms prmtopcrd/complex.pdb \
  --score prmtopcrd/xtal_plus_dock6_scored.mol2 \
  --pose -1 \
  --rmsd \
  --dir_grid grids \
  --protocol Adaptive --cool_therm_speed 15.0 --dock_therm_speed 1.0 \
  --sampler HMC \
  --solvation Desolvated \
  --seeds_per_state 10 --steps_per_seed 200 --darts_per_seed 0 \
  --sweeps_per_cycle 50 --snaps_per_cycle 10 --attempts_per_sweep 100 \
  --steps_per_sweep 50 --darts_per_sweep 0 \
  --cool_repX_cycles 3 --dock_repX_cycles 4 \
  --site Sphere \
  --site_center 1.7416 1.7416 1.7416 \
  --site_max_R 1.0 \
  --site_density 10. \
  --phases NAMD_Gas NAMD_OBC \
  --cores -1 \
  --random_seed -1 \
  --run_type timed \
  --max_time 240

# Good thermodynamic speeds for testing:
# --protocol Adaptive --cool_therm_speed 0.5 --dock_therm_speed 0.5

# Original spherical site
#  --site Sphere --site_center 1.74395 1.74395 1.74395 \
#  --site_max_R 0.6 \
