# The example is 1of6

import AlGDock.BindingPMF_plots
import os, shutil, glob

os.system('rm -rf cool dock')

import cProfile
import re
cProfile.run("self = AlGDock.BindingPMF_plots.BPMF_plots(\
  dir_dock='dock', dir_cool='cool',\
  ligand_tarball='prmtopcrd/ligand.tar.gz', \
  ligand_database='ligand.db', \
  forcefield='prmtopcrd/gaff.dat', \
  ligand_prmtop='ligand.prmtop', \
  ligand_inpcrd='ligand.trans.inpcrd', \
  receptor_tarball='prmtopcrd/receptor.tar.gz', \
  receptor_prmtop='receptor.prmtop', \
  receptor_inpcrd='receptor.trans.inpcrd', \
  receptor_fixed_atoms='receptor.pdb', \
  complex_tarball='prmtopcrd/complex.tar.gz', \
  complex_prmtop='complex.prmtop', \
  complex_inpcrd='complex.trans.inpcrd', \
  complex_fixed_atoms='complex.pdb', \
  score = 'prmtopcrd/anchor_and_grow_scored.mol2', \
  dir_grid='grids', \
  protocol='Adaptive', cool_therm_speed=1.5, dock_therm_speed=1.5,\
  sampler='NUTS', \
  MCMC_moves=1, \
  seeds_per_state=10, steps_per_seed=200, darts_per_seed=10, \
  sweeps_per_cycle=200, attempts_per_sweep=100, \
  steps_per_sweep=100, darts_per_sweep=10, \
  cool_repX_cycles=3, dock_repX_cycles=4, \
  site='Sphere', site_center=[1.74395, 1.74395, 1.74395], site_max_R=0.6, \
  site_density=10., \
  phases=['NAMD_Gas','NAMD_GBSA'], \
  cores=1, \
  rmsd=True, \
  run_type='cool', random_seed=50)","cool_stats")

import pstats
p = pstats.Stats("cool_stats")
p.strip_dirs().sort_stats('tottime').print_stats(25)

os.system('rm cool_stats')

# >>> Without Smart Darting
#          7599930 function calls (7477483 primitive calls) in 30.161 seconds

# >>> With 10 darts
#    calculated NAMD_GBSA solvation free energy of -80.063250 RT using cycles 0 to 0
#    calculated NAMD_GBSA solvation free energy of -79.206205 RT using cycles 1 to 1
#    calculated NAMD_GBSA solvation free energy of -78.535756 RT using cycles 1 to 2
#    calculated cooling free energy of 34.821994 RT using MBAR for cycles 0 to 0
#
#           9815897 function calls (9642206 primitive calls) in 53.601 seconds
#
#     Ordered by: internal time
#     List reduced from 1209 to 25 due to restriction <25>
#
#     ncalls  tottime  percall  cumtime  percall filename:lineno(function)
#       2090   17.119    0.008   47.327    0.023 BindingPMF.py:2306(_sim_one_state)
#      13631   13.497    0.001   13.497    0.001 {method 'Cartesian' of 'BAT.converter' objects}
#     413920    4.013    0.000    4.013    0.000 {method 'reduce' of 'numpy.ufunc' objects}
#    2097238    2.668    0.000    2.668    0.000 {numpy.core.multiarray.array}
#       2705    2.565    0.001    2.565    0.001 {method 'BAT' of 'BAT.converter' objects}
#      22990    1.546    0.000    4.970    0.000 SmartDarting.py:126(_closest_pose_BAT)
#        400    1.194    0.003    1.194    0.003 {repX.attempt_swaps}
#          1    1.121    1.121    5.606    5.606 BindingPMF.py:688(initial_cool)
#          5    0.764    0.153    0.764    0.153 {posix.waitpid}
#       2090    0.699    0.000   24.736    0.012 SmartDarting.py:135(__call__)
#       1822    0.661    0.000    0.661    0.000 {built-in method compress}
#  17081/17080    0.618    0.000    0.799    0.000 {apply}
#      22239    0.610    0.000    2.270    0.000 {method 'choice' of 'mtrand.RandomState' objects}
#    1245567    0.494    0.000    1.918    0.000 function_base.py:786(copy)
#      22224    0.493    0.000    1.338    0.000 numeric.py:2056(allclose)
#     689277    0.414    0.000    1.325    0.000 fromnumeric.py:1281(ravel)
#     265961    0.369    0.000    2.039    0.000 fromnumeric.py:1621(sum)
#          2    0.323    0.161    0.328    0.164 IO.py:204(read)
#     689319    0.291    0.000    0.291    0.000 {method 'ravel' of 'numpy.ndarray' objects}
#     692929    0.289    0.000    0.627    0.000 numeric.py:392(asarray)
#     288595    0.250    0.000    1.613    0.000 _methods.py:23(_sum)
#  335883/172611    0.224    0.000    0.313    0.000 {getattr}
#      16856    0.203    0.000    0.654    0.000 ForceField.py:294(__call__)
#     402763    0.193    0.000    0.193    0.000 {isinstance}

# 62.326 s with python SD
# 56.905 s with cdef/57.044 seconds with def/57.183 s with cpdef. They are about the same.

# On CCB, 55.572 s without/67.014 s 2 SD/78.530 s 5 SD/79.472 s 10 SD for cool_stats
#         33.748 s without/42.017 s 2 SD/50.633 s 5 SD/52.119 s 10 SD for sim_one_state

