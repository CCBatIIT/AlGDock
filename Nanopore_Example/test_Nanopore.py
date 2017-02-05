import AlGDock.Nanopore
from AlGDock.Nanopore import *

self = AlGDock.Nanopore.Simulation(\
  ligand_tarball='prmtopcrd/ligand.tar.gz', \
  ligand_database='ligand.db', \
  forcefield='prmtopcrd/gaff.dat', \
  frcmodList=['ligand.frcmod'], \
  ligand_prmtop='ligand.prmtop', \
  grid_LJr='grids/LJr.nc', \
  starting_conf='prmtopcrd/anchor_and_grow_scored.mol2', \
  ef=1.0E8, \
  max_trials=10000, \
  report_interval=100)

self.run()

self.view_trajectory()

