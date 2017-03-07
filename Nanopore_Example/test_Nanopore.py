import AlGDock.Nanopore
from AlGDock.Nanopore import *

# This will run a small molecule
#  self = AlGDock.Nanopore.Simulation(\
#    ligand_tarball='prmtopcrd/ligand.tar.gz', \
#    ligand_database='ligand.db', \
#    forcefield='prmtopcrd/gaff2.dat', \
#    frcmodList=['ligand.frcmod'], \
#    ligand_prmtop='ligand.prmtop', \
#    grid_LJr='grids/LJr.nc', \
#    starting_conf='prmtopcrd/anchor_and_grow_scored.mol2', \
#    ef=1.0E8, \
#    max_trials=10000, \
#    report_interval=100)

# This will run a peptide
self = AlGDock.Nanopore.Simulation(\
  # These will stay the same
  forcefield='prmtopcrd/parm10.dat', \
  frcmodList=['prmtopcrd/frcmod.ff14SB'], \
  # These will change from peptide to peptide
  ligand_database='prmtopcrd/YFF_14.db', \
  ligand_prmtop='prmtopcrd/YFF_14.prmtop', \
  starting_conf='prmtopcrd/YFF_14.inpcrd', \
  # The grids will stay the same for the same protein
  grid_LJr='grids/LJr.nc', \
  ef=1.0E8, \
  max_trials=10000, \
  report_interval=100)

self.prep_ligand()
self.run()
self.view_trajectory()
