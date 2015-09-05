import AlGDock.BindingPMF_plots
import os, shutil, glob

self = AlGDock.BindingPMF_plots.BPMF_plots(\
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
  seeds_per_state=10, steps_per_seed=200,
  sweeps_per_cycle=10, attempts_per_sweep=100, steps_per_sweep=50,
  cool_repX_cycles=2, dock_repX_cycles=3, \
  site='Sphere', site_center=[1.91650, 1.91650, 1.91650], site_max_R=0.01, \
  site_density=10., \
  phases=['NAMD_Gas','NAMD_GBSA'], \
  cores=-1, \
  rmsd=True)

(confs,Es) = self.pose_energies()

#  from AlGDock.ForceFields.Grid.TrilinearTransformGrid \
#    import TrilinearTransformGridForceField
#  for type in ['LJa','LJr']:
#    Es[type+'_TrilinearTransform'] = np.zeros((12,len(confs)),dtype=np.float)
#    for p in range(12):
#      FF = TrilinearTransformGridForceField(self._FNs['grids'][type], 1.0, \
#        'scaling_factor_'+type, grid_name='%f'%(p+1), \
#        inv_power=-float(p+1), max_val=-1)
#      self.universe.setForceField(FF)
#      for c in range(len(confs)):
#        self.universe.setConfiguration(Configuration(self.universe,confs[c]))
#        Es[type+'_transformed'][p,c] = self.universe.energy()