import os, time
import numpy as np

from desolvationGrid_util import *
# Most of the heavy lifting is done by these Cython routines
# from desolvationGrid_util:
# enumerate_SAS_points
# set_inside_sphere_to
# increment_inside_sphere
# decrement_inside_sphere
# fraction_r4inv_low_dielectric
# calc_desolvationGrid

def goldenSectionSpiral(n):
  golden_angle = np.pi * (3 - np.sqrt(5))
  theta = golden_angle * np.arange(n)
  z = np.linspace(1 - 1.0 / n, 1.0 / n - 1, n)
  radius = np.sqrt(1 - z * z)
  unit_sphere_pts = np.transpose(np.vstack(\
    [radius * np.cos(theta), radius * np.sin(theta), z]))
  return unit_sphere_pts

class desolvationGridCalculation:
  def __init__(self, **kwargs):
    ### Parse parameters
    self.FNs = {'prmtop':kwargs['prmtop_FN'], \
      'inpcrd':kwargs['inpcrd_FN'], \
      'header':kwargs['header_FN'], \
      'grid':'desolv.dx' if kwargs['grid_FN'] is None else kwargs['grid_FN']}
    if 'counts' in kwargs.keys():
      counts = kwargs['counts']
    if 'spacing' in kwargs.keys():
      spacing = kwargs['spacing']
    
    # Check that input files are available
    for FN in [self.FNs['prmtop'],self.FNs['inpcrd']]:
      if not os.path.exists(FN):
        raise Exception(FN+' missing!')
  
    # Check that output directories are available
    for FN in [self.FNs['grid']]:
      dirN = os.path.dirname(FN)
      if dirN!='' and (not os.path.isdir(dirN)):
        os.system('mkdir -p '+os.path.dirname(FN))

    ### Read header from dx file
    header = {}
    if (self.FNs['header'] is not None) and os.path.isfile(self.FNs['header']):
      print 'Reading header from '+self.FNs['header']
      headerF = open(self.FNs['header'],'r')
      headerData = headerF.read()
      headerF.close()
      
      headerLines = headerData.split('\n')
      if counts is None:
        counts = np.array([int(x) for x in headerLines.pop(0).split(' ')[-3:]])
      for name in ['origin','d0','d1','d2']:
        header[name] = [float(x) for x in headerLines.pop(0).split(' ')[-3:]]
      if spacing is None:
        spacing = np.array([header['d0'][0], header['d1'][1], header['d2'][2]])
      del headerF, headerLines

    ### Loads coordinates
    import AlGDock.IO
    IO_crd = AlGDock.IO.crd()
    self.crd = IO_crd.read(self.FNs['inpcrd'])

    ### Load Lennard-Jones radii
    IO_prmtop = AlGDock.IO.prmtop()
    prmtop = IO_prmtop.read(self.FNs['prmtop'], \
      varnames=['POINTERS','CHARGE','ATOM_TYPE_INDEX','NONBONDED_PARM_INDEX',\
                'LENNARD_JONES_ACOEF','LENNARD_JONES_BCOEF'])
    
    prmtop['CHARGE'] = prmtop['CHARGE']/18.2223 # Convert to units of electric charge

    NATOM = prmtop['POINTERS'][0]
    NTYPES = prmtop['POINTERS'][1]

    # Extract Lennard-Jones well depth and radii for each atom type
    LJ_radius = np.zeros(shape=(NTYPES), dtype=float)
    LJ_depth = np.zeros(shape=(NTYPES), dtype=float)
    for i in range(NTYPES):
      LJ_index = prmtop['NONBONDED_PARM_INDEX'][NTYPES*i+i]-1
      if prmtop['LENNARD_JONES_ACOEF'][LJ_index]<1.0e-6:
        LJ_radius[i] = 0
        LJ_depth[i] = 0
      else:
        factor = 2 * prmtop['LENNARD_JONES_ACOEF'][LJ_index] / prmtop['LENNARD_JONES_BCOEF'][LJ_index]
        LJ_radius[i] = pow(factor, 1.0/6.0) * 0.5
        LJ_depth[i] = prmtop['LENNARD_JONES_BCOEF'][LJ_index] / 2 / factor

    # Lennard Jones and SAS radii per atom
    LJ_r = np.array([LJ_radius[prmtop['ATOM_TYPE_INDEX'][atom_index]-1] \
      for atom_index in range(NATOM)])
    self.LJ_r2 = LJ_r*LJ_r
    self.SAS_r = LJ_r + kwargs['probe_radius']
    
    self.unit_sphere_pts = goldenSectionSpiral(kwargs['SAS_points'])
    
    # Outputs files and parameters
    print '*** Files and parameters ***'
    print 'Input AMBER prmtop      :\t' + self.FNs['prmtop']
    print 'Input AMBER inpcrd      :\t' + self.FNs['inpcrd']
    if self.FNs['header'] is not None:
      print 'Input grid header file  :\t' + self.FNs['header']
    print 'Output grid             :\t' + self.FNs['grid']
    print 'Grid spacing            :\t', spacing
    print 'Grid counts             :\t', counts
    print

    kwargs['counts'] = counts
    kwargs['spacing'] = spacing
    self.kwargs = kwargs

  def calc_receptor_SAS_points(self):
    print 'Finding receptor SAS points'
    startTime = time.time()

    self.receptor_SAS_points = enumerate_SAS_points(self.crd, self.crd, \
      self.unit_sphere_pts, self.SAS_r, self.LJ_r2)

    endTime = time.time()
    print ' in %3.2f s'%(endTime-startTime)

  def calc_receptor_MS(self):
    print 'Determining the number of SAS points marking each grid point'
    startTime = time.time()
    
    self.receptor_MS_grid = np.ones(shape=tuple(self.kwargs['counts']), \
      dtype=np.int)
    # Tentatively assign the grid inside the SAS to low dielectric
    for atom_index in range(len(self.SAS_r)):
      set_inside_sphere_to(self.receptor_MS_grid, self.kwargs['spacing'], \
        self.kwargs['counts'],
        self.crd[atom_index,0], self.crd[atom_index,1], self.crd[atom_index,2], \
        self.SAS_r[atom_index], 0)
    # Determine number of SAS points marking each grid point
    for SAS_point in self.receptor_SAS_points:
      increment_inside_sphere(self.receptor_MS_grid, self.kwargs['spacing'], \
        self.kwargs['counts'], SAS_point[0], SAS_point[1], SAS_point[2], \
        self.kwargs['probe_radius'])

    endTime = time.time()
    print ' in %3.2f s'%(endTime-startTime)

  def save_receptor_MS(self):
    import AlGDock.IO
    IO_Grid = AlGDock.IO.Grid()
    print 'Writing grid output'
    IO_Grid.write('receptor_MS.dx', \
      {'origin':np.array([0., 0., 0.]), 'spacing':self.kwargs['spacing'], 'counts':self.kwargs['counts'], \
       'vals':self.receptor_MS_grid.flatten()})

  def calc_desolvationGrid(self):
    print 'Calculating and saving desolvation grid'
    startTime = time.time()

    SAS_r = self.kwargs['ligand_atom_radius'] + self.kwargs['probe_radius']
    SAS_sphere_pts = SAS_r*self.unit_sphere_pts

    self.desolvationGrid = calc_desolvationGrid(self.receptor_MS_grid, \
      self.kwargs['spacing'], self.kwargs['counts'], \
      self.receptor_SAS_points, self.crd, \
      SAS_sphere_pts, self.LJ_r2, max(np.sqrt(self.LJ_r2)), \
      self.kwargs['ligand_atom_radius'], \
      self.kwargs['probe_radius'], self.kwargs['integration_cutoff'])

    import AlGDock.IO
    IO_Grid = AlGDock.IO.Grid()
    print 'Writing grid output'
    IO_Grid.write(self.FNs['grid'], \
     {'origin':np.array([0., 0., 0.]), \
      'spacing':self.kwargs['spacing'], \
      'counts':self.kwargs['counts'], \
      'vals':self.desolvationGrid.flatten()})

    endTime = time.time()
    print ' in %3.2f s'%(endTime-startTime)

if __name__ == '__main__':
  AstexDiv_Dir = '/Users/dminh/clusters/CCB/AstexDiv_xtal'
  import argparse
  parser = argparse.ArgumentParser(description='Calculates a desolvation grids')
  parser.add_argument('--prmtop_FN', \
    default = AstexDiv_Dir + '/1-build/1tow/receptor.prmtop', \
    help='Input AMBER PRMTOP file')
  parser.add_argument('--inpcrd_FN', \
    default = AstexDiv_Dir + '/3-grids/1tow/receptor.trans.inpcrd', \
    help='Input coordinates')
  parser.add_argument('--header_FN', \
    default = AstexDiv_Dir + '/3-grids/1tow/header_coarse.dx', \
    help='Input grid header (optional)')
  parser.add_argument('--grid_FN', \
    default='desolv.dx', \
    help='Output for desolvation grid')
  parser.add_argument('--probe_radius', default=1.4, \
    help='Radius of the solvent probe, in A')
  parser.add_argument('--ligand_atom_radius', default=1.4, \
    help='Radius of the ligand atom, in A')
  parser.add_argument('--SAS_points', default=1000, \
    help='Number of points on solvent accessible surface per receptor atom')
  parser.add_argument('--integration_cutoff', default=10,
    help='Numerical integration cutoff, in A')
  parser.add_argument('--spacing', nargs=3, type=float, \
    help='Grid spacing (overrides header)')
  parser.add_argument('--counts', nargs=3, type=int, \
    help='Number of point in each direction (overrides header)')
  parser.add_argument('-f')
  args = parser.parse_args()
  
  self = desolvationGridCalculation(**vars(args))
  self.calc_receptor_SAS_points()
  self.calc_receptor_MS()
  self.calc_desolvationGrid()
