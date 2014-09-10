#!/usr/bin/python

# This module calculates ele and van der Waals interactions at points on a grid

"""
alchemicalGrids

This module provides a class that calculates ele and Lennard-Jones interactions 
at points on a grid.
"""

import os, time, gzip
import numpy as np

class gridCalculation:
  def __init__(self, \
    prmtop_FN='apo.prmtop', inpcrd_FN=None, pqr_FN=None, \
    header_FN=None, site_FN=None, \
    PB_FN=None, ele_FN=None, LJa_FN=None, LJr_FN=None, \
    spacing=None, counts=None):
  
    ### Parse parameters
    self.FNs = {'prmtop':prmtop_FN, 'inpcrd':inpcrd_FN, 'header':header_FN, \
      'pqr':{True:'receptor.pqr',False:pqr_FN}[pqr_FN is None], \
      'site':{True:'../2-binding_site/measured_binding_site.py', \
              False:site_FN}[site_FN is None], \
      'PB':{True:'apbs.nc',False:PB_FN}[PB_FN is None], \
      'ele':{True:'ele.nc',False:ele_FN}[ele_FN is None], \
      'LJa':{True:'LJa.nc',False:LJa_FN}[LJa_FN is None], \
      'LJr':{True:'LJr.nc',False:LJr_FN}[LJr_FN is None]}
    del prmtop_FN, inpcrd_FN, header_FN, ele_FN, LJa_FN, LJr_FN
  
    # Check that input files are available
    for FN in [self.FNs['prmtop'],self.FNs['inpcrd']]:
      if not os.path.exists(FN):
        raise Exception(FN+' missing!')
  
    # Check that output directories are available
    for FN in [self.FNs['ele'],self.FNs['LJa'],self.FNs['LJr']]:
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
      counts = np.array([int(x) for x in headerLines.pop(0).split(' ')[-3:]])
      for name in ['origin','d0','d1','d2']:
        header[name] = [float(x) for x in headerLines.pop(0).split(' ')[-3:]]
      spacing = np.array([header['d0'][0], header['d1'][1], header['d2'][2]])
      del headerF, headerLines

    # Read binding site parameters
    if (self.FNs['site'] is not None) and os.path.isfile(self.FNs['site']):
      print 'Reading binding site parameters from '+self.FNs['site']
      F = open(self.FNs['site'],'r')
      dat = F.read()
      F.close()
      dat = dat[dat.find('half_edge_length =')+18:]
      dat = dat[:dat.find('\n')]
      half_edge_length = float(dat.strip())
      if (spacing is None):
        print 'Using default spacing of 0.25 A'
        spacing = np.array([0.25, 0.25, 0.25])
      if (counts is None):
        counts = np.array(np.ceil(np.array([ \
          2.*half_edge_length/spacing[0], \
          2.*half_edge_length/spacing[1], \
          2.*half_edge_length/spacing[2]])),dtype=int)

    # Loads coordinates
    import AlGDock.IO
    IO_crd = AlGDock.IO.crd()
    self.crd = IO_crd.read(self.FNs['inpcrd'])

    # Outputs files and parameters
    print '*** Files and parameters ***'
    print 'Input AMBER prmtop      :\t' + self.FNs['prmtop']
    print 'Input AMBER inpcrd      :\t' + self.FNs['inpcrd']
    if self.FNs['header'] is not None:
      print 'Input grid header file  :\t' + self.FNs['header']
    if self.FNs['site'] is not None:
      print 'Input binding site info :\t' + self.FNs['site']
    print 'Output Poisson-Boltzmann:\t' + self.FNs['PB']
    print 'Output electrostatics   :\t' + self.FNs['ele']
    print 'Output LJ attractive    :\t' + self.FNs['LJa']
    print 'Output LJ repulsive     :\t' + self.FNs['LJr']
    print 'Grid spacing            :\t', spacing
    print 'Grid counts             :\t', counts
    print

    if not os.path.isfile(self.FNs['PB']):
      print 'Calculating Poisson-Boltzmann grid'
      self.PB_grid(spacing*counts)
    else:
      print 'Poisson-Boltzmann grid already calculated'
    
    if not (os.path.isfile(self.FNs['ele']) and \
            os.path.isfile(self.FNs['LJa']) and \
            os.path.isfile(self.FNs['LJr'])):
      print 'Calculating direct alchemical grids'
      self.direct_grids(spacing, counts)
    else:
      print 'Direct alchemical grids already calculated'

  def direct_grids(self, spacing, counts, no_ele=False):
    """
    Calculates direct grids (Lennard Jones and electrostatic)
    """
    
    # Loads a record from AMBER parameter file
    def _loadRecord(record):
      items = []
      lines = record.split('\n')
      lines.pop(0) # Name
      FORMAT = lines.pop(0).strip()[8:-1] # Format
      if FORMAT.find('a')>-1: # Text
        w = int(FORMAT[FORMAT.find('a')+1:])
        for line in lines:
          items = items + [line[x:x+w] for x in range(0,len(line),w)]
      elif FORMAT.find('I')>-1: # Integer
        w = int(FORMAT[FORMAT.find('I')+1:])
        for line in lines:
          items = items + [int(line[x:x+w]) for x in range(0,len(line),w)]
      elif FORMAT.find('E')>-1: # Scientific
        w = int(FORMAT[FORMAT.find('E')+1:FORMAT.find('.')])
        for line in lines:
          items = items + [float(line[x:x+w]) for x in range(0,len(line),w)]
      return np.array(items)
      
    ### Loads AMBER parameter file
    prmtopF = open(self.FNs['prmtop'],'r')
    prmtopData = prmtopF.read().split('%FLAG ')
    prmtopF.close()
    del prmtopF

    varnames = ['POINTERS','CHARGE','NONBONDED_PARM_INDEX',
      'LENNARD_JONES_ACOEF','LENNARD_JONES_BCOEF','ATOM_TYPE_INDEX']

    prmtop = {}
    for record in prmtopData:
      name = record[:record.find('\n')].strip()
      if name in varnames:
        prmtop[name] = _loadRecord(record)
    del name, record, varnames, prmtopData

    prmtop['CHARGE'] = prmtop['CHARGE']/18.2223 # Convert to units of electric charge

    NATOM = prmtop['POINTERS'][0]
    NTYPES = prmtop['POINTERS'][1]

    ### Extract Lennard-Jones well depth and radii for each atom
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
    # More useful for later calculations
    root_LJ_depth = np.sqrt(LJ_depth)
    LJ_diameter = LJ_radius*2
    del i, LJ_index, factor

    ### Coordinates of grid points
    print 'Calculating grid coordinates'
    startTime = time.time()

    grid = {}
    grid['x'] = np.zeros(shape=tuple(counts), dtype=float)
    grid['y'] = np.zeros(shape=tuple(counts), dtype=float)
    grid['z'] = np.zeros(shape=tuple(counts), dtype=float)
    for i in range(counts[0]):
      for j in range(counts[1]):
        for k in range(counts[2]):
          grid['x'][i,j,k] = i*spacing[0]
          grid['y'][i,j,k] = j*spacing[1]
          grid['z'][i,j,k] = k*spacing[2]

    endTime = time.time()
    print ' in %3.2f s'%(endTime-startTime)

### Calculate ele and Lennard-Jones potential energies at grid points
# Units: kcal/mol A e

# For ele potential:
# E = 1/(4*pi*eps_o) q1 q2 / r
#
#                       1 kg m^2  | (1.60217646E-19 e)^2 | cal    s^2   |   1 kcal | 1E10 A |
# --------------------------------|----------------------|--------------|----------|--------| * 6.0221415E+23
# 4 pi 8.854187817629E-12 C^2 s^2 |                  C^2 | 4.184 kg m^2 | 1000 cal |  1 m   |

# Prefactor is:
# 1/(4*math.pi*8.85418781762E-12)*(1.60217646E-19**2)/4.184*1E10*6.0221415E+23/1000 = 332.06 kcal/mol A e^2

    print 'Calculating grid potential energies'
    startTime = time.time()

    if not no_ele:
      grid['ele'] = np.zeros(shape=tuple(counts), dtype=float)
    grid['LJr'] = np.zeros(shape=tuple(counts), dtype=float)
    grid['LJa'] = np.zeros(shape=tuple(counts), dtype=float)

    for atom_index in range(NATOM):
      dif_x = grid['x'] - self.crd[atom_index][0]
      R2 =  dif_x*dif_x
      del dif_x
      dif_y = grid['y'] - self.crd[atom_index][1]
      R2 += dif_y*dif_y
      del dif_y
      dif_z = grid['z'] - self.crd[atom_index][2]
      R2 += dif_z*dif_z
      del dif_z  
      R = np.sqrt(R2)
      del R2

      atom_type = prmtop['ATOM_TYPE_INDEX'][atom_index]-1
      
      if not no_ele:
        grid['ele'] = grid['ele'] + 332.06*prmtop['CHARGE'][atom_index]/R
      grid['LJr'] += root_LJ_depth[atom_type]*(LJ_diameter[atom_type]**6)/R**12
      grid['LJa'] += -2*root_LJ_depth[atom_type]*(LJ_diameter[atom_type]**3)/R**6

      if atom_index%100==0:
        endTime = time.time()
        print 'Completed atom %d / %d in a total of %3.2f s'%(atom_index,NATOM,endTime-startTime)

    endTime = time.time()
    print '\t%3.2f s'%(endTime-startTime)
  
    # Cap Lennard-Jones potential energies
    u_max = 10000.0
    grid['LJr'] = u_max*np.tanh(grid['LJr']/u_max)
    grid['LJa'] = u_max*np.tanh(grid['LJa']/u_max)

    ### Output grids
    import AlGDock.IO
    IO_Grid = AlGDock.IO.Grid()
    print 'Writing grid output'
    if not no_ele:
      IO_Grid.write(self.FNs['ele'], \
        {'origin':np.array([0., 0., 0.]), 'spacing':spacing, 'counts':counts, 'vals':grid['ele'].flatten()})
    IO_Grid.write(self.FNs['LJr'], \
      {'origin':np.array([0., 0., 0.]), 'spacing':spacing, 'counts':counts, 'vals':grid['LJr'].flatten()})
    IO_Grid.write(self.FNs['LJa'], \
      {'origin':np.array([0., 0., 0.]), 'spacing':spacing, 'counts':counts, 'vals':grid['LJa'].flatten()})

  def PB_grid(self, edge_length):
    """
    Calculates a Poisson-Boltzmann grid using APBS
    
    edge_length is a 3 X 1 numpy array
    """
    import inspect
    import _external_paths
    dirs = {}

    # Sets up pqr file
    if not os.path.exists(self.FNs['pqr']):
      command_paths = _external_paths.findPaths(['sander'])
      dirs['amber'] = os.path.abspath(\
        os.path.dirname(command_paths['sander'])[:-4])

      command = 'cat {0} | {1}/bin/ambpdb -p {2} -pqr > {3}'.format(\
        self.FNs['inpcrd'],dirs['amber'],self.FNs['prmtop'],self.FNs['pqr'])
      os.system(command)

    # Determine the grid parameters
    full_spacing = 1.0
    focus_spacing = 0.5
    final_spacing = focus_spacing

    #   The final grid spans the same space as the other grids
    final_dims = np.array(np.ceil(edge_length/focus_spacing),dtype=int)
    final_center = edge_length/2.

    def roundUpDime(x):
      return (np.ceil((x.astype(float)-1)/32)*32+1).astype(int)

    #   The focus grid has the same center but the dimensions are rounded up
    focus_dims = roundUpDime(final_dims)
    focus_center = final_center

    #   The full grid spans 1.5 times the molecule range
    #                             and the focus grid, whatever is larger
    min_xyz = np.array([min(self.crd[a,:]) for a in range(3)])
    max_xyz = np.array([max(self.crd[a,:]) for a in range(3)])
    mol_range = max_xyz - min_xyz
    mol_center = (min_xyz + max_xyz)/2.

    full_min = np.minimum(mol_center - mol_range/2.*1.5, \
                          focus_center - focus_dims*focus_spacing/2.*1.5)
    full_max = np.maximum(mol_center + mol_range/2.*1.5, \
                          focus_center + focus_dims*focus_spacing/2.*1.5)
    full_dims = roundUpDime((full_max-full_min)/full_spacing)
    full_center = (full_min + full_max)/2.
    
    print 'There are the grid ranges:'
    print 'Full'
    print '  Min',full_min
    print '  Max',full_max
    print '  Center',full_center
    print '  Spacing',full_spacing
    print '  Points per dimension',full_dims
    print 'Focus'
    print '  Min',focus_center - (focus_dims-1)*focus_spacing/2.
    print '  Max',focus_center + (focus_dims-1)*focus_spacing/2.
    print '  Center',focus_center
    print '  Spacing',focus_spacing
    print '  Points per dimension',focus_dims
    print 'Final'
    print '  Min',final_center - (final_dims-1)*final_spacing/2.
    print '  Max',final_center + (final_dims-1)*final_spacing/2.
    print '  Center',final_center
    print '  Spacing',final_spacing
    print '  Points per dimension',final_dims
    
    # Writes APBS script
    apbsF = open('apbs.in','w')
    apbsF.write('''READ
  mol pqr {0}
END
ELEC mg-manual # large grid centered at center of range
  bcfl mdh # multiple debye-huckel boundary condition
  chgm spl4 # quintic B-spline charge discretization
  dime {1[0]} {1[1]} {1[2]}
  gcent {2[0]} {2[1]} {2[2]}
  grid {3} {3} {3}
  lpbe # Linearized Poisson-Boltzmann
  mol 1
  pdie 2.0
  sdens 10.0
  sdie 80.0
  srad 1.4
  srfm smol # Smoothed dielectric and ion-accessibility coefficients
  swin 0.3
  temp 300.0
END
ELEC mg-manual # focus grid around ligand binding site
  bcfl focus # multiple debye-huckel boundary condition
  chgm spl4 # quintic B-spline charge discretization
  dime {4[0]} {4[1]} {4[2]}
  gcent {5[0]} {5[1]} {5[2]}
  grid {6} {6} {6}
  lpbe # Linearized Poisson-Boltzmann
  mol 1
  pdie 2.0
  sdens 10.0
  sdie 80.0
  srad 1.4
  srfm smol # Smoothed dielectric and ion-accessibility coefficients
  swin 0.3
  temp 300.0
  write pot dx apbs_focus
END'''.format(self.FNs['pqr'], \
              full_dims, full_center, full_spacing, \
              focus_dims, focus_center, focus_spacing))
    apbsF.close()

    # Execute APBS
    if not (os.path.isfile('apbs_focus.dx') or os.path.isfile(self.FNs['PB'])):
      try:
        command_paths = _external_paths.findPaths(['apbs'])
        os.system(command_paths['apbs'] + ' apbs.in > apbs.out')
      except:
        print 'APBS failure!'
        return
        
    # Truncate the grid and convert to kcal/mol
    # APBS reports electrostatic grid potential energies in kBT e_c^{-1}
    # The others are in kcal/mol e_c^{-1}
    # At 300 K, 1 kBT ~ 0.596 kcal/mol
    if not os.path.isfile(self.FNs['PB']):
      import AlGDock.IO
      IO_Grid = AlGDock.IO.Grid()
      print final_dims
      IO_Grid.truncate('apbs_focus.dx', self.FNs['PB'], \
        final_dims, multiplier=0.596)

    # Remove intermediate files
    for FN in [self.FNs['pqr'], 'io.mc', 'apbs.in', 'apbs.out']:
      if os.path.isfile(FN):
        os.remove(FN)

if __name__ == '__main__':
  import sys
  
  try:
    import argparse
    parser = argparse.ArgumentParser(description='Calculate van der Waals and ele grids')
    parser.add_argument('--prmtop_FN', help='Input AMBER PRMTOP file')
    parser.add_argument('--inpcrd_FN', help='Input coordinates')
    parser.add_argument('--pqr_FN', help='Input for APBS (optional)')
    parser.add_argument('--header_FN', help='Input grid header (optional)')
    parser.add_argument('--site_FN', help='Input binding site parameters (optional)')
    parser.add_argument('--PB_FN', help='Output for Poisson-Boltzmann grid')
    parser.add_argument('--ele_FN', help='Output for electrostatic grid')
    parser.add_argument('--LJa_FN', help='Output for attractive Lennard-Jones grid')
    parser.add_argument('--LJr_FN', help='Output for repulsive Lennard-Jones grid')
    parser.add_argument('--spacing', nargs=3, type=float, help='Grid spacing (overrides header)')
    parser.add_argument('--counts', nargs=3, type=int, help='Number of point in each direction (overrides header)')
    args = parser.parse_args()
  except:
    import optparse
    parser = optparse.OptionParser()
    parser.add_option('--prmtop_FN', help='Input AMBER PRMTOP file')
    parser.add_option('--inpcrd_FN', help='Input coordinates')
    parser.add_option('--pqr_FN', help='Input for APBS (optional)')
    parser.add_option('--header_FN', help='Input grid header (optional)')
    parser.add_option('--site_FN', help='Input binding site parameters (optional)')
    parser.add_option('--PB_FN', help='Output for Poisson-Boltzmann grid')
    parser.add_option('--ele_FN', help='Output for electrostatic grid')
    parser.add_option('--LJa_FN', help='Output for attractive Lennard-Jones grid')
    parser.add_option('--LJr_FN', help='Output for repulsive Lennard-Jones grid')
    parser.add_option('--spacing', nargs=3, type="float", help='Grid spacing')
    parser.add_option('--counts', nargs=3, type="float", help='Grid dimensions')
    (args,options) = parser.parse_args()

  calc = gridCalculation(**vars(args))
