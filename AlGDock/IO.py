import os
import numpy as np
import MMTK

class Grid:
  """
  Class to read and write alchemical grids.
  
  Data is a dictionary with
  spacing - the grid spacing, in Angstroms.
  counts - the number of points in each dimension.
  vals - the values.
  All are numpy arrays.
  """
  def __init__(self):
    pass

  def read(self, FN, multiplier=None):
    """
    Reads a grid in dx or netcdf format
    The multiplier affects the origin and spacing.
    """
    if FN is None:
      raise Exception('File is not defined')
    elif FN.endswith('.dx') or FN.endswith('.dx.gz'):
      data = self._read_dx(FN)
    elif FN.endswith('.nc'):
      data = self._read_nc(FN)
    else:
      raise Exception('File type not supported')
    if multiplier is not None:
      data['origin'] = multiplier*data['origin']
      data['spacing'] = multiplier*data['spacing']
    return data

  def _read_dx(self, FN):
    """
    Reads a grid in dx format
    """
    if FN.endswith('.dx'):
      F = open(FN,'r')
    else:
      import gzip
      F = gzip.open(FN,'r')
        
    # Read the header
    line = F.readline()
    while line.find('object')==-1:
      line = F.readline()
    header = {}
    header['counts'] = [int(x) for x in line.split(' ')[-3:]]
    for name in ['origin','d0','d1','d2']:
      header[name] = [float(x) for x in F.readline().split(' ')[-3:]]
    F.readline()
    header['npts'] = int(F.readline().split(' ')[-3])

    # Test to make sure the grid type is okay.
    # These conditions are not absolultely essential,
    #   but they reduce the number of subtraction operations.
    if not (header['d0'][1]==0 and header['d0'][2]==0 and
            header['d1'][0]==0 and header['d1'][2]==0 and
            header['d2'][0]==0 and header['d2'][1]==0):
      raise Exception('Trilinear grid must be in original basis')
    if not (header['d0'][0]>0 and header['d1'][1]>0 and header['d2'][2]>0):
      raise Exception('Trilinear grid must have positive coordinates')

    # Read the data
    vals = np.ndarray(shape=header['npts'], dtype=float)
    index = 0
    while index<header['npts']:
      line = F.readline()[:-1]
      items = [float(item) for item in line.split()]
      vals[index:index+len(items)] = items
      index = index + len(items)
    F.close()

    data = {
      'origin':np.array(header['origin']), \
      'spacing':np.array([header['d0'][0],header['d1'][1],header['d2'][2]]), \
      'counts':np.array(header['counts']), \
      'vals':vals}
    return data
    
  def _read_nc(self, FN):
    """
    Reads a grid in netcdf format
    """
    from netCDF4 import Dataset
    grid_nc = Dataset(FN,'r')
    data = {}
    for key in list(grid_nc.variables):
      data[key] = np.array(grid_nc.variables[key][:][0][:])
    grid_nc.close()
    return data

  def write(self, FN, data, multiplier=None):
    """
    Writes a grid in dx or netcdf format.
    The multiplier affects the origin and spacing.

    """
    if multiplier is not None:
      data_n = {'origin':multiplier*data['origin'],
                'counts':data['counts'],
                'spacing':multiplier*data['spacing'],
                'vals':data['vals']}
    else:
      data_n = data
    if FN.endswith('.nc'):
      self._write_nc(FN, data_n)
    elif FN.endswith('.dx') or FN.endswith('.dx.gz'):
      self._write_dx(FN, data_n)
    else:
      raise Exception('File type not supported')
  
  def _write_dx(self, FN, data):
    """
    Writes a grid in dx format
    """
    n_points = data['counts'][0]*data['counts'][1]*data['counts'][2]
    if FN.endswith('.dx'):
      F = open(FN,'w')
    else:
      import gzip
      F = gzip.open(FN,'w')
    
    F.write("""object 1 class gridpositions counts {0[0]} {0[1]} {0[2]}
origin {1[0]} {1[1]} {1[2]}
delta {2[0]} 0.0 0.0
delta 0.0 {2[1]} 0.0
delta 0.0 0.0 {2[2]}
object 2 class gridconnections counts {0[0]} {0[1]} {0[2]}
object 3 class array type double rank 0 items {3} data follows
""".format(data['counts'],data['origin'],data['spacing'],n_points))
    
    for start_n in range(0,len(data['vals']),3):
      F.write(' '.join(['%6e'%c for c in data['vals'][start_n:start_n+3]]) + '\n')

    F.write('object 4 class field\n')
    F.write('component "positions" value 1\n')
    F.write('component "connections" value 2\n')
    F.write('component "data" value 3\n')
    F.close()
  
  def _write_nc(self, FN, data):
    """
    Writes a grid in netcdf format
    """
    n_points = data['counts'][0]*data['counts'][1]*data['counts'][2]
    from netCDF4 import Dataset
    grid_nc = Dataset(FN,'w',format='NETCDF4')
    grid_nc.createDimension('one', 1)
    grid_nc.createDimension('n_cartesian', 3)
    grid_nc.createDimension('n_points', n_points)
    grid_nc.createVariable('origin','f8',('one','n_cartesian'))
    grid_nc.createVariable('counts','i8',('one','n_cartesian'))
    grid_nc.createVariable('spacing','f8',('one','n_cartesian'))
    grid_nc.createVariable('vals','f8',('one','n_points'), zlib=True)
    for key in data.keys():
      grid_nc.variables[key][:] = data[key]
    grid_nc.close()

  def truncate(self, in_FN, out_FN, counts, multiplier=None):
    """
    Truncates the grid at the origin and 
    with a limited number of counts per dimension
    
    multiplier is for the values, not the grid scaling
    """
    data_o = self.read(in_FN)
    nyz_o = data_o['counts'][1]*data_o['counts'][2]
    nz_o = data_o['counts'][2]
    
    min_i = int(-data_o['origin'][0]/data_o['spacing'][0])
    min_j = int(-data_o['origin'][1]/data_o['spacing'][1])
    min_k = int(-data_o['origin'][2]/data_o['spacing'][2])

#    vals = np.ndarray(shape=tuple(counts), dtype=float)
#    for i in range(counts[0]):
#      for j in range(counts[1]):
#        for k in range(counts[2]):
#          vals[i,j,k] = data_o['vals'][(i+min_i)*nyz_o + (j+min_j)*nz_o + (k+min_k)]

    vals = np.array(
      [[[data_o['vals'][(i+min_i)*nyz_o + (j+min_j)*nz_o + (k+min_k)]
        for k in range(counts[2])]
          for j in range(counts[1])]
            for i in range(counts[0])])

    if multiplier is not None:
      vals = vals*multiplier
    
    data_n = {'origin':np.array([0., 0., 0.]), \
      'counts':counts, 'spacing':data_o['spacing'], 'vals':vals.flatten()}
    self.write(out_FN,data_n)

class crd:
  """
  Class to read and write AMBER coordinate/restart and trajectory files.
  """
  def __init__(self):
    pass

  def read(self, FN, natoms=None, return_title=False, \
      multiplier=None, trajectory=False):
    """ 
    Reads an AMBER coordinate/restart or trajectory file.
    
    If natoms is not none, then the coordinates will be split 
      into a list of natoms X 3 arrays.
    The coordinates will be multiplied by multiplier.
    The default of 0.1 converts Angstroms into nanometers.
    """
    if not os.path.isfile(FN):
      raise Exception('Coordinate file %s does not exist!'%FN)
    if FN.endswith('.gz'):
      import gzip
      F = gzip.open(FN, 'r')
    else:
      F = open(FN,'r')
    dat = F.read().strip().split('\n')
    F.close()

    title = dat.pop(0) # Title

    if len(dat[0].split())>1:
      # VMD format (does not specify number of atoms)
      crd = []
      for line in dat:
        crd = crd + [float(x) for x in line.split()]
      crd = np.resize(crd,(len(crd)/3,3))
    else:
      # AMBER format
      file_natoms = int(dat.pop(0)) # Number of atoms
      if (natoms is not None) and (file_natoms!=natoms):
        print "Incorrect number of atoms in crd file"
        return np.array([])
      
      if trajectory:
        w = 8   # For mdcrd
      else:
        w = 12  # For inpcrd
      crd = []
      for line in dat:
        crd = crd + [float(line[x:x+w]) for x in range(0,len(line),w)]
      crd = np.resize(crd,(len(crd)/3,3))

    if multiplier is not None:
      crd = multiplier*crd
    if (natoms is not None):
      crd = np.vsplit(crd,crd.shape[0]/natoms)
      print "  read %d configurations from %s"%(len(crd), FN)

    if return_title:
      return (crd, title)
    else:
      return crd

  def write(self, FN, crd, title='', append=False, \
      multiplier=None, trajectory=False):
    """
    Writes an AMBER coordinate/restart or trajectory file
    """
    if (append and os.path.isfile(FN)):
      if FN.endswith('.gz'):
        import gzip
        F = gzip.open(FN,'a')
      else:
        F = open(FN,'a')
    else:
      if os.path.isfile(FN):
        os.rename(FN,FN+'.BAK')
      if FN.endswith('.gz'):
        import gzip
        F = gzip.open(FN,'w')
      else:
        F = open(FN,'w')
      # Write the header
      F.write(title+'\n') # Title
      if not trajectory:
        F.write('%d\n'%crd.shape[0])
  
    if not trajectory:
      flattened = np.vstack(crd).flatten()
      if multiplier is not None:
        flattened = multiplier*flattened
      for n in range(0,len(flattened),6):
        F.write(''.join(['%12.7f'%val for val in flattened[n:n+6]]) + '\n')
    else:
      for c in crd:
        flattened = c.flatten()
        if multiplier is not None:
          flattened = multiplier*flattened
        for n in range(0,len(flattened),10):
          F.write(''.join(['%8.3f'%val for val in flattened[n:n+10]]) + '\n')

    F.close()

class dcd:
  """
  Class to write DCD files
  """
  def __init__(self, molecule, ligand_atom_order=None, \
      receptorConf=None, ligand_first_atom=0):
    self.molecule = molecule
    self.receptorConf = receptorConf
    self.ligand_first_atom = ligand_first_atom
    if ligand_atom_order is None:
      self.ligand_atom_order = range(len(self.molecule.atoms))
    else:
      self.ligand_atom_order = ligand_atom_order
    pass

  def write(self, FN, confs,
      includeLigand=True, includeReceptor=False,
      factor=1.0/MMTK.Units.Ang,
      delta_t=0.1):
    """
    Writes a DCD file for a trajectory.
    If includeReceptor==True, the receptor coordinates are included.
    """
    import MMTK_DCD  # @UnresolvedImport
    from Scientific import N

    if not isinstance(confs,list):
      confs = [confs]
    
    if includeReceptor and (self.receptorConf is None):
      raise Exception("Missing receptor configuration")
    
    n_atoms = 0
    if includeReceptor:
      receptor_x0 = factor*self.receptorConf[:self.ligand_first_atom,0]
      receptor_y0 = factor*self.receptorConf[:self.ligand_first_atom,1]
      receptor_z0 = factor*self.receptorConf[:self.ligand_first_atom,2]
      receptor_x1 = factor*self.receptorConf[self.ligand_first_atom:,0]
      receptor_y1 = factor*self.receptorConf[self.ligand_first_atom:,1]
      receptor_z1 = factor*self.receptorConf[self.ligand_first_atom:,2]
      n_atoms += self.receptorConf.shape[0]
    if includeLigand:
      n_atoms += len(self.molecule.atoms)
    n_snaps = len(confs)

    fd = MMTK_DCD.writeOpenDCD(FN, n_atoms, n_snaps, 1, 1, delta_t)

    if includeReceptor and includeLigand:
      for array in confs:
        array = factor*array
        x = N.concatenate((receptor_x0,N.take(array[:,0],self.ligand_atom_order),receptor_x1)).astype(N.Float16)
        y = N.concatenate((receptor_y0,N.take(array[:,1],self.ligand_atom_order),receptor_y1)).astype(N.Float16)
        z = N.concatenate((receptor_z0,N.take(array[:,2],self.ligand_atom_order),receptor_z1)).astype(N.Float16)
        MMTK_DCD.writeDCDStep(fd, x, y, z)
      MMTK_DCD.writeCloseDCD(fd)
    elif includeLigand:
      for array in confs:
        array = factor*array
        x = N.take(array[:,0], self.ligand_atom_order).astype(N.Float16)
        y = N.take(array[:,1], self.ligand_atom_order).astype(N.Float16)
        z = N.take(array[:,2], self.ligand_atom_order).astype(N.Float16)
        MMTK_DCD.writeDCDStep(fd, x, y, z)
      MMTK_DCD.writeCloseDCD(fd)
    else:
      x = N.concatenate((receptor_x0,receptor_x1)).astype(N.Float16)
      y = N.concatenate((receptor_y0,receptor_y1)).astype(N.Float16)
      z = N.concatenate((receptor_z0,receptor_z1)).astype(N.Float16)
      MMTK_DCD.writeDCDStep(fd, x, y, z)
      MMTK_DCD.writeCloseDCD(fd)

class prmtop:
  """
  Class to read AMBER prmtop files
  """
  def __init__(self):
    pass

  def read(self, FN, varnames=['RESIDUE_LABEL','RESIDUE_POINTER']):
    """ 
    Reads an AMBER prmtop file, returning a dictionary
    """
    if not os.path.isfile(FN):
      raise Exception('prmtop file %s does not exist!'%FN)
    if FN.endswith('.gz'):
      import gzip
      F = gzip.open(FN, 'r')
    else:
      F = open(FN,'r')
    data = F.read().split('%FLAG ')
    F.close()
    
    prmtop = {}
    for record in data:
      name = record[:record.find('\n')].strip()
      if name in varnames:
        prmtop[name] = self._load_record(record)
    return prmtop

  def _load_record(self, record):
    items = []
    lines = record.split('\n')
    lines.pop(0) # Name
    FORMAT = lines.pop(0).strip()[8:-1] # Format
    if FORMAT.find('a')>-1: # Text
      w = int(FORMAT[FORMAT.find('a')+1:])
      for line in lines:
        items = items + [line[x:x+w] for x in range(0,len(line),w)]
      return np.array(items)
    elif FORMAT.find('I')>-1: # Integer
      w = int(FORMAT[FORMAT.find('I')+1:])
      for line in lines:
        items = items + [int(line[x:x+w]) for x in range(0,len(line),w)]
      return np.array(items, dtype=int)
    elif FORMAT.find('E')>-1: # Scientific
      w = int(FORMAT[FORMAT.find('E')+1:FORMAT.find('.')])
      for line in lines:
        items = items + [float(line[x:x+w]) for x in range(0,len(line),w)]
      return np.array(items, dtype=float)
