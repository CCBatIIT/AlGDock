import sys
inFN = sys.argv[-1]

if inFN.endswith('.mol2'):
  outFN = inFN[:-5]+'.nc'
elif inFN.endswith('.mol2.gz'):
  outFN = inFN[:-8]+'.nc'

import os
if os.path.isfile(outFN):
  sys.exit()

import AlGDock.IO
IO_dock6_mol2 = AlGDock.IO.dock6_mol2()
(confs, Es) = IO_dock6_mol2.read(inFN)

if confs==[]:
  F = open(inFN,'w')
  F.close()
  sys.exit()

import numpy as np
confs = np.array(confs)
confs = np.array(confs)/10. # Convert Angstroms to nanometers

from netCDF4 import Dataset
dock6_nc = Dataset(outFN,'w',format='NETCDF4')
dock6_nc.createDimension('n_poses', confs.shape[0])
dock6_nc.createDimension('n_atoms', confs.shape[1])
dock6_nc.createDimension('n_cartesian', confs.shape[2])
dock6_nc.createDimension('one',1)
dock6_nc.createVariable('confs','f4',('n_poses','n_atoms','n_cartesian'), \
  zlib=True, complevel=9, shuffle=True)
dock6_nc.variables['confs'][:,:,:] = confs
for key in Es.keys():
  datatype = 'i2' if key=='Cluster Size' else 'f4'
  dock6_nc.createVariable(key, datatype,('n_poses'), \
    zlib=True, complevel=9, shuffle=True)
  dock6_nc.variables[key][:] = np.array(Es[key])
dock6_nc.close()

os.remove(inFN)
