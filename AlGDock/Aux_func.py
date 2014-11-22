import os         #Miscellaneous operating system interfaces
from os.path import abspath
from os.path import exists
from os.path import isfile
from os.path import join
import cPickle as pickle
import gzip
import copy

import time
import numpy as np

import MMTK
import MMTK.Units
from MMTK.ParticleProperties import Configuration

import Scientific
from Scientific_vector import Vector  # @UnresolvedImport

import AlGDock as a
import pymbar.timeseries

########################
# Auxilliary functions #
########################

def random_rotate():
  """
  Return a random rotation matrix
  """
  u = np.random.uniform(size=3)

  # Random quaternion
  q = np.array([np.sqrt(1-u[0])*np.sin(2*np.pi*u[1]),
               np.sqrt(1-u[0])*np.cos(2*np.pi*u[1]),
               np.sqrt(u[0])*np.sin(2*np.pi*u[2]),
               np.sqrt(u[0])*np.cos(2*np.pi*u[2])])

  # Convert the quaternion into a rotation matrix
  rotMat = np.array([[q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3],
                     2*q[1]*q[2] - 2*q[0]*q[3],
                     2*q[1]*q[3] + 2*q[0]*q[2]],
                    [2*q[1]*q[2] + 2*q[0]*q[3],
                     q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3],
                     2*q[2]*q[3] - 2*q[0]*q[1]],
                    [2*q[1]*q[3] - 2*q[0]*q[2],
                     2*q[2]*q[3] + 2*q[0]*q[1],
                     q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3]]])
  return rotMat

#  def choice(N, weights=None):
#    """
#    Weighted random number selection of N values
#    """
#    weights /= np.sum(weights)
#    cumweights = np.cumsum(weights)
#    beginIndicies = np.zeros(N, dtype=int)
#    randVals = np.random.uniform(size=N)
#
#    for (bI_idx,randVal) in zip(range(N),randVals):
#      beginIndicies[bI_idx] = sum(cumweights < randVal)
#    beginIndicies.sort()
#    return beginIndicies

def merge_dictionaries(dicts, required_consistency=[], print_inconsistent=False):
  """
  Merges a list of dictionaries, giving priority to items in descending order.
  Items in the required_consistency list must be consistent with one another.
  """
  merged = {}
  for a in range(len(dicts)): # Loop over all dictionaries,
                              # giving priority to the first
    if not isinstance(dicts[a],dict):
      continue
    for key in dicts[a].keys():
      if dicts[a][key] is None:
        if not key in merged.keys():
          merged[key] = None
      elif isinstance(dicts[a][key],dict):
        merged[key] = merge_dictionaries(
          [dicts[n][key] for n in range(len(dicts)) if key in dicts[n].keys()],
          required_consistency=required_consistency)
      else:
        # Check for consistency with other dictionaries
        for b in (range(a) + range(a+1,len(dicts))):
          if isinstance(dicts[b],dict) and (key in dicts[b].keys()) and (dicts[b][key] is not None):
            if (isinstance(dicts[b][key],np.ndarray)):
              inconsistent_items = (dicts[b][key]!=dicts[a][key]).any()
            else:
              inconsistent_items = (dicts[b][key]!=dicts[a][key])
            if inconsistent_items:
              if print_inconsistent:
                print 'Dictionary items are inconsistent:'
                print dicts[a][key]
                print dicts[b][key]
              if key in required_consistency:
                raise Exception('Items must be consistent!')
        merged[key] = dicts[a][key]
  return merged

def convert_dictionary_relpath(d, relpath_o=None, relpath_n=None):
  """
  Converts all file names in a dictionary from one relative path to another.
  If relpath_o is None, nothing is joined to the original path.
  If relpath_n is None, an absolute path is used.
  """
  converted = {}
  for key in d.keys():
    if d[key] is None:
      pass
    elif isinstance(d[key],dict):
      converted[key] = convert_dictionary_relpath(d[key],
        relpath_o = relpath_o, relpath_n = relpath_n)
    elif isinstance(d[key],str):
      if relpath_o is not None:
        p = abspath(join(relpath_o,d[key]))
      else:
        p = abspath(d[key])
      if exists(p): # Only save file names for existent paths
        if relpath_n is not None:
          converted[key] = os.path.relpath(p,relpath_n)
        else:
          converted[key] = p
  return converted

def HMStime(s):
  """
  Given the time in seconds, an appropriately formatted string.
  """
  if s<60.:
    return '%.3f s'%s
  elif s<3600.:
    return '%d:%.3f'%(int(s/60%60),s%60)
  else:
    return '%d:%d:%.3f'%(int(s/3600),int(s/60%60),s%60)

class NullDevice():
  """
  A device to suppress output
  """
  def write(self, s):
    pass


