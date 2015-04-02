import AlGDock.BindingPMF_plots
from AlGDock.BindingPMF import *

self = AlGDock.BindingPMF_plots.BPMF_plots()

#import numpy as np
#
#for k in range(self._cool_cycle):
#  confs = self.confs['cool']['samples'][-1][k]
#  for c in range(len(confs)):
#    self.universe.setConfiguration(Configuration(self.universe,confs[c]))
#    self.universe.normalizeConfiguration()
#    self.confs['cool']['samples'][-1][k][c] = np.copy(self.universe.configuration().array)

import itertools
confs = [self.confs['cool']['samples'][-1][k] for k in range(self._cool_cycle)]
confs = np.array([conf[self.molecule.heavy_atoms,:] for conf in itertools.chain.from_iterable(confs)])

from pyRMSD.matrixHandler import MatrixHandler
rmsd_matrix = MatrixHandler().createMatrix(confs,'QCP_SERIAL_CALCULATOR')

# NOSUP_SERIAL_CALCULATOR

#GBSA_energy = [(self.cool_Es[-1][k]['LNAMD_GBSA'][:,-1]-self.cool_Es[-1][k]['LNAMD_Gas'][:,-1]) for k in range(self._cool_cycle)]
#GBSA_energy = np.array(list(itertools.chain.from_iterable(GBSA_energy)))

cum_Nk = np.cumsum([len(self.confs['cool']['samples'][-1][k]) for k in range(self._cool_cycle)])

#  # Compute distance matrix with centering
#  self._write_traj('cool.dcd',confs,moiety='L')
#  import mdtraj as md
#  traj = md.load('cool.dcd',top=self._FNs['prmtop']['L'])
#  dist_matrix = [mdtraj.rmsd(traj,traj,frame=k,atom_indices=traj.topology.select('type!=H')) for k in range(N)]
#  dist_matrix = np.array(dist_matrix)

#class RodriguezLaio:
#  def __init__(self, molecule, confs):
#
#  def _rmsd_pair(conf1, conf2):
#    return np.array([])

# Compute distance matrix without centering
import numpy as np
N = len(confs)
# This does not work in the scipy.cluster.Z() function for some reason
#dist_matrix = [0. if j>=k else np.sqrt(np.sum((confs[j][self.molecule.heavy_atoms,:]-confs[k][self.molecule.heavy_atoms,:])**2)/self.molecule.nhatoms) for j in range(N) for k in range(N)]
#dist_matrix = np.reshape(np.array(dist_matrix), (N,N))
#dist_matrix += np.transpose(dist_matrix)

#flat_dist_matrix = np.array([np.sqrt(np.sum((confs[j]-confs[k])**2)/self.molecule.nhatoms) for j in range(N) for k in range(j+1,N)])
# Hierarchical Clustering
import scipy.cluster
Z = scipy.cluster.hierarchy.complete(rmsd_matrix.get_data())
#
#nclusters = []
#rmv = []
#for dcutoff in [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7]:
#  assignments = np.array(scipy.cluster.hierarchy.fcluster(Z, dcutoff, criterion='distance'))
#  nc = len(set(assignments))
#  GBSA_cluster_variances = []
#  for n in range(nc):
#    GBSA_n = GBSA_energy[assignments==n]
#    if len(GBSA_n)>1:
#      GBSA_cluster_variances.append(np.var(GBSA_n))
#  if len(GBSA_cluster_variances)>0:
#    rmv.append(np.sqrt(np.mean(GBSA_cluster_variances)))
#  else:
#    rmv.append(0.)
#  nclusters.append(nc)

assignments = np.array(scipy.cluster.hierarchy.fcluster(Z, 0.2, criterion='distance'))
cum_Nclusters = [len(set(assignments[:cum_Nk[k]])) for k in range(self._cool_cycle)]

# Density-based clustering
# rho = sum(dist_matrix<0.1,0) # Parzen window estimate of density

# k nearest neighbors density estimate
#  sdist_matrix = np.sort(dist_matrix)
#  k = max(10,int(N/10.))
#  rho = 1/sdist_matrix[:,k]

#  # Find the minimum distance to point with higher density
#  delta = []
#  for j in range(N):
#    distances = [dist_matrix[j,k] for k in range(N) if rho[k]>rho[j]]
#    if len(distances)>0:
#      delta.append(min(distances))
#    else:
#      delta.append(np.max(dist_matrix))
#
#  import matplotlib.pyplot as plt
#  plt.ion()
#  plt.plot(rho,delta,'.')
