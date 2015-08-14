# To be run from the [TARGET]/ directory

# Parameters
from collections import OrderedDict
show = {'FF':OrderedDict([('dock6','DOCK 6'),
              ('MBAR','MBAR (none)'),
              ('NAMD_GBSA_min_Psi','Minimum $\Psi$ (GBSA)'),
              ('NAMD_GBSA_mean_Psi','Mean $\Psi$ (GBSA)'),
              ('NAMD_GBSA_inverse_FEP','Exp Mean $\Psi$ (GBSA)'),
              ('NAMD_GBSA_MBAR','MBAR (GBSA)')]),
        'lib':{'active':'Active', 'decoy':'Decoy'}}
ax_pos = (0.1,0.25,0.85,0.7)
legend_pos = (-0.1, -0.3, 1.1, .10)

# Constants
RT = 8.3144621E-3/4.184*300 # kcal/mol

# Parse input arguments
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--dock6', default='dock6/', \
  help='Directory with dock6 results')
parser.add_argument('--AlGDock', default='AlGDock/dock/', \
  help='Directory with AlGDock results')
parser.add_argument('--ligand', default='ligand/', \
  help='Directory with ligand collections')
parser.add_argument('--analysis', default='analysis/', \
  help='Directory with analysis results')
args = parser.parse_args()
del argparse

# Check for the existence of directories
import os
if not os.path.isdir(args.analysis):
  os.system('mkdir -p '+args.analysis)

# Determine the library names
import glob
library_keys = [os.path.basename(FN)[:-4] \
  for FN in glob.glob(os.path.join(args.ligand,'*.ism'))]

#######################
# Load all the scores #
#######################

import numpy as np
import pickle, gzip

# Extract and sort dock6 scores.
scores = {'dock6':{}}
dock6_scores_FN = os.path.join(args.analysis, 'dock6_scores.pkl.gz')
if os.path.isfile(dock6_scores_FN):
  F = gzip.open(dock6_scores_FN,'r')
  scores['dock6'] = pickle.load(F)
  F.close()
else:
  print 'Parsing dock6 scores...'
  nc_FNs = glob.glob(os.path.join(args.dock6,'*/*/*.nc'))
  nc_FNs = [FN for FN in nc_FNs if os.path.getsize(FN)>0]
  for nc_FN in nc_FNs:
    (dock6_dir, library, key, receptor) = nc_FN.split('/')
    library = '.'.join(library.split('.')[:-1])
    receptor = receptor[:-3]
    key = '%s.%s-%s'%(library,key,receptor)
    from netCDF4 import Dataset
    F = Dataset(nc_FN,'r')
    scores['dock6'][key] =  F.variables['Grid Score'][0]
    F.close()
  F = gzip.open(dock6_scores_FN,'w')
  pickle.dump(scores['dock6'],F)
  F.close()

scores_FN = os.path.join(args.analysis, 'scores.pkl.gz')
if os.path.isfile(scores_FN):
  F = gzip.open(scores_FN,'r')
  scores = pickle.load(F)
  F.close()
else:
  # Extract AlGDock scores
  import pickle
  f_RL_FNs = [FN for FN in glob.glob(os.path.join(\
    args.AlGDock,'*','*','*','f_RL.pkl.gz')) if os.path.getsize(FN)>0]
  for f_RL_FN in f_RL_FNs:
    print 'Reading AlGDock scores from '+f_RL_FN
    (lib_subdir,key,receptor_rep) = os.path.dirname(\
      f_RL_FN[len(args.AlGDock):]).split('/')
    library = '.'.join(lib_subdir.split('.')[:-1])
    receptor = '-'.join(receptor_rep.split('-')[:-1])
    key = '%s.%s-%s'%(library,key,receptor)
    F = gzip.open(f_RL_FN)
    (f_L, stats_RL, f_RL, B) = pickle.load(F)
    F.close()
    if 'grid_MBAR' in scores.keys():
      scores['grid_MBAR'][key] = f_RL['grid_MBAR'][-1][-1]
    else:
      scores['grid_MBAR'] = {key:f_RL['grid_MBAR'][-1][-1]}
    for FF in B.keys():
      if FF in scores.keys():
        scores[FF][key] = B[FF][-1]
      else:
        scores[FF] = {key:B[FF][-1]}
  F = gzip.open(scores_FN,'w')
  pickle.dump(scores,F)
  F.close()

for FF in show['FF'].keys():
  if not FF in scores.keys():
    scores[FF] = {}

# Counts
print '\nNumber of scores'
for FF in scores.keys():
  print '%30s: %d'%(FF, len(scores[FF]))

# Collections are groups of ligands that can be contrasted in ROC plots
collections = {}
# Create a collection of keys in each force field
for FF in show['FF'].keys():
  collections[FF] = set(scores[FF].keys())
# Create collections based on active, inactive, or decoy in the library name
for library_type in ['active','inactive','decoy']:
  collections[library_type] = \
    set([key for key in scores['dock6'].keys() if key.find(library_type)>-1])
# TODO: Read in other collections

# Histograms for each force field
print '\nPlotting histograms:'
import matplotlib.pyplot as plt

for FF in show['FF'].keys():
  for (range_name, cutoff) in [('all',float('inf')),('low',0)]:
    hist_FN = os.path.join(args.analysis,'hist_%s-%s.png'%(range_name,FF))
    hist_scores = []
    hist_legend = []
    if not os.path.isfile(hist_FN):
      for lib in show['lib'].keys():
        scores_FF = [scores[FF][key] for key in scores[FF].keys() \
          if key in collections[lib] and scores[FF][key]<cutoff]
        if len(scores_FF)>0:
          hist_scores.append(scores_FF)
          hist_legend.append(show['lib'][lib])
      if len(hist_scores)>0:
        print '  writing to '+hist_FN
        plt.clf()
        plt.hist(hist_scores)
        plt.title('Histogram of %s %s scores'%(range_name,show['FF'][FF]))
        plt.xlabel('Score')
        plt.ylabel('Count')
        plt.legend(hist_legend)
        plt.savefig(hist_FN)

import sys
def tee(string,F):
  sys.stdout.write(string)
  F.write(string)

# Receiver operating characteristic (ROC) curve
# contrasting positives and negatives
# for each non-empty set in show['FF'].keys()

# TODO: Allow ROC to contrast other collections
positives = 'active'
negatives = 'decoy'

nonempty_FF = [FF for FF in show['FF'].keys() if len(collections[FF])>0]
width_FF_field = np.max([len(show['FF'][FF]) for FF in nonempty_FF])

positive_keys = collections[positives]
negative_keys = collections[negatives]
for FF in nonempty_FF:
  positive_keys = positive_keys.intersection(collections[FF])
  negative_keys = negative_keys.intersection(collections[FF])
if len(positive_keys)==0:
  raise Exception('No positive ligands!')
if len(negative_keys)==0:
  raise Exception('No positive ligands!')

ROC_F = open(os.path.join(args.analysis,'ROC.txt'),'w')
tee('Receiver operating characteristic for %d %s versus %d %s ligands\n'%(\
  len(positive_keys), positives, len(negative_keys), negatives), ROC_F)
tee('%s\t     AUC\t    AUlC\n'%(' '*width_FF_field), ROC_F)
plt.clf()
plt.axes(ax_pos)

AUC = {}
AUlC = {}
ROC_legend = []
for (FF,symbol) in zip(nonempty_FF,('s','d','^','v','>','<','.')):
  sorted_scores_and_status = sorted(\
    [(scores[FF][key],True) for key in positive_keys] + \
    [(scores[FF][key],False) for key in negative_keys])
  isPositive = np.array([val[1] for val in sorted_scores_and_status])
  TPR = np.cumsum(isPositive)/float(sum(isPositive))
  isNegative = np.logical_not(isPositive)
  FPR = np.cumsum(isNegative)/float(sum(isNegative))

  # Area under curve
  AUC[FF] = np.trapz(TPR,FPR)
  if not (FPR[0]>0):
    FPRp = np.insert(FPR[FPR>0],0,0.001)
    TPRp = TPR[-len(FPRp):]
  else:
    FPRp = FPR[FPR>0]
    TPRp = TPR[FPR>0]
  AUlC[FF] = np.trapz(TPRp,np.log(FPRp))/(-np.log(0.001))

  tee('%s: %f\t%f\n'%(show['FF'][FF].rjust(width_FF_field), \
    AUC[FF], AUlC[FF]),ROC_F)

  plt.plot(FPR,TPR,symbol+'-')
  ROC_legend.append(show['FF'][FF]+', AUlC=%.3f'%AUlC[FF])

ROC_F.close()

plt.title('Receiver Operator Characteristic')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.legend(ROC_legend, bbox_to_anchor=legend_pos, loc=3,
  ncol=2, mode="expand", borderaxespad=0.)
plt.savefig(os.path.join(args.analysis,'ROC_%s_v_%s.png'%(positives,negatives)))
