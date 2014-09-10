# To be run from the [TARGET]/ directory

# Parameters
show = {'FF':{'dock6':'DOCK 6', \
#              'grid':'grid',
#              'Gas_MBAR_based':'MBAR-based (Gas)',
#              'Gas_min_Psi':'Minimum \Psi (Gas)',
#              'GBSA_MBAR_based':'MBAR-based (GBSA)',
#              'GBSA_min_Psi':'Minimum \Psi (GBSA)'},
#              'NAMD_Gas_MBAR_based':'MBAR-based (Gas)',
#              'NAMD_Gas_mean_Psi':'Mean $\Psi$ (Gas)',
#              'NAMD_Gas_min_Psi':'Minimum $\Psi$ (Gas)',
              'NAMD_GBSA_MBAR_based':'BPMF (GBSA)',
              'NAMD_GBSA_mean_Psi':'Mean $\Psi$ (GBSA)',
              'NAMD_GBSA_min_Psi':'Minimum $\Psi$ (GBSA)'},
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
  help='Directory with ligand libraries')
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
library_keys = ['','active','decoy'] + [os.path.basename(FN)[:-4] \
  for FN in glob.glob(os.path.join(args.ligand,'*.ism'))]

#######################
# Load all the scores #
#######################

# Extract and sort dock6 scores.
# Store them in separate libraries and all together.
import numpy as np
import gzip

scores = {'dock6':{}}
SMILES = {}
if not np.array([os.path.isfile(os.path.join(\
                  args.analysis, 'dock6-%s.txt'%lib)) \
                    for lib in library_keys]).all():
  print 'Parsing dock6 scores...'
  mol2_FNs = glob.glob(os.path.join(args.dock6,'*.mol2.gz'))
  mol2_FNs = [FN for FN in mol2_FNs if os.path.getsize(FN)>0]
  
  for mol2_FN in mol2_FNs:
    key = os.path.basename(mol2_FN)[5:-8]
    F = gzip.open(mol2_FN,'r')
    F.readline() # Blank
    SMILES_str_line = F.readline().split()
    if len(SMILES_str_line)==0:
      F.close()
      continue
    SMILES[key] = SMILES_str_line[-1]
    F.readline() # Cluster size
    scores['dock6'][key] = float(F.readline().split()[-1])
    F.close()

  dock6_scores = [(key, SMILES[key], scores['dock6'][key]) \
    for key in scores['dock6'].keys()]
  dock6_scores.sort(key=lambda record: float(record[2]))

  for lib in library_keys:
    F = open(os.path.join(args.analysis, 'dock6-%s.txt'%lib),'w')
    F.write('\n'.join(['{0[0]}\t{0[1]}\t{0[2]}'.format(score) \
                        for score in dock6_scores if score[0].find(lib)!=-1]))
    F.close()
else:
  FN = os.path.join(args.analysis, 'dock6-.txt')
  print 'Reading dock scores from '+FN
  F = open(FN,'r')
  dock6_scores = [line.split('\t') for line in F.read().strip().split('\n')]
  F.close()

  SMILES = dict([(score[0],score[1]) for score in dock6_scores])
  scores['dock6'] = dict([(score[0],float(score[2])) for score in dock6_scores])

# Extract AlGDock scores
import pickle
f_RL_FNs = [FN for FN in glob.glob(os.path.join(\
  args.AlGDock,'*','f_RL.pkl.gz')) if os.path.getsize(FN)>0]
for f_RL_FN in f_RL_FNs:
  key = '-'.join(os.path.dirname(f_RL_FN).split('/')[-1].split('-')[:-1])
  F = gzip.open(f_RL_FN)
  (f_cool, f_L_solvation, Psi, f_grid, B, \
   dock_equilibrated_cycle, dock_mean_acc) = pickle.load(F)
  F.close()
  B = B[-1]
  for FF in B.keys():
    if FF in scores.keys():
      scores[FF][key] = B[FF]
    else:
      scores[FF] = {key:B[FF]}
  del f_cool, f_L_solvation, Psi, f_grid, B
  del dock_equilibrated_cycle, dock_mean_acc

for FF in show['FF'].keys():
  if not FF in scores.keys():
    scores[FF] = {}

# Get net charges
if not os.path.isfile('nc.pkl.gz'):
  print 'Determining net charge...'
  nc = {}
else:
  print 'Loading net charge from nc.pkl.gz'
  nc_F = gzip.open('nc.pkl.gz')
  nc = pickle.load(nc_F)
  nc_F.close()

new_nc = False
ligand_db_FNs = glob.glob('ligand/AlGDock_in/*.db')
for ligand_db_FN in ligand_db_FNs:
  key = os.path.basename(ligand_db_FN)[:-3]
  if not key in nc:
    new_nc = True
    F = open(ligand_db_FN,'r')
    dat = F.read().strip().split('\n')
    F.close()
    for line in dat:
      if line.startswith('amber_charge'):
        nc[key] = int(round(np.sum([float(x) for x in line.split("'")[1::2]])))
        break
if new_nc:
  nc_F = gzip.open('nc.pkl.gz','w')
  pickle.dump(nc, nc_F)
  nc_F.close()

# Counts
print '\nNumber of scores'
for FF in scores.keys():
  print '%30s: %d'%(FF, len(scores[FF]))

# Libraries
libraries = {}
libraries['AM1BCC'] = set([key for key in scores['dock6'].keys() \
  if SMILES[key].find('Gasteiger')==-1])
libraries['Gasteiger'] = set([key for key in scores['dock6'].keys() \
  if not key in libraries['AM1BCC']])

libraries['dock6_AM1BCC'] = libraries['AM1BCC']
libraries['AlGDock_AM1BCC'] = set([key for key in scores['grid'].keys() \
  if SMILES[key].find('Gasteiger')==-1])
for FF in show['FF'].keys():
  libraries[FF] = set(scores[FF].keys())

active = [score[0] for score in dock6_scores \
  if score[0].find('active')!=-1 and score[0] in libraries['AlGDock_AM1BCC']]
decoy = [score[0] for score in dock6_scores \
  if score[0].find('decoy')!=-1 and score[0] in libraries['AlGDock_AM1BCC']]
libraries['active'] = set(active)
libraries['decoy'] = set(decoy)
libraries['top'] = set(active[:50] + decoy[:150])
del dock6_scores, active, decoy
all_AlGDock_ligands = set([key.lower().split('-')[0] \
  for key in scores['grid'].keys()])
for nc_val in set(nc.values()):
  nc_ligands = [key for key in nc.keys() if key in all_AlGDock_ligands and nc[key]== nc_val]
  libraries['nc%d'%nc_val] = set([key for key in libraries['AlGDock_AM1BCC'] if key.split('-')[0].lower() in nc_ligands])
#  show['lib']['nc%d'%nc_val] = 'Net charge %d'%nc_val

# Histograms for each force field
# TO DO: Have a histogram for each library
print '\nPlotting histograms:'
import matplotlib.pyplot as plt

for FF in show['FF'].keys():
  for (range_name, cutoff) in [('AM1BCC',float('inf')),('low',0)]:
    hist_FN = os.path.join(args.analysis,'hist_%s-%s.png'%(range_name,FF))
    hist_scores = []
    hist_legend = []
    if not os.path.isfile(hist_FN):
      for lib in show['lib'].keys():
        scores_FF = [scores[FF][key] for key in scores[FF].keys() \
          if key in libraries[lib] and scores[FF][key]<cutoff]
        if len(scores_FF)>0:
          hist_scores.append(scores_FF)
          hist_legend.append(show['lib'][lib])
      if len(hist_scores)>0:
        print '  writing to '+hist_FN
        plt.clf()
        plt.hist(hist_scores)
        plt.title('Histogram of %s %s scores'%(range_name,show['FF'][FF]))
        plt.xlabel('Score (kcal/mol)')
        plt.ylabel('Count')
        plt.legend(hist_legend)
        plt.savefig(hist_FN)

import sys
def tee(string,F):
  sys.stdout.write(string)
  F.write(string)

# Receiver operating characteristic (ROC) curve for each force field
ROC_F = open(os.path.join(args.analysis,'ROC.txt'),'w')
for libname in ['dock6_AM1BCC','AlGDock_AM1BCC','top','AM1BCC','Gasteiger'] + \
               ['nc%d'%nc_val for nc_val in set(nc.values())]:
  tee('\nReceiver operating characteristic for %s\n'%libname, ROC_F)
  tee('%35s: actives\tdecoys\ttotal\tAUC\tAUlC\n'%' ', ROC_F)

  plt.clf()
  plt.axes(ax_pos)

  AUC = {}
  AUlC = {}
  ROC_legend = []
  for FF in show['FF'].keys():
    if not libraries[libname].issubset(libraries[FF]):
      continue

    sorted_scores_and_status = sorted([(scores[FF][key],key.find('active')!=-1) for key in libraries[libname]])
    
    if len(sorted_scores_and_status)==0:
      continue
    
    isActive = [val[1] for val in sorted_scores_and_status]
    nactive = sum(isActive)
    TPR = np.cumsum(isActive)/float(nactive)
    isDecoy = [(not val[1]) for val in sorted_scores_and_status]
    ndecoy = sum(isDecoy)
    
    if nactive==0 or ndecoy==0:
      tee('%d active and %d decoy molecules\n'%(nactive,ndecoy),ROC_F)
      break
    
    FPR = np.cumsum(isDecoy)/float(ndecoy)
    AUC[FF] = np.trapz(TPR,FPR)
    
    if not (FPR[0]>0):
      FPRp = np.insert(FPR[FPR>0],0,0.001)
      TPRp = TPR[-len(FPRp):]
    else:
      FPRp = FPR[FPR>0]
      TPRp = TPR[FPR>0]
    AUlC[FF] = np.trapz(TPRp,np.log(FPRp))/(-np.log(0.001))

    tee('%35s: %6d\t%6d\t%6d\t%f\t%f\n'%(show['FF'][FF], \
      nactive, ndecoy, len(isActive), \
      AUC[FF], AUlC[FF]),ROC_F)

    plt.plot(FPR,TPR,'.-')
    ROC_legend.append(show['FF'][FF]+', AUlC=%.3f'%AUlC[FF])

  if nactive>0 and ndecoy>0:
    plt.title('Receiver Operator Characteristic')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    if len(ROC_legend)>0:
      plt.legend(ROC_legend, bbox_to_anchor=legend_pos, loc=3,
        ncol=2, mode="expand", borderaxespad=0.)
    plt.savefig(os.path.join(args.analysis,'ROC_%s.png'%libname))

ROC_F.close()
