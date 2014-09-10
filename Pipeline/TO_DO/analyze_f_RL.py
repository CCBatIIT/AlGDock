# To be run from the AstexDiv_MD/dock/ directory
# The units are kcal/mol

import os, sys
import glob
import gzip, pickle
import numpy as np

# Constants

ncycles = 10
RT = 8.3144621E-3/4.184*300 # kcal/mol
r_site = 4.0 # Angstroms
G_site = -RT*np.log((4./3.)*np.pi*(r_site**3)/1660.)

min_Bs_for_histograms = 10
min_Bs_for_ROC = 4

# Related to figures

separate_figures = False

# For many labels
# ax_pos = (0.1,0.35,0.85,0.6)
# legend_pos = (0., -0.45, 1., .10)

ax_pos = (0.1,0.30,0.85,0.65)
legend_pos = (-0.1, -0.30, 1.1, .10)

# Select the force fields to analyze
# FF_keys = sorted([key for key in f_RLs[lig][f_RLs[lig].keys()[0]].keys()])
# FF_keys = ['dock6','GBSA_min_Psi','GBSA_mean_Psi','GBSA_MBAR_based']
#  FF_keys = ['dock6', 'grid',
#   'Gas_min_Psi',
#   'Gas_mean_Psi',
#   'Gas_MBAR_based',
#   'GBSA_min_Psi',
#   'GBSA_mean_Psi',
#   'GBSA_MBAR_based']
FF_keys = ['dock6', 'GBSA_min_Psi', 'GBSA_mean_Psi', 'GBSA_MBAR_based']
# FF_keys = ['dock6', 'AGBSA_min_Psi', 'AGBSA_mean_Psi', 'AGBSA_MBAR_based']
FF_keys_legend = [key.replace('MBAR_based','binding PMF').replace('_',' ').replace('Psi','$\Psi$') for key in FF_keys]

import matplotlib.pyplot as plt
# plt.ion()

fignum = 0

######################################
# Load all the binding PMF estimates #
######################################

if os.path.isfile('f_RLs.pkl.gz'):
  F = gzip.open('f_RLs.pkl.gz','r')
  f_RLs = pickle.load(F)
  F.close()
else:
  # Load all the binding PMF estimates after ncycles cycles
  f_RL_FNs = glob.glob('*/*/*/*/0/f_RL.pkl.gz')
  f_RLs = {}
  for f_RL_FN in f_RL_FNs:
    (lig_id,snap_xtal,snap_idn) = f_RL_FN.split('/')[1:-2]
    snap_id = snap_xtal+'_'+snap_idn
    F = gzip.open(f_RL_FN,'r')
    (f_grid, B, dock_equilibrated_cycle, dock_mean_acc) = pickle.load(F)
    F.close()
    # Converts all scores from RT to kcal/mol
    scores = dict([(key,RT*val) for (key,val) in B[-1].iteritems()])
    print '%s has %d cycles'%(f_RL_FN, len(f_grid))
    if len(f_grid)==ncycles:
      # Load the dock6 score
      dock6_FN = os.path.join(os.path.dirname(f_RL_FN),'anchor_and_grow_scored.mol2.gz')
      if os.path.isfile(dock6_FN):
        print 'Loading data from '+dock6_FN
        F = gzip.open(dock6_FN,'r')
        for l in range(20):
          line = F.readline()
          if line.startswith('##########    Grid Score'):
            scores['dock6'] = float(line.split()[-1])
            break
        F.close()
      if 'dock6' not in scores.keys():
        raise Exception('No anchor and grow score found in '+dock6_FN)
      if not lig_id in f_RLs.keys():
        f_RLs[lig_id] = {snap_id:scores}
      else:
        f_RLs[lig_id][snap_id] = scores
  F = gzip.open('f_RLs.pkl.gz','w')
  pickle.dump(f_RLs,F)
  F.close()
  print 'Stored all binding PMF estimates to f_RLs.pkl.gz'
  sys.exit()

# Counts
def count_Bs(f_RLs, cutoff):
  f_RLs_n = dict([(key,val) for (key,val) in f_RLs.iteritems() if len(val)>=cutoff])
  active_f_RLs = dict([(key,val) for (key,val) in f_RLs_n.iteritems() if key.find('_')==-1])
  decoy_f_RLs  = dict([(key,val) for (key,val) in f_RLs_n.iteritems() if key.find('_')!=-1])
  print 'There are %d active, %d decoy, and %d total ligands ' \
        'with at least %d binding PMFs'%(\
    len(active_f_RLs),len(decoy_f_RLs),len(f_RLs_n),cutoff)
  return (f_RLs_n, active_f_RLs, decoy_f_RLs)

count_Bs(f_RLs, 0)
ligs_for_histograms = [key for (key,val) in count_Bs(f_RLs, min_Bs_for_histograms)[0].iteritems()]
count_Bs(f_RLs, min_Bs_for_ROC)

# Get net charges
nc = {}
ligand_db_FNs = glob.glob('../6-ligands/P00374/*/ligand.db')
for ligand_db_FN in ligand_db_FNs:
  key = ligand_db_FN.split('/')[-2]
  F = open(ligand_db_FN,'r')
  dat = F.read().strip().split('\n')
  F.close()
  for line in dat:
    if line.startswith('amber_charge'):
      nc[key] = np.sum([float(x) for x in line.split("'")[1::2]])
      break

##################################################
# Analysis of binding PMFs for a single snapshot #
##################################################

snapshot_name = '1pdb_s_1'

# Histogram of actives and decoys
for (FF_key,FF_key_legend) in zip(FF_keys,FF_keys_legend):
  hist_scores_FN = 'hist_scores_%s_%s.png'%(FF_key,snapshot_name)
  hist_low_scores_FN = 'hist_low_scores_%s_%s.png'%(FF_key,snapshot_name)
  if not (os.path.isfile(hist_scores_FN) and \
          os.path.isfile(hist_low_scores_FN)):
    active_scores = np.array([val[snapshot_name][FF_key] for (key,val) in f_RLs.iteritems() if (key.find('_')==-1) and (snapshot_name in val.keys())])
    decoy_scores = np.array([val[snapshot_name][FF_key] for (key,val) in f_RLs.iteritems() if (key.find('_')!=-1) and (snapshot_name in val.keys())])
    print 'Histograms based on %d active and %d decoy compounds'%(\
        len(list(active_scores)),len(list(decoy_scores)))
  if not os.path.isfile(hist_scores_FN):
    print 'Plotting '+hist_scores_FN
    plt.figure(fignum)
    if separate_figures:
      fignum += 1
    plt.clf()
    plt.hist([active_scores, decoy_scores])
    plt.title('Docking scores based on '+FF_key_legend)
    plt.xlabel('Score')
    plt.ylabel('Count')
    plt.legend(['Actives','Decoys'], loc='upper right')
    plt.tight_layout()
    plt.savefig(hist_scores_FN)
  if not os.path.isfile(hist_low_scores_FN):
    cutoff = 20
    if len(active_scores>cutoff)==0 and len(decoy_scores>cutoff)==0:
      continue
    print 'Plotting '+hist_low_scores_FN
    plt.figure(fignum)
    if separate_figures:
      fignum += 1
    plt.clf()
    plt.hist([active_scores[active_scores<cutoff], decoy_scores[decoy_scores<cutoff]])
    plt.title('Low docking scores based on '+FF_key_legend)
    plt.xlabel('Score')
    plt.ylabel('Count')
    plt.legend(['Actives','Decoys'], loc='upper right')
    plt.tight_layout()
    plt.savefig(hist_low_scores_FN)

# Nothing with a non-zero net charge has a binding PMF score less than 0
# The GBSA minimum interaction energy score is much
#bad_actives = [key for (key,val) in f_RLs.iteritems() if (key.find('_')==-1) and (val[snapshot_name][FF_key]>50) and (snapshot_name in val.keys())]
#  FF_key = 'GBSA_min_Psi'
#  all_scores = np.array([val[snapshot_name][FF_key] for (key,val) in f_RLs.iteritems() if (snapshot_name in val.keys())])
#  charges = np.array([nc[key] for (key,val) in f_RLs.iteritems() if (snapshot_name in val.keys())])
#  plt.plot(all_scores,charges,'.')
#  plt.show()
# If the net charge is the score, the AUlC is 0.358
# If the absolute net charge is the score, the AUlC is 0.0762
#    active_scores = np.array([val for (key,val) in nc.iteritems() if key.find('_')==-1])
#    decoy_scores = np.array([val for (key,val) in nc.iteritems() if key.find('_')!=-1])

# Receiver operator characteristic (ROC) curve for a single binding PMF
if not os.path.isfile('ROC_B_%s.png'%snapshot_name):
  print 'Plotting ROC_B_%s.png'%snapshot_name
  plt.figure(fignum)
  if separate_figures:
    fignum += 1
  plt.clf()
  plt.axes(ax_pos)

  AUC = {}
  AUlC = {}
  for FF_key in FF_keys:
    active_scores = np.array([val[snapshot_name][FF_key] for (key,val) in f_RLs.iteritems() if key.find('_')==-1])
    decoy_scores = np.array([val[snapshot_name][FF_key] for (key,val) in f_RLs.iteritems() if key.find('_')!=-1])

    nactives = float(len(active_scores))
    ndecoys = float(len(decoy_scores))
    TPR = np.array([sum(active_scores<cutoff)/nactives for cutoff in np.sort(np.concatenate([active_scores, decoy_scores]))])
    FPR = np.array([sum(decoy_scores<cutoff)/ndecoys for cutoff in np.sort(np.concatenate([active_scores, decoy_scores]))])
    AUC[FF_key] = np.trapz(TPR,FPR)

    FPRp = np.insert(FPR[FPR>0],0,0.001)
    TPRp = TPR[-len(FPRp):]
    AUlC[FF_key] = np.trapz(TPRp,np.log(FPRp))/(-np.log(0.001))

    plt.plot(FPR,TPR,'.-')

  plt.title('Receiver Operator Characteristic')
  plt.xlabel('False Positive Rate')
  plt.ylabel('True Positive Rate')
  plt.legend([legend+', AUlC=%.3f'%AUlC[key] \
    for (key,legend) in zip(FF_keys,FF_keys_legend)],
    bbox_to_anchor=legend_pos, loc=3,
    ncol=2, mode="expand", borderaxespad=0.)
  plt.savefig('ROC_B_%s.png'%snapshot_name)

  print 'Area under ROC curve',AUC
  print 'Area under semilog ROC curve',AUlC

###################################################
# Analysis of binding PMFs for multiple snapshots #
###################################################

# Binding PMF histograms
for lig in ligs_for_histograms:
  if not os.path.isfile('hist_bindingPMF_%s.png'%lig):
    print 'Plotting hist_bindingPMF_%s.png'%lig
    plt.figure(fignum)
    if separate_figures:
      fignum += 1
    plt.clf()
    plt.axes(ax_pos)
    for FF_key in FF_keys:
      B = np.array([val[FF_key] for (key,val) in f_RLs[lig].iteritems()])
      plt.hist(B)
    plt.legend(FF_keys_legend,
      bbox_to_anchor=legend_pos, loc=3,
      ncol=2, mode="expand", borderaxespad=0.)
    if lig.find('_')>-1:
      lig_title = lig[:4] + ' decoy'
    else:
      lig_title = lig
    plt.title(lig_title + ' binding PMFs')
    plt.ylabel('Count')
    plt.xlabel('Binding PMFs (kcal/mol)')
    plt.savefig('hist_bindingPMF_%s.png'%lig)

# Free energy convergence
for lig in ligs_for_histograms:
  if not os.path.isfile('convergence_FE_%s.png'%lig):
    print 'Plotting convergence_FE_%s.png'%lig
    plt.figure(fignum)
    if separate_figures:
      fignum += 1
    plt.clf()
    plt.axes(ax_pos)
    for FF_key in ['AGBSA_MBAR_based']:
      B = np.array([val[FF_key] for (key,val) in f_RLs[lig].iteritems()])
      FE_nsnaps = []
      for n in range(1,len(B)):
        B_subsample = [B[np.random.random_integers(0,len(B)-1,n)] \
          for c in range(100)]
        FE_nsnaps.append([(-(np.log(np.mean(np.exp(-B_c+min(B_c))))-min(B_c)) + \
            G_site) for B_c in B_subsample])
      FE_nsnaps = np.array(FE_nsnaps)
      plt.errorbar(range(1,len(B)), np.mean(FE_nsnaps,1), np.std(FE_nsnaps,1))
#    plt.legend(FF_keys_legend,
#      bbox_to_anchor=legend_pos, loc=3,
#      ncol=2, mode="expand", borderaxespad=0.)
    if lig.find('_')>-1:
      lig_title = lig[:4] + ' decoy'
    else:
      lig_title = lig
    plt.title(lig_title)
    plt.ylabel('Free energy estimate (kcal/mol)')
    plt.xlabel('Number of receptor snapshots')
    plt.savefig('convergence_FE_%s.png'%lig)

# Receiver operator characteristic (ROC) curve for binding free energies
if not os.path.isfile('ROC_geq%d.png'%min_Bs_for_ROC):
  print 'Plotting ROC_geq%d.png'%min_Bs_for_ROC
  plt.figure(fignum)
  if separate_figures:
    fignum += 1
  plt.clf()
  plt.axes(ax_pos)

  AUC = {}
  AUlC = {}

  active_lig_Bs = [lig_B for (lig_id,lig_B) in f_RLs.iteritems() if lig_id.find('_')==-1 and len(lig_B.keys())>min_Bs_for_ROC]
  decoy_lig_Bs = [lig_B for (lig_id,lig_B) in f_RLs.iteritems() if lig_id.find('_')!=-1 and len(lig_B.keys())>min_Bs_for_ROC]

  for FF_key in FF_keys:
    active_FF_key_Bs = [np.array([snap_B[FF_key] for (snap_id, snap_B) in active_lig_Bs[c].iteritems()]) for c in range(len(active_lig_Bs))]
    active_scores = [(-(RT*np.log(np.mean(np.exp((-B_c+min(B_c))/RT)))-min(B_c/RT)) + G_site) for B_c in active_FF_key_Bs]

    decoy_FF_key_Bs = [np.array([snap_B[FF_key] for (snap_id, snap_B) in decoy_lig_Bs[c].iteritems()]) for c in range(len(decoy_lig_Bs))]
    decoy_scores = [(-(RT*np.log(np.mean(np.exp((-B_c+min(B_c))/RT)))-min(B_c/RT)) + G_site) for B_c in decoy_FF_key_Bs]

    nactives = float(len(active_scores))
    ndecoys = float(len(decoy_scores))
    TPR = np.array([sum(active_scores<cutoff)/nactives for cutoff in np.sort(np.concatenate([active_scores, decoy_scores]))])
    FPR = np.array([sum(decoy_scores<cutoff)/ndecoys for cutoff in np.sort(np.concatenate([active_scores, decoy_scores]))])
    AUC[FF_key] = np.trapz(TPR,FPR)
    
    FPRp = np.insert(FPR[FPR>0],0,0.001)
    TPRp = TPR[-len(FPRp):]
    AUlC[FF_key] = np.trapz(TPRp,np.log(FPRp))/(-np.log(0.001))
    
    plt.plot(FPR,TPR,'.-')

  plt.title('Receiver Operator Characteristic')
  plt.xlabel('False Positive Rate')
  plt.ylabel('True Positive Rate')
  plt.legend([legend+', AUlC=%.3f'%AUlC[key] for (key,legend) in zip(FF_keys,FF_keys_legend)],
    bbox_to_anchor=legend_pos, loc=3,
    ncol=2, mode="expand", borderaxespad=0.)
  plt.savefig('ROC_geq%d.png'%min_Bs_for_ROC)
