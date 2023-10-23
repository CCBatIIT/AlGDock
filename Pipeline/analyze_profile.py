import os, inspect

script_dir = os.path.dirname(os.path.abspath(\
  inspect.getfile(inspect.currentframe())))
exec(open(os.path.join(script_dir,'_load_profile.py')).read())

if not os.path.isdir('figures'):
  os.makedirs('figures')

# Plot histogram of all sequence identities
seq_identities = [p[2] for p in profile.values()]

import matplotlib.pylab as plt
plt.clf()
plt.hist(seq_identities)
plt.xlabel('Sequence Identity')
plt.ylabel('Number of Structures')
plt.savefig('figures/hist_seq_id.png')

# Plot histogram of selected sequences, with
# sequence identities greater than min_seq_identity and
# equivalent positions greater than min_equivalent_positions
selected_seq_identities = \
  [p[2] for p in profile.values() \
    if (p[2]>=min_seq_identity) and (p[1]>=min_equivalent_positions)]
if len(selected_seq_identities)>0:
  plt.clf()
  plt.hist(selected_seq_identities)
  plt.xlabel('Sequence Identity')
  plt.ylabel('Number of Structures')
  plt.savefig('figures/hist_seq_id_selected.png')

print('%d selected chains'%len(selected_seq_identities))
print('with minimum sequence identity of %d'%min_seq_identity)
print('and at least %d equivalent positions'%min_equivalent_positions)

# Sort sequences by
# 1. Sequence Identity
# 2. Equivalent Positions

import numpy as np
seq_identities = np.array([p[2] for p in profile.values()], dtype=float)
equivalent_positions = np.array([p[1] for p in profile.values()], dtype=float)
scores = seq_identities + equivalent_positions/max(equivalent_positions)
inds = np.argsort(scores)
inds = inds[np.logical_and(\
  seq_identities[inds]>=min_seq_identity, \
  equivalent_positions[inds]>=min_equivalent_positions)]

prof_list = profile.items()
print('\n'.join(['{0[0][0]} {0[0][1]} {0[1][2]} {0[1][1]} {0[1][3]}'.format(prof_list[ind]) for ind in inds[::-1]]))
