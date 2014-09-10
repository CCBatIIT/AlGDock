import os

# Create a list of PDB ids and chains in the profile
profileF = open('profile.prf','r')
profile_lines = profileF.read().split('\n')
profileF.close()

min_seq_identity = 90
min_equivalent_positions = 1
if os.path.isfile('../search_options.py'):
  execfile('../search_options.py') # may override min_seq_identity and other defaults

profile = {}
for line in profile_lines:
  if not(len(line)>110 and line[0]!='#'):
    continue
  split_line = line.split()
  if split_line[0]=='1':
    sequence = split_line[1]
    continue

  # Keys are pdb_id, chain_id
  # 0.
  # The starting position of the target sequence in the alignment.
  # The ending position of the target sequence in the alignment.
  # The starting position of the database sequence in the alignment.
  # The ending position of the database sequence in the alignment.
  # 1. The number of equivalent positions in the alignment.
  # 2. The sequence identity between the target and database sequence.
  # 3. The sequence alignment.
  profile[(split_line[1][:-1], split_line[1][-1])] = \
    ((int(split_line[5]), int(split_line[6]),
      int(split_line[7]), int(split_line[8])), \
    int(split_line[9]), int(split_line[10][:-1]), split_line[-1])