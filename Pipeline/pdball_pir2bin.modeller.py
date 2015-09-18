# Converts pdball.pir to pdball.bin.
# Should be run using modeller from the same directory as pdball.pir

from modeller import *

log.verbose()
env = environ()

sdb = sequence_db(env)
sdb.read(seq_database_file='pdball.pir', seq_database_format='PIR', \
  chains_list='ALL', minmax_db_seq_len=(30, 4000), clean_sequences=True)
sdb.write(seq_database_file='pdball.bin', \
  seq_database_format='BINARY', chains_list='ALL')
