import optparse
parser = optparse.OptionParser()
parser.add_option('--pdb_in', default=None, help='Input PDB file')
parser.add_option('--dms_out', default=None, help='Output dms file')
(args, options) = parser.parse_args()

import os
if not os.path.exists(args.pdb_in):
  raise Exception(args.pdb_in+' does not exist')

print 'Processing '+args.pdb_in

import chimera
chimera.openModels.open(args.pdb_in)

# Molecular surface
from chimera import runCommand, openModels, MSMSModel
chimera.runCommand("surf") # generate surface using 'surf' command
surf = openModels.list(modelTypes=[MSMSModel])[0] # get the surf object
from WriteDMS import writeDMS
writeDMS(surf, args.dms_out) # write DMS
