import optparse
parser = optparse.OptionParser()
parser.add_option('--pqr_in', default=None, help='Input PQR file')
parser.add_option('--mol2_out', default=None, help='Output mol2 file')
(args, options) = parser.parse_args()

import os
if not os.path.exists(args.pqr_in):
  raise Exception(args.pqr_in+' does not exist')

print 'Converting '+args.pqr_in+' to '+args.mol2_out

import chimera
chimera.openModels.open(args.pqr_in)

# Molecular surface
chimera.runCommand("addh")
chimera.runCommand("addcharge")
chimera.runCommand("write format mol2 0 "+args.mol2_out)
