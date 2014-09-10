import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--tab',help='TAB file name',default='CD38.tab')
args = parser.parse_args()

inF = open(args.tab,'r')
lines = inF.read().strip().split('\n')
inF.close()
lines.pop(0)

def float_affinity(val):
  if val[0]=='>':
    return float("inf")
  elif val.startswith('n/a'):
    return float("nan")
  else:
    return float(val)

cutoff = 100000 # 100 micromolar

activeF = open(args.tab.split()[0]+'.active.ism','w')
inactiveF = open(args.tab.split()[0]+'.inactive.ism','w')

map = {}
nligands = 0
nactive = 0
ninactive = 0
for line in lines:
  fields = line.split('\t')
  if len(fields)>24:
    SMILES = fields[0]
    Kd = float_affinity(fields[24])
    Ki = float_affinity(fields[22])
    IC50 = float_affinity(fields[23])
    EC50 = float_affinity(fields[25])
    output_line = '%s\t%f\t%f\t%f\t%f\n'%(SMILES,Kd,Ki,IC50,EC50)
    if (Kd<cutoff) or (Ki<cutoff) or (IC50<cutoff) or (EC50<cutoff):
      activeF.write(output_line)
      nactive += 1
      ligands += 1
      # Change this for mapping
      # map[nligands] = nactive
    else:
      inactiveF.write(output_line)
      ninactive += 1
      ligands += 1
  else:
    print 'Not enough fields!'
    print line
    print len(fields)

activeF.close()
inactiveF.close()
