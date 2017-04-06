# Expands a SMILES string to a 3D structure in mol2 format with am1bcc charges

try:
  import argparse
  parser = argparse.ArgumentParser()
  parser.add_argument('out_FN', default='AAA.mol2', help='Output file name')
  parser.add_argument('SMILES', default=None, help='SMILES string')
  parser.add_argument('UseOpenEye', choices=['Y','N'], \
    help='Use OpenEye toolkit?')
  args = parser.parse_args()
except ImportError:
  import sys
  class args:
    out_FN = sys.argv[1]
    SMILES = sys.argv[2]
    UseOpenEye = sys.argv[3]

print 'Generating %s based on %s'%(args.out_FN, args.SMILES)

import os, inspect
dirs = {}
dirs['script'] = os.path.dirname(os.path.abspath(\
  inspect.getfile(inspect.currentframe())))
execfile(os.path.join(dirs['script'],'_external_paths.py'))
command_paths = findPaths(['balloon','chimera'])
threeD_FN = os.path.join('3D', os.path.basename(args.out_FN))

for dirN in [os.path.dirname(args.out_FN),'3D']:
  if (dirN!='') and not os.path.isdir(dirN):
    os.system('mkdir -p '+dirN)

if args.UseOpenEye=='Y':
  from openeye import oechem
  from openeye import oequacpac

  # Create a OpenEye molecule object from the SMILES string
  mol = oechem.OEGraphMol()
  if not oechem.OESmilesToMol(mol, args.SMILES):
    raise Exception('Invalid SMILES string', args.SMILES)
    
  oechem.OECanonicalOrderAtoms(mol)
  oechem.OECanonicalOrderBonds(mol)

  # Assign a reasonable protomer
  if not oequacpac.OEGetReasonableProtomer(mol):
    raise Exception('Failed to get a reasonable protomer at pH 7.4')
    
  oechem.OEAssignAromaticFlags(mol, oechem.OEAroModelOpenEye)
  oechem.OEAddExplicitHydrogens(mol)

  smi = oechem.OECreateSmiString(mol, oechem.OESMILESFlag_Canonical)
  print 'The canonical SMILES for a reasonably protonated state is', smi

  # Generate conformations
  from openeye import oeomega

  mol_multiconf = oechem.OEMol(mol)

  omega = oeomega.OEOmega()
  # These parameters were chosen to match http://docs.eyesopen.com/toolkits/cookbook/python/modeling/am1-bcc.html
  omega.SetMaxConfs(800)
  omega.SetIncludeInput(False)
  omega.SetCanonOrder(False)

  omega.SetStrictStereo(True)
  omega.SetStrictAtomTypes(True)

  omega.SetSampleHydrogens(True)  # Word to the wise: skipping this step can lead to significantly different charges!
  omega.SetEnergyWindow(15.0)
  omega.SetRMSThreshold(1.0)  # Word to the wise: skipping this step can lead to significantly different charges!

  omega(mol_multiconf)  # generate conformation

  # Generate am1bcc partial charges
  oequacpac.OEAssignPartialCharges(mol_multiconf, oequacpac.OECharges_AM1BCCSym)

  # Get total charge
  conf = mol_multiconf.GetConf(oechem.OEHasConfIdx(0))
  absFCharge = 0
  sumFCharge = 0
  sumPCharge = 0.0
  for atm in mol_multiconf.GetAtoms():
      sumFCharge += atm.GetFormalCharge()
      absFCharge += abs(atm.GetFormalCharge())
      sumPCharge += atm.GetPartialCharge()
  oechem.OEThrow.Info("%s: %d formal charges give total charge %d ; Sum of Partial Charges %5.4f"
                           % (mol_multiconf.GetTitle(), absFCharge, sumFCharge, sumPCharge))

  # Output file
  ofs = oechem.oemolostream(args.out_FN)
  ofs.SetFormat(oechem.OEFormat_MOL2H)
  oechem.OEWriteMolecule(ofs, conf)

if args.UseOpenEye=='N':
  if not (os.path.isfile(threeD_FN) or \
          os.path.isfile(args.out_FN)):
      # Run Balloon to convert from a SMILES string to a 3D structure
      MMFF94_FN = os.path.join(os.path.dirname(command_paths['balloon']),'MMFF94.mff')
      command = command_paths['balloon'] + ' -f ' + MMFF94_FN + \
        ' --nconfs 1 --nGenerations 300 "' + args.SMILES + '" ' + threeD_FN
      os.system(command)

  if os.path.isfile(os.path.basename(threeD_FN)[:-5]+'_bad.mol2'):
    print 'Conformer generation failed!'
    # Make the final output and empty file
    open(args.out_FN, 'a').close()

  if os.path.isfile(threeD_FN):
    # Select the first model from the mol2 file
    F = open(threeD_FN, "r")
    mol2 = F.read()
    F.close()

    if mol2.count("@<TRIPOS>MOLECULE")>1:
      print 'Keeping first configuration in'+threeD_FN
      confs = mol2.strip().split("@<TRIPOS>MOLECULE")
      if confs[0]=='':
        confs.pop(0)
      F = open(threeD_FN,"w")
      F.write("@<TRIPOS>MOLECULE"+confs[0])
      F.close()

  if not os.path.isfile(args.out_FN):
    # Run chimera to get AM1BCC charges
    prep_script = os.path.join(dirs['script'], '_prep_ligand.chimera.py')
    command = command_paths['chimera'] + " --nogui --script" + \
      " '%s --in_FN %s --out_FN %s'"%(prep_script, threeD_FN, args.out_FN)
    print command
    os.system(command)
  if os.path.isfile(threeD_FN) and \
     os.path.isfile(args.out_FN):
    os.remove(threeD_FN)
