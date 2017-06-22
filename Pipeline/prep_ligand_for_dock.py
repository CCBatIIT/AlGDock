# Expands a SMILES string to a 3D structure
# in mol2 format with am1bcc charges and sybyl atom types

try:
  import argparse
  parser = argparse.ArgumentParser()
  parser.add_argument('out_prefix', default='1hnn/ligand', \
    help='Output prefix')
  parser.add_argument('inp', default='1hnn/ligand_in.mol2', \
    help='SMILES string or input file name')
  parser.add_argument('UseOpenEye', choices=['Y','N'], \
    help='Use OpenEye toolkit?')
  parser.add_argument('--RetainProtonation', action='store_true', \
    help='Retains protonation state from input file')
  parser.add_argument('--RetainConformer', action='store_true', \
    help='Retains conformer from input file')
  args = parser.parse_args()
except ImportError:
  import sys
  class args:
    out_prefix = sys.argv[1]
    inp = sys.argv[2]
    UseOpenEye = sys.argv[3]
smi = args.inp

import os, inspect
dirs = {}
dirs['script'] = os.path.dirname(os.path.abspath(\
  inspect.getfile(inspect.currentframe())))
execfile(os.path.join(dirs['script'],'_external_paths.py'))
command_paths = findPaths(['balloon','chimera'])
# Only necessary without OpenEye
balloon_FN = os.path.abspath(args.out_prefix + '_balloon.mol2')
charged_FN = os.path.abspath(args.out_prefix + '_AM1BCC.mol2')
sybyl_FN =   os.path.abspath(args.out_prefix + '_sybyl.mol2')

def step_complete(FN):
  FN = os.path.abspath(FN)
  if os.path.isfile(FN):
    return True
  if FN==charged_FN and os.path.isfile(sybyl_FN):
    return True
  if FN==balloon_FN and os.path.isfile(sybyl_FN):
    return True
  return False

for dirN in [os.path.dirname(args.out_prefix)]:
  if (dirN!='') and not os.path.isdir(dirN):
    os.system('mkdir -p '+dirN)

if args.UseOpenEye=='Y' and not step_complete(charged_FN):
  from openeye import oechem
  from openeye import oequacpac

  mol = oechem.OEGraphMol()
  
  if os.path.isfile(args.inp):
    ifs = oechem.oemolistream(args.inp)
    oechem.OEReadMolecule(ifs, mol)
    ifs.close()
  else:
    # Create a OpenEye molecule object from the SMILES string
    if not oechem.OESmilesToMol(mol, smi):
      raise Exception('Invalid SMILES string', smi)

  oechem.OECanonicalOrderAtoms(mol)
  oechem.OECanonicalOrderBonds(mol)

  # Assign a reasonable protomer
  if args.RetainProtonation:
    for atom in mol.GetAtoms():
      atom.SetImplicitHCount(0)
  else:
    if not oequacpac.OEGetReasonableProtomer(mol):
      print 'Failed to get a reasonable protomer at pH 7.4'

  oechem.OEAssignAromaticFlags(mol, oechem.OEAroModelOpenEye)

  if not args.RetainProtonation:
    oechem.OEAddExplicitHydrogens(mol)

  smi = oechem.OECreateSmiString(mol, oechem.OESMILESFlag_Canonical)
  print 'The canonical SMILES for a reasonably protonated state is', smi

  # Generate conformations
  from openeye import oeomega

  mol_multiconf = oechem.OEMol(mol)
  oechem.OECanonicalOrderAtoms(mol_multiconf)

  omega = oeomega.OEOmega()
  # These parameters were chosen to match http://docs.eyesopen.com/toolkits/cookbook/python/modeling/am1-bcc.html
  omega.SetMaxConfs(800)
  omega.SetIncludeInput(False)
  omega.SetCanonOrder(False)

  omega.SetStrictStereo(False)
  omega.SetStrictAtomTypes(False)

  omega.SetSampleHydrogens(True)  # Word to the wise: skipping this step can lead to significantly different charges!
  omega.SetEnergyWindow(15.0)
  omega.SetRMSThreshold(1.0)  # Word to the wise: skipping this step can lead to significantly different charges!

  if omega(mol_multiconf):  # generate conformation
    # Generate am1bcc partial charges
    oequacpac.OEAssignCharges(mol_multiconf, oequacpac.OEAM1BCCELF10Charges())

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
    ofs = oechem.oemolostream(charged_FN)
    ofs.SetFormat(oechem.OEFormat_MOL2H)
    oechem.OEWriteMolecule(ofs, conf)
    ofs.close()
  else:
    # Conformer generation failed. Use Ballon + Chimera
    print 'Conformer generation with OETools failed.'

if (args.UseOpenEye=='N') or not step_complete(charged_FN):
  if not step_complete(balloon_FN):
      # Run Balloon to convert from a SMILES string to a 3D structure
      MMFF94_FN = os.path.join(os.path.dirname(command_paths['balloon']),'MMFF94.mff')
      command = command_paths['balloon'] + ' -f ' + MMFF94_FN + \
        ' --nconfs 1 --nGenerations 300 "' + smi + '" ' + balloon_FN
      os.system(command)

  if os.path.isfile(os.path.basename(balloon_FN)[:-5]+'_bad.mol2'):
    print 'Conformer generation failed!'
    # Make the final out_prefixut an empty file
    open(sybyl_FN, 'a').close()

  if step_complete(balloon_FN):
    # Select the first model from the mol2 file
    F = open(balloon_FN, "r")
    mol2 = F.read()
    F.close()

    if mol2.count("@<TRIPOS>MOLECULE")>1:
      print 'Keeping first configuration in '+balloon_FN
      confs = mol2.strip().split("@<TRIPOS>MOLECULE")
      if confs[0]=='':
        confs.pop(0)
      F = open(balloon_FN,"w")
      F.write("@<TRIPOS>MOLECULE"+confs[0])
      F.close()

  # Get the net charge based on the SMILES string
  charge = 0
  lc = ''
  for c in smi:
    if c=='+':
      if lc.isdigit():
        charge += int(lc)
      else:
        charge += 1
    elif c=='-':
      if lc.isdigit():
        charge -= int(lc)
      else:
        charge -= 1
    lc = c
  print 'Net charge is ', charge

  if not step_complete(charged_FN):
    # Run chimera to get AM1BCC charges
    prep_script = os.path.join(dirs['script'], '_prep_ligand.chimera.py')
    command = command_paths['chimera'] + " --nogui --script" + \
      " '%s --in_FN %s --out_FN %s --net_charge %d'"%(prep_script, balloon_FN, charged_FN, charge)
    os.system(command)
  if os.path.isfile(balloon_FN) and \
     os.path.isfile(charged_FN):
    os.remove(balloon_FN)

# Restore the original configuration to the charged file
if not step_complete(sybyl_FN):
  if args.RetainConformer:
    from openeye import oechem
    if not os.path.isfile(args.inp):
      raise Exception('File %s not found'%args.inp)
    mol_in = oechem.OEGraphMol()
    ifs = oechem.oemolistream(args.inp)
    oechem.OEReadMolecule(ifs, mol_in)
    ifs.close()
    oechem.OECanonicalOrderAtoms(mol_in)

    if not os.path.isfile(charged_FN):
      raise Exception('File %s not found'%charged_FN)
    mol_out = oechem.OEGraphMol()
    ifs = oechem.oemolistream(charged_FN)
    oechem.OEReadMolecule(ifs, mol_out)
    ifs.close()
    
    if mol_in.GetMaxAtomIdx() != mol_out.GetMaxAtomIdx():
      raise Exception('Number of atoms in input, %d, not equal to number in output, %d'%(\
        mol_in.GetMaxAtomIdx(), mol_out.GetMaxAtomIdx()))
    
    import numpy as np
    coords = np.zeros((mol_in.GetMaxAtomIdx(),3))
    coords_dict = mol_in.GetCoords()
    for a, ind_n in zip(mol_in.GetAtoms(), range(mol_in.GetMaxAtomIdx())):
      coords[ind_n,:] = coords_dict[a.GetIdx()]
    mol_out.SetCoords(coords.flatten())
    
    ofs = oechem.oemolostream(sybyl_FN)
    ofs.SetFormat(oechem.OEFormat_MOL2H)
    oechem.OEWriteMolecule(ofs, mol_out)
    ofs.close()
  else:
    os.system('cp %s %s'%(charged_FN, sybyl_FN))
