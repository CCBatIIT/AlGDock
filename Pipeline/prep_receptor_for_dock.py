# Prepares a PDB file for docking

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('pdb_in', default=None, help='Input PDB file')
args = parser.parse_args()

import os, inspect
dirs = {}
dirs['script'] = os.path.dirname(os.path.abspath(\
  inspect.getfile(inspect.currentframe())))
execfile(os.path.join(dirs['script'],'_external_paths.py'))
command_paths = findPaths(['pdb2pqr','sander','chimera','dock6','sphgen_cpp'])
dirs['amber'] = os.path.abspath(os.path.dirname(command_paths['sander'])[:-4])
dirs['dock6'] = os.path.abspath(os.path.dirname(command_paths['dock6'])[:-4])

if not os.path.isfile(args.pdb_in):
  raise Exception('PDB file not found!')
pdb_path = os.path.dirname(args.pdb_in)
name = os.path.basename(args.pdb_in)[:-4]

surface_script = os.path.join(dirs['script'],'_receptor_surface.chimera.py')
if not os.path.isfile(surface_script):
  raise Exception('Chimera script for dms not found!')

pqr2mol2_script = os.path.join(dirs['script'],'_pqr2mol2.chimera.py')
if not os.path.isfile(pqr2mol2_script):
  raise Exception('Chimera to convert from pqr to mol2 not found!')

for dirN in ['pqr','amber_prep','../amber_in','dock_prep','../dock_in']:
  if not os.path.isdir(dirN):
    os.makedirs(dirN)

# Add hydrogens using pdb2pqr
if not os.path.isfile('pqr/{0}.pdb2pqr_amber.pqr'.format(name)):
  print '*** Adding hydrogens ***'
  os.chdir('pqr')
  command = command_paths['pdb2pqr'] + \
    " ../{0}/{1}.pdb {1}.pdb2pqr_amber.pqr".format(pdb_path, name) + \
    " --ff amber --with-ph=7 --nodebump --ffout=amber"
  os.system(command)
  os.chdir('..')

# Create AMBER prmtop files
if not (os.path.isfile('../amber_in/{0}.prmtop'.format(name)) and \
        os.path.isfile('../amber_in/{0}.inpcrd'.format(name))):
  print '\n*** Preparing AMBER input files ***'
  tleap_F = open('amber_prep/{0}.tleap'.format(name),'w')
  tleap_F.write('''
source leaprc.ff14SB
set default PBRadii bondi

# Receptor
protein = loadpdb pqr/{0}.pdb2pqr_amber.pqr
saveamberparm protein ../amber_in/{0}.prmtop ../amber_in/{0}.inpcrd

savepdb protein amber_prep/{0}.amber.pdb
# savemol2 protein dock_prep/{0}.mol2 0 # default atom types, seems to have a bug
quit
'''.format(name))
  tleap_F.close()
  command = dirs['amber']+'/bin/tleap -f amber_prep/{0}.tleap; ' + \
    'mv leap.log amber_prep/{0}.leaplog'
  command = command.format(name)
  os.system(command)

# Convert prmtop to other formats
if not os.path.isfile('pqr/{0}.bres.pqr'.format(name)):
  print '\n*** Writing standard PQR file (Brookhaven Residue Names) ***'
  command = dirs['amber']+'/bin/ambpdb -p ../amber_in/{0}.prmtop -pqr -bres' + \
    ' < ../amber_in/{0}.inpcrd > pqr/{0}.bres.pqr'
  command = command.format(name)
  os.system(command)

if not os.path.isfile('dock_in/{0}.mol2'.format(name)):
  print '\n*** Writing mol2 file ***'
  command = command_paths['chimera'] + " --nogui --nostatus --script" + \
    " '{0} --pqr_in pqr/{1}.bres.pqr --mol2_out dock_prep/{1}.mol2'"
  command = command.format(pqr2mol2_script, name)
  os.system(command)

if not os.path.isfile('amber_prep/{0}.bres.pdb'.format(name)):
  print '\n*** Writing standard PDB file (Brookhaven Residue Names) ***'
  command = dirs['amber']+'/bin/ambpdb -p ../amber_in/{0}.prmtop -bres' + \
    ' < ../amber_in/{0}.inpcrd > amber_prep/{0}.bres.pdb'
  command = command.format(name)
  os.system(command)

# Create a molecular surface file
if not (os.path.isfile('dock_prep/'+name+'.ms') or \
        os.path.isfile('dock_prep/'+name+'.all.sph') or \
        os.path.isfile('../dock_in/'+name+'.sph')):
  print '\n*** Calculating a molecular surface ***'
  command = command_paths['chimera'] + " --nogui --nostatus --script" + \
    " '{0} --pdb_in {1}/{2}.pdb --dms_out dock_prep/{2}.ms'"
  command = command.format(surface_script, pdb_path, name)
  os.system(command)
  if not os.path.isfile('dock_prep/'+name+'.ms'):
    raise Exception('Surface generation failed!')
else:
  print 'Molecular surface already generated'

# Generate spheres
# Usage: sphgen_cpp -i inputfilename [-s surface_topology]
# [-d surfacepoints] [-l min_distance] [-m min_radius]
# [-x max_radius] -o outputfilename

#  -i	Input file name [required]
#  -o	Output file name [required]
#  -s	R: outside of receptor  L: inside of receptor [default: R]
#  -d	X: all surface points [default: X]
#  -l	Minimum distance between spheres [default: 0.0]
#  -m	Minimum sphere radius [default: 1.4]
#  -x	Maximum sphere radius [default: 4.0]

if not (os.path.isfile('dock_prep/'+name+'.all.sph') or \
        os.path.isfile('../dock_in/'+name+'.sph')):
  print '\n*** Generating all spheres ***'
  command = command_paths['sphgen_cpp'] + \
    " -i dock_prep/{0}.ms -o dock_prep/{0}.all.sph".format(name) + \
    " -s R -d X -l 0.0 -m 1.4 -x 4.0"
  os.system(command)
  if os.path.isfile('dock_prep/'+name+'.all.sph'):
    os.remove('dock_prep/{0}.ms'.format(name))
  else:
    raise Exception('Shere generation failed!')
else:
  print 'Spheres already generated'

# Load the binding site radius and half edge length
if os.path.isfile('../2-binding_site/measured_binding_site.py'):
  execfile('../2-binding_site/measured_binding_site.py')
else:
  raise Exception('No binding site information')
maxR2 = half_edge_length**2

# Keep spheres that are within the half edge length of the grid center
if not os.path.isfile('../dock_in/'+name+'.sph'):
  print '\n*** Selecting spheres ***'
  insphF = open('dock_prep/'+name+'.all.sph','r')
  insph = insphF.read().strip().split('\n')
  insphF.close()

  import numpy as np

  # Keep clusters which have spheres within the center of mass box
  osph = []
  keepCluster = False
  for line in insph:
    if line.startswith('DOCK'):
      pass
    elif line.startswith('cluster'):
      if keepCluster:
        osph += osph_c
      osph_c = []
      keepCluster = False
    elif len(line)>40:
      osph_c.append(line)
      try:
        pos = np.array([float(line[5:15]),float(line[15:25]),float(line[25:35])])
      except:
        print line
        raise Exception('Could not convert line!')
      if ((pos[0]>com_min[0]) and (pos[0]<com_max[0]) and \
          (pos[1]>com_min[1]) and (pos[1]<com_max[1]) and \
          (pos[2]>com_min[2]) and (pos[2]<com_max[2])):
        keepCluster = True

  osphF = open('../dock_in/'+name+'.sph','w')
  osphF.write(insph[0]+'\n')
  osphF.write('cluster     1   number of spheres in cluster   %d\n'%len(osph))
  osphF.write('\n'.join(osph))
  osphF.flush()
  osphF.close()
  if os.path.isfile('../dock_in/'+name+'.sph'):
    pass
    # os.remove('dock_prep/'+name+'.all.sph')
  else:
    raise Exception('Sphere selection failed!')
else:
  print 'Spheres already selected'

# Run showbox
if not os.path.isfile('dock_prep/'+name+'.box.pdb'):
  print '\n*** Generating the box ***'
  showboxF = open('dock_prep/'+name+'.showbox.in','w')
  showboxF.write('''N
U
{0} {0} {0}
{1} {1} {1}
dock_prep/{2}.box.pdb
'''.format(half_edge_length, 2*half_edge_length, name))
  showboxF.close()
  command = '{0}/bin/showbox < dock_prep/{1}.showbox.in'.format(dirs['dock6'], name)
  os.system(command)
  if os.path.isfile('dock_prep/'+name+'.box.pdb'):
    os.remove('dock_prep/{0}.showbox.in'.format(name))
  else:
    raise Exception('Box generation failed!')
else:
  print 'Box already generated'

# Calculate the grid
if not os.path.isfile('../dock_in/'+name+'.nrg'):
  print '\n*** Calculating the DOCK6 grid ***'
  F = open('dock_prep/'+name+'.grid.in','w')
  F.write('''compute_grids                  yes
grid_spacing                   0.25
output_molecule                no
contact_score                  no
energy_score                   yes
energy_cutoff_distance         9999
atom_model                     a
attractive_exponent            6
repulsive_exponent             12
distance_dielectric            yes
dielectric_factor              4
bump_filter                    yes
bump_overlap                   0.75
receptor_file                  dock_prep/{0}.mol2
box_file                       dock_prep/{0}.box.pdb
vdw_definition_file            {1}/parameters/vdw_AMBER_parm99.defn
score_grid_prefix              ../dock_in/{0}
  '''.format(name, dirs['dock6']))
  F.close()
  command = '{0}/bin/grid -i dock_prep/{1}.grid.in'.format(dirs['dock6'], name)
  os.system(command)
  if os.path.isfile('../dock_in/'+name+'.nrg'):
    os.remove('dock_prep/{0}.grid.in'.format(name))
  else:
    raise Exception('Grid calculation failed!')
else:
  print 'Grid already calculated'