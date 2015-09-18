import os, inspect
dir_external_paths = os.path.dirname(os.path.abspath(\
  inspect.getfile(inspect.currentframe())))

def findPath(locations):
  """
  Parses a list of locations, returning the first file that exists.
  If none exist, then None is returned.
  """
  import os.path
  for location in locations:
    if location is not None and os.path.exists(location):
      return os.path.abspath(location)
  return None

def findPaths(keys):
  paths = dict([(key,findPath(search_paths[key])) \
    for key in keys])
  for key in paths.keys():
    if paths[key] is None:
      # Download file and install program if available
      if key in download_paths.keys():
        (FN,command,path) = download_paths[key]
        # Check that it has not already been downloaded
        import os
        if os.path.isfile(path):
          paths[key] = os.path.abspath(path)
        else:
          import time
          download_start_time = time.time()
          print 'Downloading and installing '+key
          os.system('wget --no-verbose --no-check-certificate http://stash.osgconnect.net/+daveminh/%s'%(FN))
          os.system('tar xzf %s'%FN)
          if command != '':
            os.system(command)
          if os.path.isfile(path):
            print '  ' + key + ' downloaded and installed in %f s'%(\
              time.time() - download_start_time)
            paths[key] = os.path.abspath(path)
          else:
            print 'Could not download '+key
            raise Exception('Could not download '+key)
      else:
        raise Exception('Missing file for '+key)
  return paths

# Define search paths for external programs and files
# Defined for
# David's IIT MacBook Pro, WH210, and the CCB cluster
search_paths = {
# These files/programs are used in the pipeline
     'balloon':['/Users/dminh/Installers/Balloon-1.5.0.1143/balloon',
                '/Applications/Darwin_64bit-Balloon/balloon',
                '/share/apps/balloon/1.5.0.1143/balloon'],
               # For adding hydrogens and charges to receptors
     'pdb2pqr':['/Users/dminh/Applications/pdb2pqr-osx-bin-1.9.0/pdb2pqr',
                '/share/apps/pdb2pqr/1.9.0/pdb2pqr'],
               # For adding hydrogens and charges to ligands
     'chimera':['/Applications/Chimera.app/Contents/MacOS/chimera', # same in WH210
                '/share/apps/chimera/1.9/bin/chimera'],
               # For initial ligand pose
       'dock6':['/Users/dminh/Installers/dock6/bin/dock6',
                '/Applications/dock6/dock6',
                '/share/apps/dock/6/bin/dock6'],
               # For homology modelling
    'modeller':['/Library/modeller-9.15/bin/mod9.15', # same in WH210
                '/share/apps/modeller/9.13/bin/mod9.14'],
               # Submits a command to the queue
'qsub_command':['/Users/dminh/scripts/qsub_command.py',
                '~/scripts/qsub_command.py',
                '/home/dminh/scripts/qsub_command.py',
                '/home/daveminh/scripts/qsub_command.py'],
               # Spheres for UCSF DOCK 6
  'sphgen_cpp':['/Users/dminh/Applications/sphgen_cpp.1.2/sphgen_cpp',
                '/Applications/dock6/sphgen',
                '/share/apps/sphgen_cpp/1.2/sphgen_cpp'],
               # Preparing the system for AMBER
      'sander':['/Users/dminh/Installers/amber14/bin/sander',
                '/Applications/amber14/bin/sander',
                '/share/apps/amber/14/bin/sander'],
               # Calculating a Poisson-Boltzmann Grid
        'apbs':['/Users/dminh/Applications/APBS-1.4.1/APBS.app/Contents/MacOS/apbs',
                '/share/apps/apbs/1.4/bin/apbs'],
               # Generalized AMBER force field
    'gaff.dat':['/Users/dminh/Installers/AlGDock-0.0.1/Data/gaff.dat',
                '/home/dminh/Installers/AlGDock-0.0.1/Data/gaff.dat',
                '/home/daveminh/algdock_data/gaff.dat'],
               # AlGDock
     'algdock':['/Users/dminh/Library/Enthought/Canopy_64bit/User/lib/python2.7/site-packages/AlGDock/BindingPMF.py',
                '/share/apps/canopy/1.5.0/Canopy_64bit/User/lib/python2.7/site-packages/AlGDock/BindingPMF.py']}

algdock_setup = '''
# Modify paths
echo "
search_paths = {
  'gaff.dat':[None],
  'catdcd':[None],
  'namd':[None],
  'sander':[None],
  'MMTK':['$WORK_DIR/AlGDock/MMTK'],
  'vmd':['$WORK_DIR/vmd/bin/vmd']}
" | cat AlGDock/AlGDock/_external_paths.py - > paths.py
mv paths.py AlGDock/AlGDock/_external_paths.py
export ALGDOCK=$WORK_DIR/AlGDock/BindingPMF
'''

# The file to download, special installation commands, and the final location
download_paths = {
  'balloon':('balloon.tar.gz','','Balloon-1.5.0.1143/balloon'),
  'chimera':('chimera.tar.gz', \
    './chimera-1.9-linux_x86_64.bin < chimera_install.in', \
    'chimera-1.9/bin/chimera'),
  'dock6':('dock6.tar.gz','','dock6/bin/dock6'),
  'apbs':('APBS-1.4-linux-static-x86_64.tar.gz','','APBS-1.4-linux-static-x86_64/bin/apbs'),
  'algdock':('algdock.tar.gz',algdock_setup,'AlGDock/BindingPMF')}
