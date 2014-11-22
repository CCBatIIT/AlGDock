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
        if os.path.isfile(path):
          paths[key] = os.path.abspath(path)
        else:
          import time
          download_start_time = time.time()
          print 'Downloading and installing '+key
          os.system('wget --no-check-certificate http://stash.osgconnect.net/+daveminh/%s'%(FN))
          os.system('tar -xvf %s'%FN)
          if command != '':
            os.system(command)
          if os.path.isfile(path):
            print key + ' downloaded and installed in %d s'%(\
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
# David's IIT MacBook Pro, DSCR cluster, and CCB cluster
search_paths = {
# These files/programs are used in HREX.py
             # Generalized AMBER force field
  'gaff.dat':['/Users/dminh/Installers/AlGDock-0.0.1/data/gaff.dat',
              '/home/dbchem/dm225/.local/installers/AlGDock-0.0.1/Data/gaff.dat',
              '/home/dminh/Installers/AlGDock-0.0.1/Data/gaff.dat'],
             # For postprocessing snapshots
      'namd':['/home/xin/DevelopmentTool/NAMD_2.9_Linux-x86_64-multicore/namd2',
	     							 #'/Users/dminh/Installers/NAMD_2.9_Source/MacOSX-x86_64-g++/namd2',
              '/home/dbchem/dm225/.local/bin/namd2',
              '/share/apps/namd/2.9/Linux-x86_64-g++/namd2'],
             # For postprocessing snapshots
    'sander':['/Users/dminh/Installers/amber14/bin/sander',
              '/home/dbchem/dm225/.local/installers/amber14/bin/sander',
              '/share/apps/amber/14/bin/sander'],
             # HREX.py is built on MMTK
      'MMTK':['/home/xin/DevelopmentTool/MMTK-2.7.9',
              '/home/dbchem/dm225/.local/installers/MMTK-2.7.9',
              '/home/dminh/Installers/MMTK-2.7.9'],
             # For visualizing (not essential)
       'vmd':['/usr/local/bin/vmd'
					#'/Applications/VMD 1.9.1.app/Contents/Resources/VMD.app/Contents/MacOS/VMD',
              '/home/dbchem/dm225/.local/bin/vmd',
              '/share/apps/vmd/1.9.1/bin/vmd']}

download_paths = {
  'namd':('namd.tar.gz','','namd2')}
