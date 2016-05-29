import os, inspect
dir_external_paths = os.path.dirname(os.path.abspath(\
  inspect.getfile(inspect.currentframe())))

def findPath(locations):
  """
  Parses a list of locations, returning the first file that exists.
  If none exist, then None is returned.
  """
  import os.path
  locations = [l for l in locations if l is not None]
  for location in locations:
    if os.path.exists(location):
      return os.path.abspath(location)
  if len(locations)>0:
    print '  failed to find file in: '+', '.join(locations)
  return None

def findPaths(keys):
  """
  Finds a path for each key.
  """
  paths = {}
  for key in keys:
    paths[key] = findPath(search_paths[key]) \
      if key in search_paths.keys() else None
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

def loadModules(programs):
  import os
  import numpy as np
  if os.path.isdir('/share/apps/amber/16') and np.array([p in programs \
      for p in ['sander','elsize','gbnsr6','ambpdb','molsurf']]).any():
    os.environ['PATH'] = '/share/apps/amber/16/bin:' + os.environ['PATH']
    os.environ['LD_LIBRARY_PATH'] = '/share/apps/netcdf/4.3.0/lib:/share/apps/amber/16/lib:' + os.environ['LD_LIBRARY_PATH']
    os.environ['PYTHONPATH'] = '/share/apps/amber/16/lib/python2.6/site-packages' + os.environ['PYTHONPATH']
    os.environ['AMBERHOME'] = '/share/apps/amber/16'

# Define search paths for external programs and files
# Defined for
# David's IIT MacBook Pro, DSCR cluster, and CCB cluster
search_paths = {
# These files/programs are used in BindingPMF.py
             # Generalized AMBER force field
  'gaff.dat':['/Users/dminh/Installers/AlGDock-0.0.1/data/gaff.dat',
              '/home/dminh/Installers/AlGDock-0.0.1/Data/gaff.dat'],
             # For postprocessing snapshots with NAMD
      'namd':['/Users/dminh/Applications/NAMD_2.10/namd2',
              '/share/apps/namd/2.9/Linux-x86_64-g++/namd2',
              'namd2'],
             # For postprocessing snapshots with sander
    'sander':['/Users/dminh/Installers/amber16/bin/sander',
              '/share/apps/amber/16/bin/sander',
              'sander'],
    'elsize':['/Users/dminh/Installers/amber16/bin/elsize',
              '/share/apps/amber/16/bin/elsize',
              'elsize'],
             # For postprocessing snapshots with gbnsr6
    'gbnsr6':['/Users/dminh/Installers/amber16/bin/gbnsr6',
              '/share/apps/amber/16/bin/gbnsr6',
              'gbnsr6'],
             # For postprocessing snapshots with APBS
      'apbs':['/Users/dminh/Applications/APBS-1.4.1/APBS.app/Contents/MacOS/apbs',
              '/share/apps/apbs/1.4/bin/apbs',
              'apbs'],
    'ambpdb':['/Users/dminh/Installers/amber16/bin/ambpdb',
              '/share/apps/amber/16/bin/ambpdb',
              'ambpdb'],
   'molsurf':['/Users/dminh/Installers/amber16/bin/molsurf',
              '/share/apps/amber/16/bin/molsurf',
              'molsurf'],
             # BindingPMF.py is built on MMTK
      'MMTK':['/Users/dminh/Installers/MMTK-2.7.9',
              '/home/dminh/Installers/MMTK-2.7.9'],
             # For visualizing (not essential)
       'vmd':['/Applications/VMD.app/Contents/Resources/VMD.app/Contents/MacOS/VMD',
              '/share/apps/vmd/1.9.1/bin/vmd'],
             # For concatenating and stripping dcd trajectories
    'catdcd':['/Users/dminh/Applications/catdcd',
              '/share/apps/catdcd/4.0/catdcd'],
             # For making movies (not essential)
   'convert':['/Applications/ImageMagick-6.9.1/bin/convert',
              '/share/apps/imagemagick/6.9.1-6/bin/convert'],
             # For labeling images (not essential)
      'font':['/Library/Fonts/Microsoft/Arial.ttf',
              '/home/dminh/shared/TTF/Arial.ttf']}

download_paths = {
  'namd':('namd.tar.gz','','namd2'),
  'sander':('sander.tar.gz','','sander'),
  'elsize':('elsize.tar.gz','','elsize'),
  'gbnsr6':('gbnsr6.tar.gz','','gbnsr6'),
  'ambpdb':('ambpdb.tar.gz','','ambpdb'),
  'molsurf':('molsurf.tar.gz', \
    'export LD_LIBRARY_PATH=./molsurf:$LD_LIBRARY_PATH', \
    'molsurf/molsurf'),
  'apbs':('APBS-1.4-linux-static-x86_64.tar.gz','','APBS-1.4-linux-static-x86_64/bin/apbs')}
