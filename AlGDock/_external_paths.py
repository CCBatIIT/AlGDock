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

def loadModules(programs):
  import os
  for program in programs:
    if program in modules.keys() and \
        os.path.isfile('/share/apps/modules/'+modules[program]):
      os.system('module load '+modules[program])

# Define search paths for external programs and files
# Defined for
# David's IIT MacBook Pro, DSCR cluster, and CCB cluster
search_paths = {
# These files/programs are used in BindingPMF.py
             # Generalized AMBER force field
  'gaff.dat':['/Users/dminh/Installers/AlGDock-0.0.1/data/gaff.dat',
              '/home/dminh/Installers/AlGDock-0.0.1/Data/gaff.dat'],
             # For postprocessing snapshots
      'namd':['/Users/dminh/Applications/NAMD_2.10/namd2',
              '/share/apps/namd/2.9/Linux-x86_64-g++/namd2',
              'namd2'],
             # For postprocessing snapshots
    'sander':['/Users/dminh/Installers/amber14/bin/sander',
              '/share/apps/amber/14/bin/sander',
              'sander'],
             # For postprocessing snapshots
      'apbs':['/Users/dminh/Applications/APBS-1.4.1/APBS.app/Contents/MacOS/apbs',
              '/share/apps/apbs/1.4/bin/apbs',
              'apbs'],
    'ambpdb':['/Users/dminh/Installers/amber14/bin/ambpdb',
              '/share/apps/amber/14/bin/ambpdb',
              'ambpdb'],
   'molsurf':['/Users/dminh/Installers/amber14/bin/molsurf',
              '/share/apps/amber/14/bin/molsurf',
              'molsurf'],
             # BindingPMF.py is built on MMTK
      'MMTK':['/Users/dminh/Installers/MMTK-2.7.9',
              '/home/dminh/Installers/MMTK-2.7.9'],
             # For visualizing (not essential)
       'vmd':['/Applications/VMD.app/Contents/Resources/VMD.app/Contents/MacOS/VMD',
              '/share/apps/vmd/1.9.1/bin/vmd'],
             # For making movies (not essential)
   'convert':['/Applications/ImageMagick-6.9.1/bin/convert',
              '/share/apps/imagemagick/6.9.1-6/bin/convert'],
             # For labeling images (not essential)
      'font':['/Library/Fonts/Microsoft/Arial.ttf',
              '/home/dminh/shared/TTF/Arial.ttf']}

download_paths = {
  'namd':('namd.tar.gz','','namd2'),
  'sander':('sander.tar.gz','','sander'),
  'ambpdb':('ambpdb.tar.gz','','ambpdb'),
  'molsurf':('molsurf.tar.gz', \
    'export LD_LIBRARY_PATH=./molsurf:$LD_LIBRARY_PATH', \
    'molsurf/molsurf'),
  'apbs':('APBS-1.4-linux-static-x86_64.tar.gz','','APBS-1.4-linux-static-x86_64/bin/apbs')}

modules = {'namd':'namd/2.9', \
  'sander':'ambertools/14', \
  'ambpdb':'ambertools/14', \
  'openmm':'openmm/6.2'}