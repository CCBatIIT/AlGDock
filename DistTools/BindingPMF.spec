# -*- mode: python -*-

# Creates an executable for AlGDock
# pyinstaller AlGDocker.spec

import os, sys, glob
import Scientific, MMTK, AlGDock, simtk

extension_library_files = glob.glob(os.path.join(os.path.dirname(Scientific.__file__), sys.platform, '*'))
extension_library_files += glob.glob(os.path.join(os.path.dirname(MMTK.__file__), sys.platform, '*'))
extension_library_files += glob.glob(os.path.join(os.path.dirname(AlGDock.__file__), sys.platform, '*')) 

lib_paths = []
for path in sys.path:
  lib_ind = path.find('lib')
  if lib_ind>0:
    lib_path = path[:lib_ind]+'lib'
    if not lib_path in lib_paths:
      lib_paths.append(lib_path)

for path in lib_paths:
  extension_library_files += glob.glob(os.path.join(path,'*netcdf*'))
  extension_library_files += glob.glob(os.path.join(path,'*hdf5*'))
  extension_library_files += glob.glob(os.path.join(path,'*libmkl*'))
  extension_library_files += glob.glob(os.path.join(path,'*libiomp5*'))

print 'Extension Libraries'
print extension_library_files
extension_libraries = [(os.path.basename(F),F,'BINARY') for F in extension_library_files]

init_script_files = glob.glob(os.path.join(os.path.dirname(AlGDock.__file__), '__p*.py'))
init_scripts = [(os.path.join('AlGDock',os.path.basename(F)),F,'DATA') for F in init_script_files]

data_files = []
code_extensions = ['pyc','so']
for pkg in ['Scientific','MMTK','AlGDock','simtk']:
  pkgdir = os.path.dirname(locals()[pkg].__file__)
  for cdir, subdirs, files in os.walk(pkgdir):
    rp = os.path.relpath(cdir,os.path.join(pkgdir,'..'))
    data_files += [(os.path.join(rp,F),os.path.join(cdir,F),'DATA') for F in files if F.split('.')[-1] not in code_extensions]

del os, sys, glob, Scientific, MMTK, AlGDock

a = Analysis(['../AlGDock/BindingPMF.py'],
             pathex=['../AlGDock'],
             hiddenimports=['scipy.special._ufuncs_cxx','netCDF4_utils','netcdftime'],
             hookspath=None,
             runtime_hooks=None)

a.binaries += TOC(extension_libraries)

# Removes data files that cause security alerts
newdatas = []
for ind in reversed(range(len(a.datas))):
  if a.datas[ind][0].startswith('..'):
    olddata = a.datas.pop(ind)
    print 'Removing %s from a.datas'%(olddata[0])
a.datas += TOC(init_scripts)
a.datas += TOC(data_files)

pyz = PYZ(a.pure)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='BindingPMF',
          debug=False,
          strip=None,
          upx=False,
          console=True )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=None,
               upx=False,
               name='AlGDock')
