#!/usr/bin/env python

package_name = "AlGDock"

from distutils.core import setup, Command, Extension
from distutils.command.build import build
from distutils.command.sdist import sdist
from distutils.command.install_data import install_data
from distutils import dir_util
from distutils.filelist import FileList, translate_pattern
import distutils.sysconfig
sysconfig = distutils.sysconfig.get_config_vars()

import os, sys, types
import ctypes, ctypes.util
from glob import glob

class Dummy:
    pass
pkginfo = Dummy()
execfile('AlGDock/__pkginfo__.py', pkginfo.__dict__)
execfile('AlGDock/_external_paths.py')

from site import USER_BASE as userbase

# Check module dependencies

# Check for Cython
try:
  from Cython.Distutils import build_ext
  cython_ok = True  
except ImportError:
  cython_ok = False
if not cython_ok:
  print 'AlGDock requires Cython'
  raise SystemExit

# Check that we have Scientific 2.6 or higher
try:
    from Scientific import __version__ as scientific_version
    if scientific_version[-2:] == 'hg':
        scientific_version = scientific_version[:-2]
    scientific_version = scientific_version.split('.')
    scientific_ok = int(scientific_version[0]) >= 2 and \
                    int(scientific_version[1]) >= 6
except ImportError:
    scientific_ok = False
if not scientific_ok:
    print "AlGDock needs ScientificPython 2.6 or higher"
    raise SystemExit

# Check that we have MMTK 2.6 or higher
try:
  from MMTK import __version__ as mmtk_version
  mmtk_version = mmtk_version.split('.')
  mmtk_ok = int(mmtk_version[0]) >= 2 and \
            int(mmtk_version[1]) >= 6
except ImportError:
  mmtk_ok = False
if not mmtk_ok:
  print "AlGDock requires MMTK version 2.6 or higher"
  raise SystemExit

# Configure compile arguments and include directories 
compile_args = []
include_dirs = [os.path.join(findPath(search_paths['MMTK']),'Include'),
                os.path.join(findPath(search_paths['MMTK']),'include',
                             'python%d.%d'%sys.version_info[:2],'MMTK')]

from Scientific import N
assert N.package == "NumPy"

compile_args.append("-DNUMPY=1")
import numpy.distutils.misc_util
include_dirs.extend(numpy.distutils.misc_util.get_numpy_include_dirs())

if (int(scientific_version[1]) >= 8 or \
    (int(scientific_version[1]) == 7 and int(scientific_version[2]) >= 8)):
    netcdf_h = os.path.join(sys.prefix, 'include',
                            'python%d.%d' % sys.version_info[:2],
                            'Scientific', 'netcdf.h')
    if os.path.exists(netcdf_h):
        compile_args.append("-DUSE_NETCDF_H_FROM_SCIENTIFIC=1")
        include_dirs.append(os.path.join(sys.prefix, 'include',
                            'python%d.%d' % sys.version_info[:2]))

    userbase_netcdf_h = os.path.join(userbase, 'include',
                            'python%d.%d' % sys.version_info[:2],
                            'Scientific', 'netcdf.h')
    if os.path.exists(userbase_netcdf_h):
        compile_args.append("-DUSE_NETCDF_H_FROM_SCIENTIFIC=1")
        include_dirs.append(os.path.join(sys.prefix, 'include',
                            'python%d.%d' % sys.version_info[:2]))
else:
    # Take care of the common problem that netcdf is in /usr/local but
    # /usr/local/include is not on $CPATH.
    if os.path.exists('/usr/local/include/netcdf.h'):
        include_dirs.append('/usr/local/include')
        # netcdf is in /opt/local/include/netcdf.h
        # include_dirs.append('/opt/local/include')
    if ('NETCDF_PREFIX' in os.environ):
      include_dirs.append(os.path.join(os.environ['NETCDF_PREFIX'],'include'))

for user_include_path in [
    os.path.join(userbase, 'include', 
       'python%d.%d' % sys.version_info[:2])]:
  if os.path.exists(user_include_path):
    include_dirs.append(user_include_path)

headers = []
paths = [os.path.join('AlGDock', 'ForceFields', 'Cylinder'),
         os.path.join('AlGDock', 'ForceFields', 'Sphere'),
         os.path.join('AlGDock', 'ForceFields', 'Grid'),
         os.path.join('AlGDock', 'Integrators', 'VelocityVerlet'),
         os.path.join('AlGDock', 'Integrators', 'NUTS')]
data_files = []
for dir in paths:
    files = []
    for f in glob(os.path.join(dir, '*')):
        if f[-3:] != '.py' and f[-4:-1] != '.py' and os.path.isfile(f):
            files.append(f)
    data_files.append((dir, files))

class ModifiedFileList(FileList):

    def findall(self, dir=os.curdir):
        from stat import ST_MODE, S_ISREG, S_ISDIR, S_ISLNK
        list = []
        stack = [dir]
        pop = stack.pop
        push = stack.append
        while stack:
            dir = pop()
            names = os.listdir(dir)
            for name in names:
                if dir != os.curdir:
                    fullname = os.path.join(dir, name)
                else:
                    fullname = name
                stat = os.stat(fullname)
                mode = stat[ST_MODE]
                if S_ISREG(mode):
                    list.append(fullname)
                elif S_ISDIR(mode) and not S_ISLNK(mode):
                    list.append(fullname)
                    push(fullname)
        self.allfiles = list


class modified_build(build):

    def has_sphinx(self):
        if sphinx is None:
            return False
        setup_dir = os.path.dirname(os.path.abspath(__file__))
        return os.path.isdir(os.path.join(setup_dir, 'Doc'))

    sub_commands = build.sub_commands + [('build_sphinx', has_sphinx)]

class modified_sdist(sdist):

    def run (self):

        self.filelist = ModifiedFileList()
        self.check_metadata()
        self.get_file_list()
        if self.manifest_only:
            return
        self.make_distribution()

    def make_release_tree (self, base_dir, files):
        self.mkpath(base_dir)
        dir_util.create_tree(base_dir, files,
                             verbose=self.verbose, dry_run=self.dry_run)
        if hasattr(os, 'link'):         # can make hard links on this system
            link = 'hard'
            msg = "making hard links in %s..." % base_dir
        else:                           # nope, have to copy
            link = None
            msg = "copying files to %s..." % base_dir
        if not files:
            self.warn("no files to distribute -- empty manifest?")
        else:
            self.announce(msg)
        for file in files:
            if os.path.isfile(file):
                dest = os.path.join(base_dir, file)
                self.copy_file(file, dest, link=link)
            elif os.path.isdir(file):
                dir_util.mkpath(os.path.join(base_dir, file))
            else:
                self.warn("'%s' not a regular file or directory -- skipping"
                          % file)

class modified_install_data(install_data):

    def run(self):
        install_cmd = self.get_finalized_command('install')
        self.install_dir = getattr(install_cmd, 'install_lib')
        return install_data.run(self)

class test(Command):

    user_options = []
    def initialize_options(self):
        self.build_lib = None
    def finalize_options(self):
        self.set_undefined_options('build',
                                   ('build_lib', 'build_lib'))

    def run(self):
        import sys, subprocess
        self.run_command('build_py')
        self.run_command('build_ext')
        ff = sum((fns for dir, fns in data_files if 'ForceFields' in dir), [])
        for fn in ff:
            self.copy_file(fn,
                           os.path.join(self.build_lib, fn),
                           preserve_mode=False)
#        subprocess.call([sys.executable, 'Tests/all_tests.py'],
#                        env={'PYTHONPATH': self.build_lib,
#                             'MMTKDATABASE': 'MMTK/Database'})

cmdclass = {
    'build' : modified_build,
    'sdist': modified_sdist,
    'install_data': modified_install_data,
    'build_ext': build_ext,
    'test': test
}

# Build the sphinx documentation if Sphinx is available
try:
    import sphinx
except ImportError:
    sphinx = None

if sphinx:
    from sphinx.setup_command import BuildDoc as _BuildDoc

    class BuildDoc(_BuildDoc):
        def run(self):
            # make sure the python path is pointing to the newly built
            # code so that the documentation is built on this and not a
            # previously installed version
            build = self.get_finalized_command('build')
            sys.path.insert(0, os.path.abspath(build.build_lib))
            ff = sum((fns for dir, fns in data_files if 'ForceFields' in dir),
                     [])
            for fn in ff:
                self.copy_file(fn,
                               os.path.join(build.build_lib, fn),
                               preserve_mode=False)
            try:
                sphinx.setup_command.BuildDoc.run(self)
            except UnicodeDecodeError:
                print >>sys.stderr, "ERROR: unable to build documentation because Sphinx do not handle source path with non-ASCII characters. Please try to move the source package to another location (path with *only* ASCII characters)."            
            sys.path.pop(0)

    cmdclass['build_sphinx'] = BuildDoc

#################################################################
# Check various compiler/library properties

libraries = []
if sysconfig['LIBM'] != '':
    libraries.append('m')

macros = []
try:
    from Scientific.MPI import world
except ImportError:
    world = None
if world is not None:
    if type(world) == types.InstanceType:
        world = None
if world is not None:
    macros.append(('WITH_MPI', None))

if hasattr(ctypes.CDLL(ctypes.util.find_library('m')), 'erfc'):
    macros.append(('LIBM_HAS_ERFC', None))

if sys.platform != 'win32':
    if ctypes.sizeof(ctypes.c_long) == 8:
        macros.append(('_LONG64_', None))

if sys.version_info[0] == 2 and sys.version_info[1] >= 2:
    macros.append(('EXTENDED_TYPES', None))

#################################################################
# System-specific optimization options

low_opt = []
if sys.platform != 'win32' and 'gcc' in sysconfig['CC']:
    low_opt = ['-O0']
low_opt.append('-g')

high_opt = []
if sys.platform[:5] == 'linux' and 'gcc' in sysconfig['CC']:
    high_opt = ['-O3', '-ffast-math', '-fomit-frame-pointer',
                '-fkeep-inline-functions']
if sys.platform == 'darwin' and 'gcc' in sysconfig['CC']:
    high_opt = ['-O3', '-ffast-math', '-fomit-frame-pointer',
                '-fkeep-inline-functions', '-falign-loops=16']
if sys.platform == 'aix4':
    high_opt = ['-O4']
if sys.platform == 'odf1V4':
    high_opt = ['-O2', '-fp_reorder', '-ansi_alias', '-ansi_args']

high_opt.append('-g')

#################################################################

setup (name = package_name,
       version = pkginfo.__version__,
       description = "Molecular docking with an adaptively alchemical interaction grid",
       long_description=
"""
AlGDock is an Open Source program for molecular docking. In addition to 
low-energy poses, AlGDock also provides an estimate of the binding potential 
of mean force, or free energy of binding between a (flexible) ligand and 
rigid receptor.
""",
       author = "David Minh",
       author_email = "dminh@iit.edu",
       url = "TBA",
       license = "CeCILL-C",

       packages = ['AlGDock', 'AlGDock.ForceFields', 
                   'AlGDock.ForceFields.Cylinder',
                   'AlGDock.ForceFields.Sphere',
                   'AlGDock.ForceFields.Grid',
                   'AlGDock.Integrators',
                   'AlGDock.Integrators.VelocityVerlet',
                   'AlGDock.Integrators.HamiltonianMonteCarlo',
                   'AlGDock.Integrators.NUTS',
                   'AlGDock.Integrators.SmartDarting',
                   # 'AlGDock.Integrators.ExternalMC',
                   ],
       ext_package = 'AlGDock.'+sys.platform,

       ext_modules = [Extension('MMTK_cylinder',
                                ['AlGDock/ForceFields/Cylinder/MMTK_cylinder.pyx'],
                                extra_compile_args = compile_args,
                                include_dirs=include_dirs),
                      Extension('MMTK_sphere',
                                ['AlGDock/ForceFields/Sphere/MMTK_sphere.pyx'],
                                extra_compile_args = compile_args,
                                include_dirs=include_dirs),
                      Extension('MMTK_trilinear_grid',
                                ['AlGDock/ForceFields/Grid/MMTK_trilinear_grid.pyx'],
                                extra_compile_args = compile_args,
                                include_dirs=include_dirs),
                      Extension('MMTK_trilinear_thresh_grid',
                                ['AlGDock/ForceFields/Grid/MMTK_trilinear_thresh_grid.pyx'],
                                extra_compile_args = compile_args,
                                include_dirs=include_dirs),
                      Extension('MMTK_trilinear_isqrt_grid',
                                ['AlGDock/ForceFields/Grid/MMTK_trilinear_isqrt_grid.pyx'],
                                extra_compile_args = compile_args,
                                include_dirs=include_dirs),
                      Extension('MMTK_trilinear_transform_grid',
                                ['AlGDock/ForceFields/Grid/MMTK_trilinear_transform_grid.pyx'],
                                extra_compile_args = compile_args,
                                include_dirs=include_dirs),
                      Extension('MMTK_BSpline_grid',
                                ['AlGDock/ForceFields/Grid/MMTK_BSpline_grid.pyx'],
                                extra_compile_args = compile_args,
                                include_dirs=include_dirs),
                      Extension('MMTK_BSpline_transform_grid',
                                ['AlGDock/ForceFields/Grid/MMTK_BSpline_transform_grid.pyx'],
                                extra_compile_args = compile_args,
                                include_dirs=include_dirs),
                      Extension('MMTK_CatmullRom_grid',
                                ['AlGDock/ForceFields/Grid/MMTK_CatmullRom_grid.pyx'],
                                extra_compile_args = compile_args,
                                include_dirs=include_dirs),
                      Extension('MMTK_CatmullRom_transform_grid',
                                ['AlGDock/ForceFields/Grid/MMTK_CatmullRom_transform_grid.pyx'],
                                extra_compile_args = compile_args,
                                include_dirs=include_dirs),
                      Extension('NUTS',
                                ['AlGDock/Integrators/NUTS/NUTS.pyx'],
                                extra_compile_args = compile_args,
                                include_dirs=include_dirs),
                      ],

       data_files = data_files,
       scripts = [],

       cmdclass = cmdclass,

       command_options = {
           'build_sphinx': {
               'source_dir' : ('setup.py', 'Doc')}
           },   
       )
