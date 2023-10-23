#!/usr/bin/env python

# Run
#
# python setup.py build
#
# to build this extension in place.
# To use it, the extension should be copied to the appropriate AlGDock
# subdirectory for extensions.

package_name = "MMTK"

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

if os.environ.get('MMTKHOME'):
  mmtk_home = os.environ['MMTKHOME']
else:
  mmtk_home = '/home/lspirido/Installers/0Work/0David/PoseFF/MMTK-2.7.9/'
  mmtk_home = '/Users/dminh/Installers/MMTK-2.7.9/'
poseff_src = os.getcwd() + '/'
#poseff_src = '/home/lspirido/tmp/PoseFF/'
#sys.path.insert(1, mmtk_home)

class Dummy:
    pass
pkginfo = Dummy()
execfile(mmtk_home + 'MMTK/__pkginfo__.py', pkginfo.__dict__)

# Check for Cython and use it if the environment variable
# MMTK_USE_CYTHON is set to a non-zero value.
use_cython = int(os.environ.get('MMTK_USE_CYTHON', '0')) != 0
if use_cython:
    try:
        from Cython.Distutils import build_ext
        use_cython = True
    except ImportError:
        use_cython = False
if not use_cython:
    from distutils.command.build_ext import build_ext
src_ext = 'pyx' if use_cython else 'c'

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
    print "MMTK needs ScientificPython 2.6 or higher"
    raise SystemExit

compile_args = []
include_dirs = [mmtk_home + 'Include']

if (int(scientific_version[1]) >= 8 or \
    (int(scientific_version[1]) == 7 and int(scientific_version[2]) >= 8)):
    netcdf_h = os.path.join(sys.prefix, 'include',
                            'python%d.%d' % sys.version_info[:2],
                            'Scientific', 'netcdf.h')
    print "netcdf.h path", netcdf_h
    if os.path.exists(netcdf_h):
        compile_args.append("-DUSE_NETCDF_H_FROM_SCIENTIFIC=1")
        include_dirs.append(os.path.join(sys.prefix, 'include', 'python%d.%d' % sys.version_info[:2])) # EU
else:
    # Take care of the common problem that netcdf is in /usr/local but
    # /usr/local/include is not on $CPATH.
    if os.path.exists('/usr/local/include/netcdf.h'):
        include_dirs.append('/usr/local/include')

from Scientific import N
try:
    num_package = N.package
except AttributeError:
    num_package = "Numeric"
if num_package == "NumPy":
    compile_args.append("-DNUMPY=1")
    import numpy.distutils.misc_util
    include_dirs.extend(numpy.distutils.misc_util.get_numpy_include_dirs())

headers = glob(os.path.join ("Include", "MMTK", "*.h"))

paths = [os.path.join(mmtk_home + 'MMTK', 'ForceFields'),
         os.path.join(mmtk_home + 'MMTK', 'ForceFields', 'Amber'),
         os.path.join(mmtk_home + 'MMTK', 'Database', 'Atoms'),
         os.path.join(mmtk_home + 'MMTK', 'Database', 'Groups'),
         os.path.join(mmtk_home + 'MMTK', 'Database', 'Molecules'),
         os.path.join(mmtk_home + 'MMTK', 'Database', 'Complexes'),
         os.path.join(mmtk_home + 'MMTK', 'Database', 'Proteins'),
         os.path.join(mmtk_home + 'MMTK', 'Database', 'PDB'),
         os.path.join(mmtk_home + 'MMTK', 'Tools', 'TrajectoryViewer')]
data_files = []
for dir in paths:
    files = []
    for f in glob(os.path.join(dir, '*')):
        if f[-3:] != '.py' and f[-4:-1] != '.py' and os.path.isfile(f):
            files.append(f)
    data_files.append((dir, files))


class ModifiedFileList(FileList):

    #def findall(self, dir=os.curdir):
    def findall(self, dir=mmtk_home):
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
        subprocess.call([sys.executable, 'Tests/all_tests.py'],
                        env={'PYTHONPATH': self.build_lib,
                             'MMTKDATABASE': 'MMTK/Database'})

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
       description = "Molecular Modelling Toolkit",
       long_description=
"""
The Molecular Modelling Toolkit (MMTK) is an Open Source program
library for molecular simulation applications. It provides the most
common methods in molecular simulations (molecular dynamics, energy
minimization, normal mode analysis) and several force fields used for
biomolecules (Amber 94, Amber 99, several elastic network
models). MMTK also serves as a code basis that can be easily extended
and modified to deal with non-standard situations in molecular
simulations.
""",
       author = "Konrad Hinsen",
       author_email = "hinsen@cnrs-orleans.fr",
       url = "http://dirac.cnrs-orleans.fr/MMTK/",
       license = "CeCILL-C",

       package_dir = {'' : mmtk_home},
       #packages = ['MMTK', 'MMTK.ForceFields', 'MMTK.ForceFields.Amber',
       #            'MMTK.NormalModes', 'MMTK.Tk', 'MMTK.Tools',
       #            'MMTK.Tools.TrajectoryViewer'],
       headers = headers,
       ext_package = 'MMTK.'+sys.platform,
       ext_modules = [Extension('MMTK_pose',
                                [poseff_src + 'MMTK_pose.c', poseff_src + 'pose.c',
                                 mmtk_home + 'Src/bonded.c', mmtk_home + 'Src/nonbonded.c',
                                 mmtk_home + 'Src/ewald.c', mmtk_home + 'Src/sparsefc.c'],
                                extra_compile_args = compile_args + high_opt,
                                include_dirs=include_dirs + ['Src'],
                                define_macros = [('SERIAL', None),
                                                 ('VIRIAL', None),
                                                 ('MACROSCOPIC', None)]
                                                + macros,
                                libraries=libraries),
                      ],

       data_files = data_files,
       #scripts = ['tviewer'],

       cmdclass = cmdclass,

       #command_options = {
       #    'build_sphinx': {
       #        'source_dir' : ('setup.py', 'Doc')}
       #    },   
       )
