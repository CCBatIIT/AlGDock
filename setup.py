#!/usr/bin/env python

package_name = "AlGDock"

from setuptools import setup, Command, Extension
from setuptools.command.build import build
from setuptools.command.sdist import sdist
try:
    from setuptools.command.install_data import install_data
except ModuleNotFoundError:
    from setuptools._distutils.command.install_data import install_data

from distutils import dir_util
import sysconfig as python_sysconfig

import os
import sys
import ctypes
import ctypes.util
from glob import glob

config_vars = python_sysconfig.get_config_vars()

class Dummy:
    pass
pkginfo = Dummy()
#execfile('AlGDock/__pkginfo__.py', pkginfo.__dict__)
#execfile('AlGDock/path_tools.py')

try:
    with open("AlGDock/__pkginfo__.py") as f:
        exec(f.read(), pkginfo.__dict__)
    with open("AlGDock/path_tools.py") as f:
        exec(f.read())
except FileNotFoundError as e:
    print(f"Error: Required file not found: {e}. Make sure AlGDock/__pkginfo__.py and AlGDock/path_tools.py exist.")
    sys.exit(1)

from site import USER_BASE as userbase

# Check module dependencies

# Check for Cython
try:
  from Cython.Distutils import build_ext
  cython_ok = True
except ImportError:
  cython_ok = False
if not cython_ok:
  print('AlGDock requires Cython')
  sys.exit(1) # Changed to sys.exit(1) for a non-zero exit code

compile_args = ['-fdiagnostics-color=always']
include_dirs = ['Include'] # Base include directory

import numpy
include_dirs_config = ['Include']
numpy_include_path = numpy.get_include()
include_dirs_config.append(numpy_include_path)

if os.path.exists('/usr/local/include/netcdf.h'):   #Generic NetCDF path checks
    include_dirs.append('/usr/local/include')

if ('NETCDF_PREFIX' in os.environ):
  include_dirs.append(os.path.join(os.environ['NETCDF_PREFIX'],'include'))

py_version_short = f"python{sys.version_info.major}.{sys.version_info.minor}"
user_include_path_python = os.path.join(userbase, 'include', py_version_short)
if os.path.exists(user_include_path_python):
    include_dirs.append(user_include_path_python)

# Remove duplicates from include_dirs
include_dirs = sorted(list(set(include_dirs)))

headers = []
paths = [os.path.join('AlGDock', 'ForceFields', 'Cylinder'),
         os.path.join('AlGDock', 'ForceFields', 'Sphere'),
         os.path.join('AlGDock', 'ForceFields', 'Grid'),
         os.path.join('AlGDock', 'ForceFields', 'OBC'),
         os.path.join('AlGDock', 'ForceFields', 'OpenMM'),
         os.path.join('AlGDock', 'ForceFields', 'Pose'),
         os.path.join('AlGDock', 'ForceFields', 'ElectricField'),
         os.path.join('AlGDock', 'Integrators', 'ExternalMC'),
         os.path.join('AlGDock', 'Integrators', 'HamiltonianMonteCarlo'),
         os.path.join('AlGDock', 'Integrators', 'NUTS'),
         os.path.join('AlGDock', 'Integrators', 'SmartDarting'),
         os.path.join('AlGDock', 'Integrators', 'VelocityVerlet')]

data_files = []
for p_dir in paths:
    files_in_dir = []
    if os.path.isdir(p_dir):
        for f_path in glob(os.path.join(p_dir, '*')):
            if os.path.isfile(f_path) and not f_path.endswith(('.py', '.pyc', '.so', '.dylib', '.pyd')):
                files_in_dir.append(f_path)
    if files_in_dir:
        data_files.append((p_dir, files_in_dir))

from setuptools._distutils.filelist import FileList as DistutilsFileList # Use the one from distutils
class ModifiedFileList(DistutilsFileList):  # Changed inheritance
    def findall(self, dir=os.curdir):
        from stat import ST_MODE, S_ISREG, S_ISDIR, S_ISLNK
        file_list = []
        stack = [dir]
        pop = stack.pop
        push = stack.append
        while stack:
            current_dir = pop()
            try:
                names = os.listdir(current_dir)
            except OSError:
                continue
            for name in names:
                if current_dir != os.curdir:
                    fullname = os.path.join(current_dir, name)
                else:
                    fullname = name

                try:
                    stat_info = os.stat(fullname)
                except OSError:
                    continue

                mode = stat_info[ST_MODE]
                if S_ISREG(mode):
                    file_list.append(fullname)
                elif S_ISDIR(mode) and not S_ISLNK(mode):
                    file_list.append(fullname)
                    push(fullname)
        self.allfiles = file_list

class modified_sdist(sdist):
    def run(self):
        self.filelist = ModifiedFileList()
        self.check_metadata()
        self.get_file_list()
        if getattr(self, 'manifest_only', False):
            return
        self.make_distribution()

    def make_release_tree(self, base_dir, files):
        self.mkpath(base_dir)
        dir_util.create_tree(base_dir, files, dry_run=self.dry_run)  # Use distutils dir_util

        if not files:
            self.warn("no files to distribute -- empty manifest?")
            return

        self.announce(f"making hard links in {base_dir}..."
                      if hasattr(os, 'link') else f"copying files to {base_dir}...")

        for file_path in files:
            dest = os.path.join(base_dir, file_path)
            if os.path.isfile(file_path):
                parent_dest_dir = os.path.dirname(dest)
                if not os.path.exists(parent_dest_dir):
                    self.mkpath(parent_dest_dir)
                self.copy_file(file_path, dest, link='hard' if hasattr(os, 'link') else None)
            elif os.path.isdir(file_path):
                # If ModifiedFileList includes dirs, ensure they are created if not already by copy_file's dest
                # dir_util.create_tree should generally handle this, but explicit mkpath can be a fallback.
                self.mkpath(dest)
            else:
                self.warn(f"'{file_path}' not a regular file or directory -- skipping")


class modified_build(build):
    def has_sphinx(self):
        setup_dir = os.path.dirname(os.path.abspath(__file__))
        return os.path.isdir(os.path.join(setup_dir, 'Doc'))
    sub_commands = build.sub_commands + [('build_sphinx', has_sphinx)]


class modified_install_data(install_data):
    def run(self):
        install_cmd = self.get_finalized_command('install')
        self.install_dir = getattr(install_cmd, 'install_lib')
        return super().run()  # Use super() for Python 3


class test(Command):
    user_options = []

    def initialize_options(self):
        self.build_lib = None

    def finalize_options(self):
        self.set_undefined_options('build', ('build_lib', 'build_lib'))

    def run(self):
        import subprocess  # sys is already imported
        self.run_command('build_py')
        self.run_command('build_ext')

        # Flatten the list of force field data files
        ff_files_to_copy = []
        for _, fns in data_files:  # Iterate through (dir, files_in_dir)
            if 'ForceFields' in _:  # Check if the directory string contains 'ForceFields'
                ff_files_to_copy.extend(fns)

        for fn in ff_files_to_copy:
            dest_file = os.path.join(self.build_lib, fn)
            dest_dir = os.path.dirname(dest_file)
            if not os.path.exists(dest_dir):
                os.makedirs(dest_dir)
            self.copy_file(fn, dest_file, preserve_mode=False)

cmdclass = {
    'build': modified_build,
    'sdist': modified_sdist,
    'install_data': modified_install_data,
    'build_ext': build_ext,
    'test': test
}

# Build the sphinx documentation if Sphinx is available
try:
    import sphinx
    from sphinx.setup_command import BuildDoc as _BuildDoc


    class BuildDoc(_BuildDoc):
        def run(self):
            build = self.get_finalized_command('build')
            sys.path.insert(0, os.path.abspath(build.build_lib))

            ff_files_to_copy = []
            for _, fns in data_files:
                if 'ForceFields' in _:
                    ff_files_to_copy.extend(fns)

            for fn in ff_files_to_copy:
                dest_file = os.path.join(build.build_lib, fn)
                dest_dir = os.path.dirname(dest_file)
                if not os.path.exists(dest_dir):
                    os.makedirs(dest_dir)
                self.copy_file(fn, dest_file, preserve_mode=False)

            try:
                super().run()
            except UnicodeDecodeError:
                print(
                    "ERROR: unable to build documentation because Sphinx may not handle source paths with non-ASCII characters. Please try moving the source package to a path with only ASCII characters.")
            sys.path.pop(0)

    cmdclass['build_sphinx'] = BuildDoc
except ImportError:
    sphinx = None

libraries = []
if config_vars.get('LIBM') and config_vars['LIBM'] != '':
    libraries.append('m')

macros = []

try:
    libm_path = ctypes.util.find_library('m')
    if libm_path:
        libm = ctypes.CDLL(libm_path)
        if hasattr(libm, 'erfc'):
            macros.append(('LIBM_HAS_ERFC', None))
except (OSError, AttributeError):
    print("Warning: Could not check for erfc in libm.")

if sys.platform != 'win32':
    if ctypes.sizeof(ctypes.c_long) == 8:  # Check for 64-bit long type
        macros.append(('_LONG64_', None))

low_opt = []
cc_name = os.path.basename(config_vars.get('CC', 'gcc'))  # Default to 'gcc' if CC not found
if sys.platform != 'win32' and 'gcc' in cc_name or 'clang' in cc_name:
    low_opt = ['-O0']
low_opt.append('-g')  # Debug symbols

high_opt = []
if sys.platform.startswith('linux') and ('gcc' in cc_name or 'clang' in cc_name):
    high_opt = ['-O3', '-ffast-math', '-fomit-frame-pointer']
if sys.platform == 'darwin' and ('gcc' in cc_name or 'clang' in cc_name):
    high_opt = ['-O3', '-ffast-math', '-fomit-frame-pointer']
if sys.platform == 'aix4':
    high_opt = ['-O4']
if sys.platform == 'odf1V4':
    high_opt = ['-O2', '-fp_reorder', '-ansi_alias', '-ansi_args']

high_opt.append('-g')

ext_module_name_and_path = [
    ('MMTK_sphere', ['AlGDock/ForceFields/Sphere/MMTK_sphere.pyx']),
    ('MMTK_cylinder', ['AlGDock/ForceFields/Cylinder/MMTK_cylinder.pyx']),
    ('MMTK_trilinear_grid', ['AlGDock/ForceFields/Grid/MMTK_trilinear_grid.c']),
    ('MMTK_trilinear_one_fourth_grid', ['AlGDock/ForceFields/Grid/MMTK_trilinear_one_fourth_grid.c']),
    ('MMTK_OBC', ['AlGDock/ForceFields/OBC/MMTK_OBC.c',
                  'AlGDock/ForceFields/OBC/ObcParameters.cpp',
                  'AlGDock/ForceFields/OBC/ObcWrapper.cpp',
                  'AlGDock/ForceFields/OBC/ReferenceForce.cpp',
                  'AlGDock/ForceFields/OBC/ReferenceObc.cpp']),
    ('MMTK_OBC_desolv', ['AlGDock/ForceFields/OBC/MMTK_OBC_desolv.c',
                         'AlGDock/ForceFields/OBC/ObcParameters.cpp',
                         'AlGDock/ForceFields/OBC/ObcWrapper.cpp',
                         'AlGDock/ForceFields/OBC/ReferenceForce.cpp',
                         'AlGDock/ForceFields/OBC/ReferenceObc.cpp']),
    ('MMTK_pose', ['AlGDock/ForceFields/Pose/MMTK_pose.c',
                   'AlGDock/ForceFields/Pose/pose.c']),
    ('MMTK_electric_field', ['AlGDock/ForceFields/ElectricField/MMTK_electric_field.c']),
    ('MMTK_electric_field_z', ['AlGDock/ForceFields/ElectricField/MMTK_electric_field_z.c']),
    ('NUTS', ['AlGDock/Integrators/NUTS/NUTS.pyx']),
    ('SmartDarting', ['AlGDock/Integrators/SmartDarting/SmartDarting.pyx']),
    ('BAT', ['Src/BAT.pyx']),
    ('repX', ['Src/repX.pyx'])
]

# Conditional compilation blocks for unused/obsolete modules are kept as is (all False).

setup(name=package_name,
      version=pkginfo.__version__,
      description="Molecular docking with an adaptive alchemical interaction grid",
      long_description=
      """
       AlGDock is an Open Source program for molecular docking. In addition to
       low-energy poses, AlGDock also provides an estimate of the binding potential
       of mean force, or free energy of binding between a (flexible) ligand and
       rigid receptor.
       """,
      author="David Minh",
      author_email="dminh@iit.edu",
      url="TBA",  # Consider updating this
      license="MIT",

      packages=['AlGDock',
                'AlGDock.ForceFields',
                'AlGDock.ForceFields.Cylinder',
                'AlGDock.ForceFields.Sphere',
                'AlGDock.ForceFields.Grid',
                'AlGDock.ForceFields.OBC',
                'AlGDock.ForceFields.OpenMM',
                'AlGDock.ForceFields.Pose',
                'AlGDock.ForceFields.ElectricField',
                'AlGDock.Integrators',
                'AlGDock.Integrators.ExternalMC',
                'AlGDock.Integrators.HamiltonianMonteCarlo',
                'AlGDock.Integrators.NUTS',
                'AlGDock.Integrators.SmartDarting',
                'AlGDock.Integrators.VelocityVerlet'],
      # ext_package defines where the .so files will be placed within the package structure
      ext_package='AlGDock',
      ext_modules=[Extension(f"AlGDock.{name}",  # Prepend AlGDock for proper packaging.
                             path,
                             extra_compile_args=compile_args + high_opt,
                             include_dirs=include_dirs,
                             define_macros=[('SERIAL', None), ('VIRIAL', None), ('MACROSCOPIC', None)] + macros,
                             libraries=libraries,
                             language='c++' if any(s.endswith('.cpp') for s in path) else 'c'
                             )
                   for (name, path) in ext_module_name_and_path],
      data_files=data_files,
      scripts=[],

      cmdclass=cmdclass,

      command_options={
          'build_sphinx': {
              'source_dir': ('setup.py', 'Doc'),
              'build_dir': ('setup.py', os.path.join('build', 'sphinx'))
          }
      },
      python_requires='>=3.8',
      install_requires=[
          'numpy',
          'cython',
      ],
      )
