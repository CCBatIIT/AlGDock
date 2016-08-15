import sys
import os
import subprocess as sproc
import glob
import itertools
import argparse
from os.path import expanduser

# Options - Not necessary for now
parser = argparse.ArgumentParser()
parser.add_argument('--include_packages_dirs', nargs='+', default=None, \
  help='Adds a list of directories to be searched for packages')
args_in = parser.parse_args()

home = expanduser("~")
packages_dirs = [os.path.join(home, "Installers"), os.path.join(home, "Downloads"), "/usr/local", "usr"]
if args_in.include_packages_dirs:
  packages_dirs = [os.path.abspath(i) for i in args_in.include_packages_dirs] + packages_dirs

# User defined paths to required packages
execfile('paths.py')

# Compiler flags
incl_flags = {}
rdirs_flags = {}
ldirs_flags = {}
libs_flags = {}

incl_flags['CURR']=["-I./"]
incl_flags['GMOL']=["-Igmolmodel"]
incl_flags['SIMTK']=["-I" + os.path.join(SIMBODYHOME, "include"), 
  "-I" + os.path.join(SIMBODYHOME, "include/molmodel/internal"), 
  "-I" + SIMBODYSRC, "-I" + os.path.join(MOLMODELHOME, "src"),  "-I" + os.path.join(MOLMODELHOME, "include")]
incl_flags['MMTK']=["-I" + os.path.join(PYTHONHOME, "include/python2.7"), "-I" + PYTHON2_6_INCLUDE]

incl_flags['NUMPY']=["-I" + os.path.join(PYTHONHOME, "lib/python2.7/site-packages/numpy/core/include")]
incl_flags['EIGEN']=["-I" + os.path.join(EIGENHOME)]
incl_flags['BOOST']=["-I" + os.path.join(BOOSTHOME)]

rdirs_flags['CURR']=["-Wl,-rpath,./"]
rdirs_flags['GEN']=["-Wl,-rpath,/usr/lib64"]
rdirs_flags['GMOL']=["-Wl,-rpath,gmolmodel"]
rdirs_flags['SIMTK']=["-Wl,-rpath," + os.path.join(SIMBODYHOME, "lib64")]
rdirs_flags['MMTK']=["-Wl,-rpath," + os.path.join(MMTKHOME, "MMTK/linux2")]
rdirs_flags['NUMPY']=["-Wl,-rpath," + os.path.join(PYTHONHOME, "lib/python2.7/site-packages/numpy")]

ldirs_flags['CURR']=["-L./"]
ldirs_flags['GEN']=["-L" + "/usr/lib64"]
ldirs_flags['SIMTK']=["-L" + os.path.join(SIMBODYHOME, "lib64")]
ldirs_flags['GMOL']=["-Lgmolmodel/"]
ldirs_flags['MMTK']=["-L" + os.path.join(MMTKHOME, "MMTK/linux2")]
ldirs_flags['NUMPY']=["-L" + os.path.join(PYTHONHOME, "lib/python2.7/site-packages/numpy")]

libs_flags['GMOL']=["-lgmolmodel"]
libs_flags['SIMTK']=["-lSimTKsimbody", "-lSimTKmath", "-lSimTKcommon"]
libs_flags['MOLMO']=["-lSimTKmolmodel"]
libs_flags['NUMPY']=["-l_compiled_base"]
#libs_flags['BOOST']=["-lboost_python"]
libs_flags['DEPS']=['-llapack', '-lblas', '-lpthread', '-ldl', '-lrt']

cxx_comm = ['g++']
obj_flags = ['-Wall', '-fPIC', '-g', '-c']

# Compile objects for Gmolmodel library
os.chdir('gmolmodel')
a = 0
srcFNs = [s for s in glob.glob("*.cpp") if (s != 'simmain.cpp' and s != 'MidVVIntegrator.cpp')]
FNroots = [os.path.splitext(os.path.basename(s))[0] for s in srcFNs]

if os.path.isfile('libgmolmodel.so'):
  sproc.call(['rm', 'libgmolmodel.so', 'libgmolmodel.so.1', 'libgmolmodel.so.1.0'])
if os.path.isfile('compil.log'):
  sproc.call(['rm', 'compil.log'])
sproc.call(['touch', 'compil.log'])

for FNroot in FNroots:
  srcFN = FNroot + '.cpp'
  hdrFN = FNroot + '.hpp'
  objFN = FNroot + '.o'
  if FNroot != 'bSystem':
    command = cxx_comm + obj_flags + [srcFN] + incl_flags['SIMTK'] + incl_flags['BOOST'] + ['-o', objFN]
    if os.path.isfile(objFN): # Check existence
      if (os.path.getmtime(srcFN) > os.path.getmtime(objFN)) or \
        (os.path.getmtime(hdrFN) > os.path.getmtime(objFN)): # Check modification time
        with open('compil.log', 'a') as f:
          sproc.call(['echo'] + command)
          f.write(' '.join(command))
          a += sproc.call(command, stdout=f, stderr=f)
    else: # Object does not exist
      with open('compil.log', 'a') as f:
        sproc.call(['echo'] + command)
        f.write(' '.join(command))
        a += sproc.call(command, stdout=f, stderr=f)
  else:
    command = cxx_comm + ['-std=gnu++0x'] + obj_flags + [srcFN, '-DNUMPY=1'] + \
      incl_flags['SIMTK'] + incl_flags['BOOST'] + incl_flags['MMTK'] + incl_flags['EIGEN'] + incl_flags['NUMPY'] + \
      ['-l_compiled_base', '-llapack', '-o', objFN]
    if os.path.isfile(objFN): # Check existence
      if (os.path.getmtime(srcFN) > os.path.getmtime(objFN)) or \
        (os.path.getmtime(hdrFN) > os.path.getmtime(objFN)) or \
        (os.path.getmtime('MidVVIntegrator.cpp') > os.path.getmtime(objFN)) or \
        (os.path.getmtime('MidVVIntegrator.hpp') > os.path.getmtime(objFN)): # Check modification time
        with open('compil.log', 'a') as f:
          sproc.call(['echo'] + command)
          f.write(' '.join(command))
          a += sproc.call(command, stdout=f, stderr=f)
    else: # Object does not exist
        with open('compil.log', 'a') as f:
          sproc.call(['echo'] + command)
          f.write(' '.join(command))
          a += sproc.call(command, stdout=f, stderr=f)  
os.chdir('../')  

if a != 0:
  print "\n===", a, "objects failed to compile.  ===\n"
  sys.exit(1)
else:
  print "\n=== Gmolmodel objects compiled successfully. ===\n"

# Compile Gmolmodel library
if os.path.isfile('compil.log'):
  sproc.call(['rm', 'compil.log'])
sproc.call(['touch', 'compil.log'])
objs = ['gmolmodel/' + s + '.o' for s in FNroots]

needs_compiling = False
if os.path.isfile('libgmolmodel.so.1.0'):
  objs_times = [os.path.getmtime(o) for o in objs]
  lib_time   = os.path.getmtime('libgmolmodel.so.1.0')
  if [t for t in objs_times if t > lib_time]: 
    needs_compiling = True
else:
  needs_compiling = True

if needs_compiling:
  add_flags = ['-std=gnu++0x', '-g', '-fPIC', '-shared', '-Wl,-soname,libgmolmodel.so.1.0', '-o', 'libgmolmodel.so.1.0']
  command = cxx_comm + add_flags + objs + incl_flags['SIMTK'] + ['-Wl,-Bdynamic'] + \
    ldirs_flags['SIMTK'] + ldirs_flags['GEN'] + \
    rdirs_flags['SIMTK'] + rdirs_flags['GEN'] + \
     libs_flags['SIMTK']  + libs_flags['MOLMO'] + libs_flags['DEPS']

  with open('gmolmodel/compil.log', 'a') as f:
    sproc.call(['echo'] + command)
    f.write(' '.join(command))
    a = sproc.call(command, stdout=f, stderr=f)
    if a != 0:
      print "\n=== libgmolmodel.so.1.0 compilation failed. Exit code:", a, "===\n"
      sys.exit(1)
    else:
      print "\n=== Gmolmodel library compiled successfully. ===\n"
      sproc.call(['ln', '-sf', 'libgmolmodel.so.1.0', 'libgmolmodel.so.1'])
      sproc.call(['ln', '-sf', 'libgmolmodel.so.1', 'libgmolmodel.so'])


# Compile GCHMC.so (C++ Python extension)
# Boost uses bjam program to build Python extension
# Bjam requires boost-build.jam & Jamroot files
# Write boost-build.jam
bb_script="""
boost-build {0}/tools/build/v2 ;
""".format(BOOSTHOME)
bb_FN = 'boost-build.jam'
bb_F = open(bb_FN,'w')
bb_F.write(bb_script)
bb_F.close()

# Write Jamroot
jam_script = """
# Copyright David Abrahams 2006. Distributed under the Boost
# Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
import python ;
if ! [ python.configured ]
{{
    ECHO "notice: no Python configured in user-config.jam" ;
    ECHO "notice: will use default configuration" ;
    using python ;
}}

# Specify the path to the Boost project.  If you move this project,
# adjust this path to refer to the Boost root directory.
use-project boost
  : {0} ;

# Set up the project-wide requirements that everything uses the
# boost_python library from the project whose global ID is
# /boost/python.
project
  : requirements <library>/boost/python//boost_python
                 <implicit-dependency>/boost//headers

  : usage-requirements <implicit-dependency>/boost//headers
  ;
lib gmolmodel
    :
    : <file>./libgmolmodel.so <variant>release
    ;

# Declare the extension modules.  You can specify multiple
# source files after the colon separated by spaces.
python-extension GCHMC : ./gmolmodel/simmain.cpp ./libgmolmodel.so ;

# Put the extension and Boost.Python DLL in the current directory, so
# that running script by hand works.
install convenient_copy
  : GCHMC
  : <install-dependencies>on <install-type>SHARED_LIB <install-type>PYTHON_EXTENSION
    <location>.
  ;

# A little "rule" (function) to clean up the syntax of declaring tests
# of these extension modules.
local rule run-test ( test-name : sources + )
{{
    import testing ;
    testing.make-test run-pyd : $(sources) : : $(test-name) ;
}}

# Declare test targets
run-test GCHMCIntegrator : GCHMC test_GCHMC.py ;
""".format(BOOSTHOME)

jam_FN = 'Jamroot'
jam_F = open(jam_FN,'w')
jam_F.write(jam_script)
jam_F.close()

# Compile GCHMC.so using bjam
debug_flags = ['-Wall', "-g"]
optim_flags = ["-O2"]
add_flags = ['-Xlinker', '-export-dynamic', '-std=gnu++0x']
mmtk_add_flags = []
mmtk_add_flags.append("-pthread")
mmtk_add_flags.append("-fno-strict-aliasing")
mmtk_add_flags.append("-DNDEBUG")
mmtk_add_flags.append("-DLIBM_HAS_ERFC")
mmtk_add_flags.append("-D_LONG64_")
mmtk_add_flags.append("-DEXTENDED_TYPES")
mmtk_add_flags.append("-DUSE_NETCDF_H_FROM_SCIENTIFIC=1")
mmtk_add_flags.append("-DNUMPY=1")
mmtk_add_flags.append("-ffast-math")
mmtk_add_flags.append("-fomit-frame-pointer")
#mmtk_add_flags.append("-keep-inline-functions")


bjam = [os.path.join(BOOSTHOME, 'bjam')] 
bjam_arg1_list =  debug_flags + add_flags + mmtk_add_flags + \
  [l for l in list(itertools.chain.from_iterable(incl_flags.values()))] + \
  [l for l in list(itertools.chain.from_iterable(rdirs_flags.values()))] + \
  [l for l in list(itertools.chain.from_iterable(ldirs_flags.values()))] + \
  [l for l in list(itertools.chain.from_iterable(libs_flags.values()))]
bjam_arg1_list[0] = "cxxflags=" + bjam_arg1_list[0]
bjam_arg1_list[-1] += ""
bjam_arg1 = [" ".join(bjam_arg1_list)]

command = bjam + bjam_arg1
with open('compil.log', 'a') as f:
  sproc.call(['echo'] + command)
  f.write(' '.join(command))
  a = sproc.call(command, stdout=f, stderr=f)

if a != 0:
  print a, "\n=== GCHMC extension failed to compile or to be tested. ===\n"
  sys.exit(1)
else:
  print "\n=== GCHMC extension compiled and tested successfully. ===\n"







