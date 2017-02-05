# Use this script to compile the C module by running
#
#         python setup.py build_ext --inplace
#

from distutils.core import setup, Extension
import os, sys

compile_args = ['-fdiagnostics-color=always']
include_dirs = [] # ['../../../../Include']

from Scientific import N
try:
    num_package = N.package
except AttributeError:
    num_package = "Numeric"
if num_package == "NumPy":
    compile_args.append("-DNUMPY=1")
    import numpy.distutils.misc_util
    include_dirs.extend(numpy.distutils.misc_util.get_numpy_include_dirs())

setup (name = "TrilinearGrid",
       version = "1.0",
       description = "Trilinear Grid term for MMTK",

       py_modules = ['TrilinearGrid'],
       ext_modules = [Extension('MMTK_trilinear_grid',
                                ['MMTK_trilinear_grid.c'],
                                extra_compile_args = compile_args,
                                include_dirs=include_dirs)]
       )
