# Use this script to compile the C module by running
#
#        cpython setup.py build_ext --inplace
#
# Then run "python test.py"

from distutils.core import setup
from Cython.Build import cythonize

import os, sys

compile_args = []
include_dirs = ['.'] # '../../../../Include']

from Scientific import N
try:
    num_package = N.package
except AttributeError:
    num_package = "Numeric"
if num_package == "NumPy":
    compile_args.append("-DNUMPY=1")
    import numpy.distutils.misc_util
    include_dirs.extend(numpy.distutils.misc_util.get_numpy_include_dirs())

setup (name = "OBC",
       version = "1.0",
       description = "Onufriev-Bashford-Case Generalized Born term for MMTK",

       py_modules = ['OBC'],
       ext_modules = cythonize(
          "MMTK_OBC.pyx",
          sources = ['MMTK_OBC.c', 'ObcParameters.cpp', 'ReferenceForce.cpp', 'ReferenceObc.cpp'],
          language="c++",
          extra_compile_args = compile_args,
          include_dirs=include_dirs))
