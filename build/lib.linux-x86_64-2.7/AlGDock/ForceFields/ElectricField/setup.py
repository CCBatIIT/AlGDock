# Use this script to compile the C module by running
#
#         python setup.py build_ext --inplace
#
# Then run "python test.py"

from distutils.core import setup, Extension
import os, sys

compile_args = []
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

setup (name = "ElectricField",
       version = "1.0",
       description = "Electric field term for MMTK",

       py_modules = ['ElectricField'],
       ext_modules = [Extension('MMTK_electric_field',
                                ['MMTK_electric_field.c'],
                                extra_compile_args = compile_args,
                                include_dirs=include_dirs),
                      Extension('MMTK_electric_field_z',
                                ['MMTK_electric_field_z.c'],
                                extra_compile_args = compile_args,
                                include_dirs=include_dirs)                                ]
       )
