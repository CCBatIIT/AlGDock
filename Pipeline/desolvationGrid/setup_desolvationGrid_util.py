# python setup_desolvationGrid_util.py  build_ext --inplace
from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

extensions = [
  Extension("desolvationGrid_util", ["desolvationGrid_util.pyx"],
    include_dirs = [numpy.get_include()])]

setup(
    ext_modules = cythonize(extensions)
)
