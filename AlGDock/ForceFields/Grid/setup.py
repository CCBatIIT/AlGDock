from distutils.core import setup
from Cython.Build import cythonize

setup(
  name = 'Tricubic grid',
  ext_modules = cythonize("MMTK_Tricubic_grid.pyx"),
)
