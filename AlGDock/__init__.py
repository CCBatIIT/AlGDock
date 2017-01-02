# Add shared library path to sys.path

import os, sys
pkg_loc = os.path.split(__file__)[0]
sys.path.append(os.path.join(pkg_loc, sys.platform))
# sys.path.append(os.path.join(pkg_loc,'Integrators','CDHMC'))

# Load external path search routine

execfile(os.path.join(pkg_loc,'__pkginfo__.py'))
execfile(os.path.join(pkg_loc,'_external_paths.py'))

del os
del sys
