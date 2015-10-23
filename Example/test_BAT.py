
import os
import sys
sys.path.append('/Users/dminh/canopy/lib/python2.7/site-packages/AlGDock/darwin')

import numpy as np
from Scientific.Geometry.Objects3D import Sphere, Cone, Plane, Line, \
                                          rotatePoint
from Scientific.Geometry import Vector
import MMTK

MMTK.Database.molecule_types.directory = os.path.dirname(os.path.abspath('ligand.db'))
universe = MMTK.Universe.InfiniteUniverse()
molecule = MMTK.Molecule(os.path.basename('ligand.db'))
universe.addObject(molecule)

original_xyz = np.copy(universe.configuration().array)
import BAT

self = BAT.converter(universe, molecule)

# This tests a conversion to BAT coordinates and back
BAT = self.BAT(original_xyz, extended=True)
new_xyz = self.Cartesian(BAT)
print sum(sum(new_xyz - original_xyz))
