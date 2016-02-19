import sys, os, subprocess
import time
import copy
import numpy as np
from Scientific import N
import random
from MMTK import Vector

sys.path.insert(1, '/home/lspirido/Installers/AlGDock-0.0.1.01/Test_1s3vA/Integrators/ICHMC/')
sys.path.insert(1, '../../')

from MMTK import Database
from MMTK import Molecule
from MMTK.ParticleProperties import Configuration
from MMTK.ForceFields import Amber12SBForceField
from MMTK import Universe, InfiniteUniverse
from MMTK import Units
from bfuncs import *

from NUTS import NUTSIntegrator

import MMTK_DCD
from MMTK.Trajectory import Trajectory

import GCHMC

mol_name = "2butanol"

R = 8.3144621*Units.J/Units.mol/Units.K

#mol_dir = '../../prmtopcrd/'
mol_dir = '/home/lspirido/Installers/AlGDock-0.0.1.01/Test_1s3vA/prmtopcrd_2butanol/'
#mol_dir = '/home/lspirido/Installers/AlGDock-0.0.1.01/Test_1s3vA/prmtopcrd_1s3v/'
#mol_dir = '/home/lspirido/github/AlGDock/Example/prmtopcrd/'
#mol_dir = '/home/lspirido/AstexDiv_xtal/1-build/1of6/'
#mol_dir = '/home/lspirido/AstexDiv_xtal/1-build/1hnn/'
#Database.molecule_types.directory = mol_dir
newmol = Molecule(mol_name)

#mol_frcmod = os.path.join(mol_dir, 'ligand.frcmod')
mol_frcmod = mol_dir + 'ligand.frcmod'
universe = InfiniteUniverse(Amber12SBForceField(parameter_file='../../../Data/gaff.dat', mod_files=['frcmod.ff12SB', mol_frcmod]))

universe.addObject(newmol)

universe.configuration();
configuration =  universe.configuration().array

### Initialize ###
GCintegrator = GCHMC.GCHMCIntegrator(universe, mol_dir, '../../../Data/')

(confs, Es_MM, acc, ntrials, dt) = GCintegrator.Call(30, 10, 300, 0.0015, random.randint(1,300))
print "GC: ", Es_MM

GCintegrator.Clear()






