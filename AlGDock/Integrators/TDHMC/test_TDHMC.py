import sys, os, subprocess
import copy
import numpy as np
import random

mol_FN = "ligand.db"
mol_dir = os.path.abspath('2but')

import MMTK

MMTK.Database.molecule_types.directory = mol_dir
from MMTK import Database
from MMTK import Molecule
from MMTK.ParticleProperties import Configuration
from MMTK.ForceFields import Amber12SBForceField
from MMTK import Universe, InfiniteUniverse
from MMTK import Units

R = 8.3144621*Units.J/Units.mol/Units.K

import TDHMC

parm_dir = '/share/apps/amber/16/dat/leap/parm/' # Amber parameters directory
newmol = Molecule(mol_FN)
mol_frcmod = os.path.join(mol_dir, 'ligand.frcmod')
gaff_FN = os.path.join(mol_dir,'gaff.dat')
universe = InfiniteUniverse(Amber12SBForceField(parameter_file=gaff_FN, mod_files=['frcmod.ff12SB', mol_frcmod]))
universe.addObject(newmol)

universe.configuration();
configuration =  universe.configuration().array

GCintegrator = TDHMC.TDHMCIntegrator(universe, mol_dir, parm_dir)
(confs, Es_MM, acc, ntrials, dt) = GCintegrator.Call(30, 10, 300, 0.0015, random.randint(1,300), 0, 1, 0.5)
print "GC: ", Es_MM
GCintegrator.Clear()






