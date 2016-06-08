import os
import sys
from copy import deepcopy

import MMTK
MMTK.Database.molecule_types.directory = os.getcwd()
from MMTK.Trajectory import Trajectory, TrajectoryOutput, StandardLogOutput
from MMTK.Dynamics import VelocityVerletIntegrator, TranslationRemover, RotationRemover

import AlGDock.IO
from AlGDock.Integrators.VelocityVerlet.VelocityVerlet import VelocityVerletIntegrator
from PoseFF import PoseForceField
import numpy as np
import BAT
sys.path.append('/home/lspirido/tests/temp_AlGDock_distrib_test/AlGDock/AlGDock/Integrators/GCHMC/')
import random
from GCHMC import GCHMCIntegrator

amberFF = MMTK.ForceFields.Amber12SBForceField(\
  parameter_file='/export/apps/amber/14/dat/leap/parm/gaff.dat', mod_files=['1s3v.frcmod'])
#universe = MMTK.InfiniteUniverse(amberFF)
universe = MMTK.InfiniteUniverse(MMTK.ForceFields.Amber12SBForceField(\
  parameter_file='/export/apps/amber/14/dat/leap/parm/gaff.dat', mod_files=['1s3v.frcmod']))


molecule = MMTK.Molecule('1s3v.db')
universe.addObject(molecule)

universe.configuration();
configuration = universe.configuration().array

cart2bat = BAT.converter(universe, molecule)

XYZ = {}
for a in universe.atomList():
  XYZ[a.index] = np.array([a.position()[0], a.position()[1], a.position()[2]])
  #print "atoms", a.name, a.index, XYZ[a.index]

torsions = []
for t in cart2bat._torsionIndL:
  torsions.append([t[0], t[1], t[2], t[3], 0, (BAT.BAT4(XYZ[t[0]], XYZ[t[1]], XYZ[t[2]], XYZ[t[3]])[2]), 0.03490659/2, 1000])

rootIs = np.array([cart2bat.rootInd[0], cart2bat.rootInd[1], cart2bat.rootInd[2]])
bat = cart2bat.BAT(XYZ, extended=True)

poseInp = []
extPoseInp = np.concatenate(( rootIs, XYZ[rootIs[0]], XYZ[rootIs[1]], XYZ[rootIs[2]], bat[3:6]))
poseInp = [extPoseInp, torsions]
print poseInp

poseFF = PoseForceField(poseInp)
#compoundFF = amberFF + poseFF
compoundFF = amberFF
universe.setForceField(compoundFF)
#
universe.writeToFile("1s3v.pdb")

mol_dir = '/home/lspirido/AstexDiv_xtal/1-build/1s3v/'

integrator = VelocityVerletIntegrator(universe, delta_t = 1.*MMTK.Units.fs)
GCintegrator = GCHMCIntegrator(universe, mol_dir,  '/export/apps/amber/14/dat/leap/parm/')

allconfs = []
for i in range(10000):
  (confs, energies, H, acc, delta_t) = GCintegrator.Call(50, 5, 300, 0.0015, random.randint(0, 254))
  (confs, energies, H, acc, delta_t) = integrator(steps=20, T=300, steps_per_trial =1)
  allconfs.append(confs[-1])
#

traj_FN = "1s3v.dcd"
IO_dcd = AlGDock.IO.dcd(molecule)
IO_dcd.write(traj_FN, allconfs)

GCintegrator.Clear()
