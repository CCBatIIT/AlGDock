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
import AlGDock.RigidBodies

import random

amberFF = MMTK.ForceFields.Amber12SBForceField(\
  parameter_file='../../../Example/prmtopcrd/gaff.dat', mod_files=['1s3v.frcmod'])
universe = MMTK.InfiniteUniverse()

molecule = MMTK.Molecule('1s3v.db')
universe.addObject(molecule)

conf = universe.configuration().array

rb = AlGDock.RigidBodies.identifier(universe, molecule)

poseFF = PoseForceField(rb.poseInp())
compoundFF = amberFF + poseFF
universe.setForceField(compoundFF)

universe.energyTerms()

#
universe.writeToFile("1s3v.pdb")

integrator = VelocityVerletIntegrator(universe, delta_t = 1.*MMTK.Units.fs)

allconfs = []
for i in range(10000):
  (confs, energies, H, acc, delta_t) = integrator(steps=20, T=300, steps_per_trial =1)
  allconfs.append(confs[-1])

traj_FN = "1s3v.dcd"
IO_dcd = AlGDock.IO.dcd(molecule)
IO_dcd.write(traj_FN, allconfs)
