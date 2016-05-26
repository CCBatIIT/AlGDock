import os
import sys
from copy import deepcopy

from MMTK import *
from MMTK.ForceFields import Amber12SBForceField
from MMTK.Trajectory import Trajectory, TrajectoryOutput, StandardLogOutput
from MMTK.Dynamics import VelocityVerletIntegrator, TranslationRemover, RotationRemover

import AlGDock.IO
from AlGDock.Integrators.VelocityVerlet.VelocityVerlet import VelocityVerletIntegrator
from PoseFF import PoseForceField

#  if os.environ.get('MMTKDATABASE', 'Not Set'): 
#    os.environ['MMTKDATABASE'] = \
#    '/export/apps/algdock/AlGDock/MMTK/Database/ /home/lspirido/Installers/0Work/0David/PoseFF/Database/'
#
#  mol_dir = '/home/lspirido/Installers/AlGDock-0.0.1.01/Test_1s3vA/prmtopcrd_1s3v/'
#  mol_frcmod = mol_dir + 'ligand.frcmod'

universe = InfiniteUniverse(Amber12SBForceField(parameter_file='../../../../../Example/prmtopcrd/gaff.dat', mod_files=['1s3v.frcmod']))

molecule = Molecule(os.path.basename('1s3v.db'))

universe.addObject(molecule)
universe.configuration();
configuration = universe.configuration().array

# External dofs parameters
extDist = 0.15 # 1.5 A
extAngl = 1.9024 # 109 degrees
extDihe = 0.6108 # 35 degrees
extIs = [0, 1, 2]
for i in range(3): print universe.atomList()[extIs[i]].name,
print
extXYZ = [2.5, 2.5, 2.5] # 20 A
extBK = 1000
extAK = 600
extDK = 500
bbot = 0.1 # 1 A
abot = 0.3490659 # 20 degrees
dbot = 0.523598776 # 30 degrees
intK = 500
poseInp = [[] for i in range(6)]

# I1, I2, I3, X, Y, Z, 0, 0
poseInp[0] = extIs
poseInp[0] = poseInp[0] + extXYZ
poseInp[0] = poseInp[0] + [bbot, abot]
# [b0, kb, theta0, kd, n, gamma = 35 degrees, bottom, kd]
poseInp[1] = [extDist, extBK, extAngl, extAK, 0, extDihe, dbot, extDK]

poseInpLine = [34,13,11,53]
gamma = 1.5708 # 90 deg
poseInpLine = poseInpLine + [0, gamma, dbot, intK]
poseInp[2] = deepcopy(poseInpLine)
for i in range(4): print universe.atomList()[poseInpLine[i]].name,
print

poseInpLine = [45,47,1,48]
gamma = -1.047 # -60 deg
poseInpLine = poseInpLine + [0, gamma, dbot, intK]
poseInp[3] = deepcopy(poseInpLine)
for i in range(4): print universe.atomList()[poseInpLine[i]].name,
print

poseInpLine = [43,18,50,0]
gamma = 2.0943 # 120 deg
poseInpLine = poseInpLine + [0, gamma, dbot, intK]
poseInp[4] = deepcopy(poseInpLine)
for i in range(4): print universe.atomList()[poseInpLine[i]].name,
print

poseInpLine = [20,3,51,4]
gamma = -2.0943 # -120 deg
poseInpLine = poseInpLine + [0, gamma, dbot, intK]
poseInp[5] = deepcopy(poseInpLine)
for i in range(4): print universe.atomList()[poseInpLine[i]].name,
print

print poseInp

poseFF = PoseForceField(poseInp)
universe.setForceField(poseFF)

universe.writeToFile("1s3v.pdb")
integrator = VelocityVerletIntegrator(universe, delta_t = 1.*Units.fs)
(confs, energies, H, acc, delta_t) = integrator(steps=2000, T=300, steps_per_trial=1)

traj_FN = "1s3v.dcd"
IO_dcd = AlGDock.IO.dcd(molecule)
IO_dcd.write(traj_FN, confs)
