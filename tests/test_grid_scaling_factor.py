import os, sys
import numpy as np

try:
    import MMTK
    from MMTK.ForceFields import Amber12SBForceField
except ImportError:
    print('MMTK not found.')
try:
    from openmm.app import AmberPrmtopFile, AmberInpcrdFile, Simulation, NoCutoff
    from openmm import *
    from openmm.unit import *
except ImportError:
    print('OpenMM not found.')
    openmm = None

def generate_benchmark(ligand_database, gaff_forcefield_file, ligand_frcmod_file):
    MMTK.Database.molecule_types.directory = os.path.dirname(os.path.abspath(ligand_database))
    molecule = MMTK.Molecule(os.path.basename(ligand_database))
    universe = MMTK.Universe.InfiniteUniverse()
    universe.addObject(molecule)
    forcefields_gaff = Amber12SBForceField(parameter_file = gaff_forcefield_file, mod_files = [ligand_frcmod_file])

    fdata = open('benchmark_MMTK_scaling_factors.txt','w')
    fdata.write('name    electrostatic   LJr         index\n')
    for a in molecule.atomList():
        full_name = molecule.getAtomProperty(a, 'name')
        name = full_name.split('i')[0]
        idx = int(full_name.split('i')[1])
        ele = float(molecule.getAtomProperty(a, 'scaling_factor_electrostatic'))
        ljr = float(molecule.getAtomProperty(a, 'scaling_factor_LJr'))
        fdata.write('{:<{}}'.format(name, 4) + "{:15.6f}".format(round(ele,6))+ "{:15.6f}".format(round(ljr,6))+ "{:4d}".format(idx)+'\n')
    fdata.close()
    print('Benchmark data generated successfully.')
    return

ligand_database = 'tests/input/ligand.db'
gaff_forcefield_file = 'tests/input/gaff2.dat'
ligand_frcmod_file = 'tests/input/ligand.frcmod'

ligand_prmtop = 'tests/input/ligand.prmtop'
ligand_inpcrd = 'tests/input/ligand.trans.inpcrd'

#------------  MMTK results -------#
if not os.path.isfile('benchmark_MMTK_scaling_factors.txt'):
    try:
        generate_benchmark(ligand_database, gaff_forcefield_file, ligand_frcmod_file)
    except:
        print('benchmark data generation failed. Terminating process.')
        sys.exit(1)

ele_data = dict()
ljr_data = dict()
idxes = dict()
fdata = open('benchmark_MMTK_scaling_factors.txt','r')
lines = fdata.read().split('\n')
for line in lines[1:]:
    if line != '':
        name = line.split()[0]
        ele = float(line.split()[1])
        ljr = float(line.split()[2])
        idx = int(line.split()[3])
        ele_data[name] = ele
        ljr_data[name] = ljr
        idxes[name] = idx
fdata.close()

#-----------------------------------#
# Use OpenMM to reproduce potentials
#-----------------------------------#
if not openmm:
    print('OpenMM not found. Failed to compare. Terminating process.')
    sys.exit(1)

prmtop = AmberPrmtopFile(ligand_prmtop)
inpcrd = AmberInpcrdFile(ligand_inpcrd)
topology = prmtop.topology
positions = inpcrd.positions # angstrom

r_distance = dict()
for bond in topology.bonds():
    atom1, atom2 = bond
    atom1_position = list(positions[atom1.index].value_in_unit(nanometer))
    atom2_position = list(positions[atom2.index].value_in_unit(nanometer))
    dist = ((atom1_position[0] - atom2_position[0])**2 +
            (atom1_position[1] - atom2_position[1])**2 +
            (atom1_position[2] - atom2_position[2])**2)**0.5
    dist = round(dist,6)
    r_distance[atom1.name] = (atom2.name, dist)
    r_distance[atom2.name] = (atom1.name, dist)

    
atoms_name = [a.name for a in topology.atoms()]
atoms_elements = dict()
for a in topology.atoms():
    atoms_elements[a.name] = a.element.symbol
    
system = prmtop.createSystem(nonbondedMethod=NoCutoff, nonbondedCutoff=1.0*nanometer, constraints=None)
force = system.getForce(3)
assert force.getName() == 'NonbondedForce'

scaling_factors_LJr = []
scaling_factors_ELE = []
MMTK_ele_data = []
for i in range(system.getNumParticles()):
    atom_name = atoms_name[i]
    charge, sigma, epsilon = force.getParticleParameters(i)
    charge_val = charge.value_in_unit(elementary_charge)
    sigma_val = round(sigma.value_in_unit(nanometer),6)
    epsilon = epsilon.value_in_unit(kilojoule/mole)
    if epsilon == 0:
        continue
    r = r_distance[atom_name][1]
    lj_potential = 4 * epsilon * ((sigma_val / r) ** 12 - (sigma_val / r) ** 6)
    #print(atoms_name[i], '4 *',epsilon, '* ((',sigma_val,'/r ))**12 -', sigma_val, '/r ))**6 =',ljr_data[atom_name])
    #print(atoms_name[i], round(charge_val*4.184, 6), ele_data[atom_name])
    MMTK_ele_data.append(ele_data[atom_name])

    scaling_factors_LJr.append(lj_potential)
    scaling_factors_ELE.append(round(charge_val*4.184, 6))

assert scaling_factors_ELE == MMTK_ele_data

#TODO: Need to figure out scaling_factors_LJr
#left of the equal sign: openmm items
#right of the equal sign: mmtk value
#r is the estimation
#H1 4 * 0.06276000026869928 * (( 0.259964 /r ))**12 - 0.259964 /r ))**6 = 316.3363
#r = 0.1421
#H2 4 * 0.06276000026869928 * (( 0.259964 /r ))**12 - 0.259964 /r ))**6 = 316.3363
#r = 0.1421
#H3 4 * 0.06276000026869928 * (( 0.259964 /r ))**12 - 0.259964 /r ))**6 = 316.3363
#r = 0.1421
#H4 4 * 0.06276000026869928 * (( 0.259964 /r ))**12 - 0.259964 /r ))**6 = 316.3363
#r = 0.1421
#C1 4 * 0.35982400053705343 * (( 0.339967 /r ))**12 - 0.339967 /r ))**6 = 3788.707674
#r = 0.1762
#C2 4 * 0.35982400053705343 * (( 0.339967 /r ))**12 - 0.339967 /r ))**6 = 3788.707674
#r = 0.1762
#C3 4 * 0.35982400053705343 * (( 0.339967 /r ))**12 - 0.339967 /r ))**6 = 3788.707674
#r = 0.1762
#C4 4 * 0.35982400053705343 * (( 0.339967 /r ))**12 - 0.339967 /r ))**6 = 3788.707674
#r = 0.1762
#C5 4 * 0.35982400053705343 * (( 0.339967 /r ))**12 - 0.339967 /r ))**6 = 3788.707674
#r = 0.1762
#C6 4 * 0.35982400053705343 * (( 0.339967 /r ))**12 - 0.339967 /r ))**6 = 3788.707674
#r = 0.1762
#C7 4 * 0.35982400053705343 * (( 0.339967 /r ))**12 - 0.339967 /r ))**6 = 3788.707674
#r = 0.1762
#O1 4 * 0.8786399993381852 * (( 0.295992 /r ))**12 - 0.295992 /r ))**6 = 2578.771323
#r = 0.1706
#O2 4 * 0.8786399993381852 * (( 0.295992 /r ))**12 - 0.295992 /r ))**6 = 2578.771323
#r = 0.1706
#O3 4 * 0.8803136010405576 * (( 0.306647 /r ))**12 - 0.306647 /r ))**6 = 3191.388968
#r = 0.1737
