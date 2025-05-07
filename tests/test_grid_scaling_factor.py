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

def _load_record(record):
    """
    Note: This function originates from the 'prmtop' class in the AlGDock.IO module.
    """
    items = []
    lines = record.split('\n')
    lines.pop(0)  # Name
    FORMAT = lines.pop(0).strip()[8:-1]  # Format
    if FORMAT.find('a') > -1:  # Text
        w = int(FORMAT[FORMAT.find('a') + 1:])
        for line in lines:
            items = items + [line[x:x + w] for x in range(0, len(line), w)]
        return np.array(items)
    elif FORMAT.find('I') > -1:  # Integer
        w = int(FORMAT[FORMAT.find('I') + 1:])
        for line in lines:
            items = items + [int(line[x:x + w]) for x in range(0, len(line), w)]
        return np.array(items, dtype=int)
    elif FORMAT.find('E') > -1:  # Scientific
        w = int(FORMAT[FORMAT.find('E') + 1:FORMAT.find('.')])
        for line in lines:
            items = items + [float(line[x:x + w]) for x in range(0, len(line), w)]
        return np.array(items, dtype=float)
        
def readPrmtop(FN, varnames=['RESIDUE_LABEL', 'RESIDUE_POINTER']):
    """
    Reads an AMBER prmtop file, returning a dictionary
    Note: This function originates from the 'prmtop' class in the AlGDock.IO module.
    """
    if not os.path.isfile(FN):
        raise Exception('prmtop file %s does not exist!' % FN)
    if FN.endswith('.gz'):
        import gzip
        F = gzip.open(FN, 'r')
    else:
        F = open(FN, 'r')
    data = F.read().split('%FLAG ')
    F.close()
    prmtop = {}
    for record in data:
        name = record[:record.find('\n')].strip()
        if name in varnames:
            prmtop[name] = _load_record(record)
    return prmtop


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

varnames = ['POINTERS','TITLE','ATOM_NAME','AMBER_ATOM_TYPE','CHARGE','MASS',\
            'NONBONDED_PARM_INDEX','LENNARD_JONES_ACOEF','LENNARD_JONES_BCOEF',\
            'ATOM_TYPE_INDEX','BONDS_INC_HYDROGEN','BONDS_WITHOUT_HYDROGEN',\
            'RADII','SCREEN']
try:
    import AlGDock
    import AlGDock.IO
    prmtop_IO = AlGDock.IO.prmtop()
    prmtop_alg = prmtop_IO.read(ligand_prmtop, varnames)
except:
    prmtop_alg = readPrmtop(ligand_prmtop, varnames)
    
NATOM = prmtop_alg['POINTERS'][0]
NTYPES = prmtop_alg['POINTERS'][1]
LJ_radius = np.ndarray(shape=(NTYPES), dtype=float)
LJ_depth = np.ndarray(shape=(NTYPES), dtype=float)
for i in range(NTYPES):
    LJ_index = prmtop_alg['NONBONDED_PARM_INDEX'][NTYPES*i+i]-1
    if prmtop_alg['LENNARD_JONES_ACOEF'][LJ_index]<1.0e-6:
        LJ_radius[i] = 0
        LJ_depth[i] = 0
    else:
        factor = 2 * prmtop_alg['LENNARD_JONES_ACOEF'][LJ_index] / prmtop_alg['LENNARD_JONES_BCOEF'][LJ_index]
        LJ_radius[i] = pow(factor, 1.0/6.0) * 0.5
        LJ_depth[i] = prmtop_alg['LENNARD_JONES_BCOEF'][LJ_index] / 2 / factor
root_LJ_depth = np.sqrt(LJ_depth)
LJ_diameter = LJ_radius*2
atom_type_indicies = [prmtop_alg['ATOM_TYPE_INDEX'][atom_index]-1 for atom_index in range(NATOM)]
scaling_factor_LJr_dict = dict()
for (name,type_index) in zip(prmtop_alg['ATOM_NAME'],atom_type_indicies):
    scaling_factor_LJr_dict[name.strip()] = round(4.184 * root_LJ_depth[type_index]*(LJ_diameter[type_index]**6),6)

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
MMTK_ljr_data = []
for i in range(system.getNumParticles()):
    atom_name = atoms_name[i].strip()
    charge, sigma, epsilon = force.getParticleParameters(i)
    charge_val = charge.value_in_unit(elementary_charge)
    sigma_val = round(sigma.value_in_unit(nanometer),6)
    epsilon = epsilon.value_in_unit(kilojoule/mole)
    if epsilon == 0:
        continue
    MMTK_ele_data.append(ele_data[atom_name])
    MMTK_ljr_data.append(ljr_data[atom_name])
    scaling_factors_ELE.append(round(charge_val*4.184, 6))
    scaling_factors_LJr.append(scaling_factor_LJr_dict[atom_name])



assert scaling_factors_ELE == MMTK_ele_data
assert scaling_factors_LJr == MMTK_ljr_data
