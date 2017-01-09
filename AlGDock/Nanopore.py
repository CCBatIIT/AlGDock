#!/usr/bin/env python

Simulation_arguments = {
  # Input files
  'ligand_tarball':{'help':'tarball with the other ligand files inside'},
  'ligand_database':{'help':'MMTK molecule definition for ligand'},
  'forcefield':{'help':'AMBER force field file'},
  'frcmodList':{'help':'AMBER force field modifications file(s)', 'nargs':'+'},
  'ligand_prmtop':{'help':'AMBER prmtop for the ligand'},
  'grid_LJr':{'help':'DX file for Lennard-Jones repulsive grid'},
  'grid_LJa':{'help':'DX file for Lennard-Jones attractive grid'},
  'grid_ELE':{'help':'DX file for electrostatic grid'},
  'starting_conf':{'help':'Starting configuration in inpcrd or mol2 format'},
  # Input parameters
  'ef':{'help':'Electric field strength'},
  'max_trials':{'help':'Maximum number of Hamiltonian Monte Carlo trials'},
  'report_interval':{'help':'Number of Hamiltonian Monte Carlo trials' + \
    ' between saved configurations'},
  # Output files
  'output_dcd':{'help':'data file to store configurations'}
}

necessary_FN_keys = [\
  'ligand_database','forcefield','ligand_prmtop',\
  'grid_LJr','starting_conf']

# 'grid_LJa','grid_ELE'

###########
# Imports #
###########

import os
import sys
import time
import numpy as np

import MMTK
import MMTK.Units
from MMTK.ParticleProperties import Configuration
from MMTK.ForceFields import ForceField

import Scientific
try:
  from Scientific._vector import Vector
except:
  from Scientific.Geometry.VectorModule import Vector

R = 8.3144621 * MMTK.Units.J / MMTK.Units.mol / MMTK.Units.K

####################
# Simulation class #
####################

class Simulation:
  def __init__(self, **args):
    # Check and store arguments, setting undefined values to None
    for key in Simulation_arguments.keys():
      if not key in args.keys():
        args[key] = None
    for key in necessary_FN_keys:
      if args[key] is None:
        raise Exception('No input for necessary file '+key)
    self.args = args

    # Check for necessary files, including in the ligand tarball if provided
    for key in necessary_FN_keys:
      args[key] = os.path.abspath(args[key])
    args['frcmodList'] = [os.path.abspath(FN) for FN in args['frcmodList']]
    
    necessary_FNs = [args[file_type] for file_type in necessary_FN_keys] + \
      args['frcmodList']

    self._toClear = []
    if args['ligand_tarball'] is not None:
      seekFNs = [FN for FN in necessary_FNs if not os.path.isfile(FN)]
      import tarfile
      tarF = tarfile.open(args['ligand_tarball'],'r')
      for member in tarF.getmembers():
        for seekFN in seekFNs:
          if member.name.endswith(os.path.basename(seekFN)):
            tarF.extract(member)
            self._toClear.append(os.path.abspath(seekFN))
            print '  extracted '+seekFN

    for (FN,file_type) in zip(necessary_FNs,necessary_FN_keys+['ligand_frcmod']):
      if not os.path.isfile(FN):
        raise Exception('Necessary file for %s is not at %s'%(file_type,FN))
        
    if self.args['output_dcd'] is None:
      self.args['output_dcd'] = 'output.dcd'

    self._setup_universe()

  def _setup_universe(self):
    # Create molecule object
    MMTK.Database.molecule_types.directory = \
      os.path.dirname(self.args['ligand_database'])
    self.molecule = MMTK.Molecule(\
      os.path.basename(self.args['ligand_database']))

    # Hydrogen mass repartitioning to 4.0 amu
    from AlGDock.HMR import hydrogen_mass_repartitioning
    self.molecule = hydrogen_mass_repartitioning(self.molecule, 4.0)

    # Helpful variables for referencing and indexing atoms in the molecule
    self.molecule.heavy_atoms = [ind for (atm,ind) in \
      zip(self.molecule.atoms,range(self.molecule.numberOfAtoms())) \
      if atm.type.name!='hydrogen']

    self.molecule.prmtop_atom_order = np.array([atom.number \
      for atom in self.molecule.prmtop_order], dtype=int)
    self.molecule.inv_prmtop_atom_order = np.zeros(shape=len(self.molecule.prmtop_atom_order), dtype=int)
    for i in range(len(self.molecule.prmtop_atom_order)):
      self.molecule.inv_prmtop_atom_order[self.molecule.prmtop_atom_order[i]] = i
    
    # Create universe and add molecule to universe
    self.universe = MMTK.Universe.InfiniteUniverse()
    self.universe.addObject(self.molecule)

    # Determine the net charge on the ligand
    net_charge = 0.
    for o in self.universe:
      for a in o.atomList():
        net_charge += float(o.getAtomProperty(a, 'scaling_factor_electrostatic'))
    net_charge = net_charge/4.184
    if abs(net_charge)<0.1:
      raise Exception('The net charge on the ligand is too low' + \
        ' for the electric field to have an effect')
    
    # Force fields
    self._forceFields = {}
    from MMTK.ForceFields import Amber12SBForceField
    self._forceFields['gaff'] = Amber12SBForceField(
      parameter_file=self.args['forcefield'],mod_files=self.args['frcmodList'])

    from AlGDock.ForceFields.OBC.OBC import OBCForceField
    self._forceFields['OBC'] = OBCForceField(self.args['ligand_prmtop'],
      self.molecule.prmtop_atom_order,self.molecule.inv_prmtop_atom_order)

    from AlGDock.ForceFields.Grid.Interpolation import InterpolationForceField
    for grid_type in ['LJa','LJr','ELE']:
      if (self.args['grid_'+grid_type] is not None) and \
          os.path.isfile(self.args['grid_'+grid_type]):
        self._forceFields[grid_type] = InterpolationForceField(\
          self.args['grid_'+grid_type], name=grid_type, \
          interpolation_type='Trilinear', strength=1.0, \
          scaling_property='scaling_factor_' + \
            {'LJr':'LJr','LJa':'LJa','ELE':'electrostatic'}[grid_type], \
          inv_power=4 if grid_type=='LJr' else None)
      else:
        print '  grid type %s not found'%grid_type

    # Set electric field strength to move the ligand to larger z
    # If the net charge is positive, then a negative electric field will favor larger z
    self.args['ef'] = abs(self.args['ef'])
    if net_charge>0:
      self.args['ef'] = -1*self.args['ef']
      
    from AlGDock.ForceFields.ElectricField.ElectricField import ElectricField
    self._forceFields['electric_field'] = ElectricField(\
      self.args['ef']*(MMTK.Units.V/MMTK.Units.m), 'scaling_factor_electrostatic')
    # scaling_factor_electrostatic is amber_charge multiplied by 4.184, 
    # converting kcal/mol (AMBER units) to kJ/mol (MMTK Units).
    # 1 V m−1 = 1 kg m s−3 A−1.      
    
    FFkeys = self._forceFields.keys()
    compoundFF = self._forceFields[FFkeys[0]]
    for FFkey in FFkeys:
      compoundFF += self._forceFields[FFkey]
    self.universe.setForceField(compoundFF)
    
    # Set the ligand starting coordinates
    # Load the file
    import AlGDock.IO
    if self.args['starting_conf'].endswith('.inpcrd'):
      reader = AlGDock.IO.crd()
      lig_crd = IO_crd.read(self.args['starting_conf'], multiplier=0.1)
    elif self.args['starting_conf'].endswith('.mol2') or \
       self.args['starting_conf'].endswith('.mol2.gz'):
      reader = AlGDock.IO.dock6_mol2()
      lig_crd = reader.read(self.args['starting_conf'])[0][0]
    else:
      raise Exception('Unknown file extension')
    lig_crd = lig_crd[self.molecule.inv_prmtop_atom_order,:]
    
    # Randomly rotate the ligand
    from AlGDock.Integrators.ExternalMC.ExternalMC import random_rotate
    lig_crd = np.dot(lig_crd, np.transpose(random_rotate()))
    self.universe.setConfiguration(Configuration(self.universe,lig_crd))
    
    # Translate the ligand to the middle of the grid on the x and y axes 
    # and 10% of the distance on the z axis
    gd = None
    for grid_type in ['ELE','LJr','LJa']:
      if grid_type in self._forceFields.keys():
        gd = self._forceFields[grid_type].grid_data
    if gd is None:
      raise Exception('No grids loaded')
    starting_position = np.copy(gd['origin'])
    starting_position[-1] = (gd['counts'][-1]*gd['spacing'][-1]/10.)
    starting_position[:2] += (gd['counts']*gd['spacing']/2)[:2]
    self.universe.translateTo(MMTK.Vector(starting_position))
    self.z_i = starting_position[-1]
    self.z_f = gd['origin'][-1] + (gd['counts'][-1]*gd['spacing'][-1]*9./10.)
  
  def run(self):
    """
    Performs a molecular dynamics until the system crosses 
    the desired position on the z axis.
    """
    print 'Simulating motion from %.5f to %.5f'%(self.z_i, self.z_f)

    confs = []
    for n in range(self.args['max_trials']):
      from AlGDock.Integrators.HamiltonianMonteCarlo.HamiltonianMonteCarlo \
        import HamiltonianMonteCarloIntegrator
      sampler = HamiltonianMonteCarloIntegrator(self.universe)
      sampler(T=300.0*MMTK.Units.K, steps=100)
      if (n%100==0):
        print 'Trial %6d, COM:'%n, self.universe.centerOfMass()
        confs.append(np.copy(self.universe.configuration().array))
      if self.universe.centerOfMass()[-1]>self.z_f:
        break
    print 'Number of trials before pore crossing: ', n
    print 'Number of configurations: ', len(confs)
    
    dcdFN = 'test.dcd'

    import AlGDock.IO
    IO_dcd = AlGDock.IO.dcd(self.molecule,
      ligand_atom_order = self.molecule.prmtop_atom_order)
    IO_dcd.write(self.args['output_dcd'], confs)
    
    self.confs = confs

  def view_trajectory(self):
    """
    Output the trajectory as a DCD file and views it in VMD
    """
    if self.args['grid_LJr'].endswith('.nc') and (not os.path.isfile(self.args['grid_LJr'][:-2]+'dx')):
      import AlGDock.IO
      IO_Grid = AlGDock.IO.Grid()
      IO_Grid.read(self.args['grid_LJr'])
      IO_Grid.write(self.args['grid_LJr'][:-2]+'dx')
    
    scriptFN = 'script.vmd'

    script = ''
    script +=  'set ligand [mol new '+self.args['ligand_prmtop']+']\n'
    script += 'mol addfile '+self.args['output_dcd']+' type dcd waitfor all\n'
    script += 'mol new {%s} type {dx} first 0 last -1 step 1 waitfor 1\n'%(self.args['grid_LJr'][:-2]+'dx')
    script += 'mol modstyle 0 1 Isosurface 10.000000 0 2 0 1 1\n'

    scriptF = open(scriptFN,'w')
    scriptF.write(script)
    scriptF.close()

    import AlGDock
    vmdFN = AlGDock.findPaths(['vmd'])['vmd']
    vmd_args = [vmdFN] + ['-e', scriptFN]
    import subprocess
    p = subprocess.call(vmd_args)

    # p = subprocess.Popen(vmd_args, \
    #   stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # (vmd_stdout, vmd_stderr) = p.communicate()
    # p.wait()

  def __del__(self):
    if len(self._toClear)>0:
      print '\n>>> Clearing files'
      for FN in self._toClear:
        if os.path.isfile(FN):
          os.remove(FN)
          print '  removed '+FN

if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(
    description='Simulation of a small molecule in a nanopore stochastic sensor')
  
  for key in Simulation_arguments.keys():
    parser.add_argument('--'+key, **Simulation_arguments[key])
  args = parser.parse_args()

  self = Simulation(**vars(args))
