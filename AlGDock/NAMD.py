#!/usr/local/bin/env python

"""

A python wrapper for NAMD 2.9

"""

#=============================================================================================
# IMPORTS
#=============================================================================================

import os, sys, gzip

try:
  from AlGDock import findPath
  from AlGDock import search_paths
except:
  def findPath(locations):
    """
      Parses a list of locations, returning the first file that exists.
      If none exist, then None is returned.
      """
    import os.path
    for location in locations:
      if location is not None and os.path.exists(location):
        return os.path.abspath(location)
    if not locations == [None]:
      print 'File not found!'
      print 'Searched:'
      print locations
    return None

  # Define search paths for external programs
  # Defined for David's IIT MacBook Pro, DSCR cluster, and CCB cluster
  search_paths = {
    'namd':['/Users/dminh/Installers/NAMD_2.9_Source/MacOSX-x86_64-g++/namd2',
            '/home/dbchem/dm225/.local/bin/namd2',
            '/share/apps/namd/2.9/Linux-x86_64-g++/namd2',
            None]}

class GRID:
  def __init__(self,
      forceFN, LJR_FN, LJA_FN, ELE_FN,
      lambdaVal=1.0):
    """
    Class to define the grid portion of a NAMD configuration file.
    
    forceFN - a pdb file with the following in the columns
      X - partial charges
      Y - Lennard-Jones repulsive interaction multiplier
      Z - Lennard-Jones attractive interaction multiplier
      O - anything
      B - ones
      
    LJR_FN - a dx file with the Lennard-Jones repulsive potential grid
    LJA_FN - a dx file with the Lennard-Jones attractive potential grid
    ELE_FN - a dx file with the electrostatic potential grid
    
    lambdaVal - the value of lambda
    """
    
    self.set_scale(lambdaVal)
    
    self.LJ_conf = '''
mgridforcefile        LJR {%s}
mgridforcecol         LJR B # Ones
mgridforcechargecol   LJR Y # Lennard-Jones repulsive interactions
mgridforcepotfile     LJR {%s}
mgridforcecont1       LJR no
mgridforcecont2       LJR no
mgridforcecont3       LJR no
mgridforcelite        LJR yes

mgridforcefile        LJA {%s}
mgridforcecol         LJA B # Ones
mgridforcechargecol   LJA Z # Lennard-Jones attractive interactions
mgridforcepotfile     LJA {%s}
mgridforcecont1       LJA no
mgridforcecont2       LJA no
mgridforcecont3       LJA no
mgridforcelite        LJA yes

'''%(forceFN,LJR_FN,forceFN,LJA_FN)
  
    self.ELE_conf = '''
mgridforcefile        ELE {%s}
mgridforcecol         ELE B # Ones
mgridforcechargecol   ELE X # Charges
mgridforcepotfile     ELE {%s}
mgridforcecont1       ELE no
mgridforcecont2       ELE no
mgridforcecont3       ELE no
mgridforcelite        ELE yes

'''%(forceFN,ELE_FN)

  def set_scale(self,lambdaVal=1.0):
    """
    Based on lambdaVal, sets class variables scale_LJ and scale_ELE
    """
    self.lambdaVal = lambdaVal
    if self.lambdaVal==0:
      self.scale_LJ = 0
      self.scale_ELE = 0
    elif self.lambdaVal<0.5:
      self.scale_LJ = self.lambdaVal/0.5
      self.scale_ELE = 0
    else:
      self.scale_LJ = 1
      self.scale_ELE = (self.lambdaVal-0.5)/0.5

  def _scale_conf(self,type,scale):
    """
    Returns a string that describes grid force scaling
    """
    return 'mgridforcescale       %3s %8.6e %8.6e %8.6e\n'%(type,scale,scale,scale)

  def script(self):
    """
    Returns the grid force portion of a NAMD configuration script, with scaling according to lambdaVal.
    """
    conf = ''
    if self.lambdaVal > 0:
      conf = conf + 'mgridforce            on\n' + self.LJ_conf
      conf = conf + self._scale_conf('LJR',self.scale_LJ) + self._scale_conf('LJA',self.scale_LJ)
    if self.lambdaVal > 0.5:
      conf = conf + self.ELE_conf + self._scale_conf('ELE',self.scale_ELE)
    return conf

  def script_LJ(self):
    """
    Returns a fully scaled Lennard-Jones grid force portion of a NAMD configuration script
    """
    return 'mgridforce            on\n' + self.LJ_conf + self._scale_conf('LJR',1) + self._scale_conf('LJA',1)

  def script_ELE(self):
    """
    Returns a fully scaled electrostatic grid force portion of a NAMD configuration script
    """
    return 'mgridforce            on\n' + self.ELE_conf + self._scale_conf('ELE',1)

class COLVARS_BINDING_SITE_SPHERE:
  def __init__(self,
      radius,
      atomFN, atomsCol='B', atomsColValue=1.0,
      center=[0.0,0.0,0.0], centerFN=None,
      colvarsFN='sphere.colvars'):
    """
    Class containing collective variables script and configuration file options for a spherical binding site.
    """

    self.radius = radius
    self.center = center
    self.colvarsFN = colvarsFN

    if not (centerFN==None):
      centerF = open(centerFN,'r')
      center = [float(item) for item in centerF.read().strip().split(' ')]
      self.center = center
      centerF.close()
    
    self.colvars_script = '''
colvar {
  name DistanceFromSite
  upperWall %.2f
  upperWallConstant 10.0
  distance {
    group1 {
      atomsFile {%s}
      atomsCol {%s}
      atomsColValue {%.2f}
    }
    group2 {
      dummyAtom (%8.6f, %8.6f, %8.6f)
    }
  }
}
'''%(self.radius,atomFN,atomsCol,atomsColValue,self.center[0],self.center[1],self.center[2])

    self.conf = '\ncolvars              on\ncolvarsConfig        {%s}\n'%colvarsFN
    
  def write(self):
    if not os.path.exists(self.colvarsFN):
      colvarsF = open(self.colvarsFN,'w')
      colvarsF.write(self.colvars_script)
      colvarsF.close()

class COLVARS_BINDING_SITE_CYLINDER:
  def __init__(self,
      radius,
      minZ,
      maxZ,
      atomFN, atomsCol='B', atomsColValue=1.0,
      axis=[0.0,0.0,0.0,1.0,0.0,0.0], axisFN=None,
      colvarsFN='cylinder.colvars'):
    """
    Class containing collective variables script and configuration file options for a cylindrical binding site.
    """

    self.radius = radius
    self.minZ = minZ
    self.maxZ = maxZ
    
    self.colvarsFN = colvarsFN

    if not (axisFN==None):
      axisF = open(axisFN,'r')
      axis = [float(item) for item in axisF.read().strip().split(' ')]
      axisF.close()
      self.axis = axis
    
    groups = '''
    main {
      atomsFile {%s}
      atomsCol %s
      atomsColValue %d
    }
    ref {
      dummyAtom (%.6f, %.6f, %.6f)
    }
    axis (%.6f, %.6f, %.6f)
    '''%(atomFN,atomsCol,atomsColValue,axis[0],axis[1],axis[2],axis[3],axis[4],axis[5])
    
    self.colvars_script = '''
colvar {
  name DistXY
  upperWall %.2f
  upperWallConstant 10.0
  distanceXY { %s }
}

colvar {
  name DistZ
  lowerWall %.2f
  upperWall %.2f
  lowerWallConstant 10.0
  upperWallConstant 10.0
  distanceZ { %s }
}
'''%(self.radius,groups,self.minZ,self.maxZ,groups)

    self.conf = '\ncolvars              on\ncolvarsConfig        {%s}\n'%colvarsFN
    
  def write(self):
    if not os.path.exists(self.colvarsFN):
      colvarsF = open(self.colvarsFN,'w')
      colvarsF.write(self.colvars_script)
      colvarsF.close()

class NAMD:
  def __init__(self,
      namd_command='namd2',NPROCS=1,
      prmtop=None, inpcrd=None, bincoordinates=None, binvelocities=None,
      xsc=None, fixed=None,
      solvent='Gas', useCutoff=True, grid=None, colvars=None, alchemical=None,
      seed=None, finishBy=None,
      debug=False):
    """
    Wrapper for NAMD 2.9
    
    Required arguments:
    namd_command - the location of NAMD 2.9 [Default is 'namd2']
    NPROCS - the number of processors [Default is 1]
    prmtop - AMBER parameters and topology for the system
    inpcrd - AMBER coordinates for the system
    
    Optional arguments:
    bincoordinates - binary coordinates to restart a NAMD run [Default is None]
    binvelocities - binary velocities to restart a NAMD run [Default is None]
    xsc - extended system file [Default is None]
    fixed - pdb file with fixed atoms labeled with 1 in the occupancy column
    solvent - either 'TIP3P', 'GBSA', 'GB', or 'Gas' [Default is Gas]
      TIP3P uses a shorter cutoff than the others
      GBSA and GB use Generalized Born solvation
      GBSA also includes a surface area term      
    useCutoff - use a cutoff for nonbonded terms
    grid - GRID class, or None [Default is None]
    colvars - class containing conf and colvars objects, or None [Default is None]
    alchemical - ALCHEMICAL class (Doesn't do anything yet)
    seed - a random number seed to start NAMD simulation [Default is None]
    finishBy - the time (in seconds) any NAMD instance should be finished by [Default is none, meaning that there is no time limit.]
    debug - does not remove any files [Default = False]
    """
    
    self.namd_command = namd_command
    self.NPROCS = NPROCS
    self.prmtop = prmtop
    self.inpcrd = inpcrd
    self.bincoordinates = bincoordinates
    self.binvelocities = binvelocities
    self.xsc = xsc
    self.fixed = fixed
    self.solvent = solvent
    self.useCutoff = useCutoff
    self.grid = grid
    self.colvars = colvars
    self.alchemical = alchemical
    self.seed = seed
    self.finishBy = finishBy
    self.debug = debug
    
    if self.prmtop==None:
      raise Exception("AMBER prmtop file required")
    if self.inpcrd==None:
      raise Exception("AMBER inpcrd file required")

    # Find NAMD
    self.namd_command = findPath([self.namd_command] + search_paths['namd'])

    if self.namd_command==None:
      raise Exception("NAMD not found!")

  def _readEnergyDatGZ(self, file):
    """
    Reads energies from a gzip file
    """
    energyF = gzip.open(file,'r')
    lines = energyF.read().strip().split('\n')
    energyF.close()
    energies = []
    for line in lines:
      energies.append([float(item) for item in line.split('\t')])
    return energies

  def _writeEnergyDatGZ(self, file, energies):
    """
    Writes energies to a gzip file
    """
    energyF = gzip.open(file,'w')
    for line in energies:
      energyF.write('\t'.join('%.4f'%item for item in line)+'\n')
    energyF.close()

  def _removeFile(self,filename):
    """
    If the file exists and we are not debugging, remove it
    """
    if os.path.exists(filename) and (not self.debug):
      os.remove(filename)
  
  def _removeFiles(self,searchString,forceRemove=False):
    if (not self.debug) or forceRemove:
      import glob
      list = glob.glob(searchString)
      for FN in list:
        os.remove(FN)

  def _writeConfiguration(self, outputname, temperature,
      integrator_script, output_script, execution_script,
      grid_script=''):
    """
    Writes a NAMD configuration script based on class variables
    """
    conf = '''# Variables
set outputname       %s
set temperature      %d

# Input files
amber                on
parmfile             {%s}
ambercoor            {%s}
'''%(outputname,temperature,self.prmtop,self.inpcrd)

    if not (self.bincoordinates==None):
      conf = conf + 'bincoordinates       {%s}\n'%self.bincoordinates
    if not (self.binvelocities==None):
      conf = conf + 'binvelocities        {%s}\n'%self.binvelocities
    if not (self.xsc==None):
      conf = conf + 'extendedSystem       {%s}\n'%self.xsc

    conf = conf + '''
# Force field parameters
exclude              scaled1-4
1-4scaling           0.833333  # =1/1.2, default for AMBER.  NAMD default is 1.0
scnb                 2  # This is AMBER and NAMD default
switching            on
'''

    if not self.useCutoff:
      conf = conf + 'cutoff              999.0\nswitchdist          999.0\npairlistdist        999.0\n\n'
    elif self.solvent=='TIP3P':
      conf = conf + 'cutoff               10.0\nswitchdist            9.0\npairlistdist         11.0\n\n'
    else:
      conf = conf + 'cutoff               16.0\nswitchdist           15.0\npairlistdist         18.0\n\n'

    if self.solvent=='GB' or self.solvent=='GBSA':
      conf = conf + 'GBIS                 on\nionConcentration     0.0\n'
    if self.solvent=='GBSA':
      conf = conf + 'sasa                 on\nsurfaceTension       0.006\n'

    conf = conf + grid_script

    if self.colvars is not None:
      conf = conf + self.colvars.conf

    if self.fixed is not None:
      conf = conf + '''
# Do not calculate fixed-atom energies
fixedAtoms           on
fixedAtomsFile       '''+self.fixed+'''
fixedAtomsCol        O
'''

    if self.binvelocities is None:
      tLine = 'temperature          $temperature'
    else:
      tLine = ''

    conf = conf + '''
# Temperature control
%s
langevin             on
langevinDamping      1
langevinTemp         $temperature
langevinHydrogen     off ;# Don't couple bath to hydrogens
'''%tLine

    if not (self.seed==None):
      conf = conf + 'seed                 %d\n'%self.seed

    conf = conf + integrator_script

    conf = conf + '''
# Output parameters
outputName           $outputname
binaryoutput         yes
'''+output_script

    conf = conf + execution_script

    confF = open(outputname+'.namd','w')
    confF.write(conf)
    confF.close()

  def _execute(self,
      outputname, temperature,
      integrator_script, output_script, execution_script, grid_script='',
      energyFields=[12], writeEnergyDatGZ=False,
      keepScript=False, keepOutput=False, keepCoor=False,
      prmtop=None, inpcrd=None, bincoordinates=None, binvelocities=None,
      xsc=None, solvent=None, grid=None, colvars=None, alchemical=None,
      seed=None, finishBy=None, totalSteps=None,
      debug=None, retry=True):
    """
    Executes a NAMD instance, returning the energies.
    
    Required arguments:
    outputname - prefix for NAMD output
    temperture - simulation temperature
    integrator_script - part of the NAMD configuration file that defines the integration. [Default: '']
    output_script - part of the NAMD configuration file that dictates output.  [Default: '']
    execution_script - part of the NAMD configuration file that determines how the simulation is executed. [Default: '']
    
    Optional Arguments:
    grid_script - part of the NAMD configuration file that defines the interaction grids.  [Default: '']
    energyFields - a list of fields to keep from the ENERGY output lines [Default: [12], which is the total potential energy]
    writeEnergyDatGZ - writes the energies into outputname.dat.gz
    
    finishBy - the time, in seconds, by which the MD simulation should be complete.  If it is defined and totalSteps is defined, NAMD will abort if the projected simulation length (for totalSteps) is longer than the allotted time.  [Default is none, meaning that there is no time limit.]
    totalSteps - the total number of simulation steps.  Only relevant if finishBy is defined.

    keepScript - keeps the NAMD configuration script
    keepOutput - keeps the NAMD output file
    keepCoor - keeps the final coordinate set from the simulation
    
    All optional arguments in the initialization function.
    """
    
    # Parse the arguments
    if not (prmtop==None):
      self.prmtop = prmtop
    if not (inpcrd==None):
      self.inpcrd = inpcrd
    if not (bincoordinates==None):
      self.bincoordinates = bincoordinates
    if not (binvelocities==None):
      self.binvelocities = binvelocities
    if not (xsc==None):
      self.xsc = xsc
    if not (solvent==None):
      self.solvent = solvent
    if not (grid==None):
      self.grid = grid
    if not (colvars==None):
      self.colvars = colvars
    if not (alchemical==None):
      self.alchemical = alchemical
    if not (seed==None):
      self.seed = seed
    if not (finishBy==None):
      self.finishBy = finishBy
    if not (debug==None):
      self.debug = debug
    
    del prmtop, inpcrd, bincoordinates, binvelocities, xsc
    del solvent, seed, finishBy, debug

    # Check to make sure remaining arguments make sense
    if execution_script==None:
      raise Exception("Execution script required")

    original_dir = os.getcwd()
    execution_dir = os.path.dirname(outputname)
    if execution_dir != '':
      os.chdir(execution_dir)
    outputname = os.path.basename(outputname)
    
    # Checks that the process isn't already complete
    if os.path.exists(outputname+'.coor'):
      return []
  
    # Writes the configuration file
    self._writeConfiguration(outputname, temperature, integrator_script, output_script,execution_script, grid_script)

    # Gets the original parameters in the input scripts
    def getParm(keyword,config):
      ind = config.find(keyword)
      if not (ind==-1):
        val = config[ind+len(keyword):]
        if val.find('\n')>-1:
          val = val[:val.find('\n')]
        if val.find(';')>-1:
          val = val[:val.find(';')]
        return float(val)
      return 0.0

    original = {}
    original['timestep'] = getParm('timestep',integrator_script)
    original['dcdfreq'] = getParm('dcdfreq',output_script)
    original['outputEnergies'] = getParm('outputEnergies',output_script)
    original['run'] = getParm('run',execution_script) 
 
    ##############################
    # Executes the NAMD instance #
    ##############################
    
    # Restarts with a new random number seed if necessary
    attempts = 0
    noRunError = True
    while (not os.path.exists(outputname+'.coor')) and noRunError:
      energies = []
      if (not self.finishBy==None) and (not totalSteps==None):
        import time
        startTime = time.time()
        elapsedTimeSteps = 0
        timePerStep = 0

      # Start NAMD
      import subprocess
      proc = subprocess.Popen([self.namd_command,'+p','%d'%self.NPROCS,outputname+'.namd'],stdout=subprocess.PIPE)
    
      outF = open(outputname+'.out','w')
      for line in iter(proc.stdout.readline,''):
        outF.write(line)
        # Store energy output
        if line.find('ENERGY:')==0:
          ENERGY = line[9:].split()
          energyLine = []
          for energyField in energyFields:
            energyLine.append(float(ENERGY[energyField]))
          energies.append(energyLine)
        # End NAMD if there is an error
        if line.find('ERROR:')>-1:
          if line.find('velocity')>-1:
            print 'Atoms moving too fast'
            break
          else:
            print line
            noRunError = False
            raise Exception('Error in NAMD')
          try:
            proc.kill() # After python 2.6
          except AttributeError:
            print 'NAMD not killed.'
        # If there is a time limit, end NAMD is there is insufficient time to complete the instance
        if (not self.finishBy==None) and (not totalSteps==None):
#          if (line.find('ENERGY:')==0) or (line.find('Benchmark time')>-1):
#            if (line.find('ENERGY:')==0):
#              elapsedTimeSteps = int(line[7:15])
#              print '%d / %d steps completed'%(elapsedTimeSteps,totalSteps)
#              if elapsedTimeSteps > 0:
#                timePerStep = (time.time()-startTime)/elapsedTimeSteps
#                print 'Observed rate of %.7f s/step'%timePerStep
#            elif line.find('Benchmark time')>-1:
          if line.find('Benchmark time')>-1:
            timePerStep = float(line[line.find('CPUs')+4:line.find('s/step')])
            print 'Benchmark rate of %.7f s/step'%timePerStep
            projectedCompletion = (totalSteps-elapsedTimeSteps)*timePerStep
            remainingTime = self.finishBy - time.time()
            print 'Projected completion in %3.2f s. %3.2f s remaining.'%(projectedCompletion,remainingTime)
            if projectedCompletion>remainingTime:
              print 'Insufficient time remaining for cycle.'
              outF.write('Insufficient time remaining for cycle.\n')
              noRunError = False
              try:
                proc.kill() # After python 2.6.  Will raise error otherwise.
              except:
                print 'NAMD not terminated normally.'
                sys.exit()
      if noRunError:
        proc.wait() # Let NAMD finish
      outF.close()
      # Restart NAMD if there was only a velocity error
      if (not os.path.exists(outputname+'.coor')) and noRunError:
        if not retry:
          raise Exception('Error in NAMD')
          break
        attempts = attempts + 1
        if attempts <= 5:
          retry_string = 'Retrying'
          if not (self.seed==None):
            self.seed = self.seed + 1
            self._writeConfiguration(outputname, temperature,
              integrator_script, output_script, execution_script, grid_script)
            retry_string = retry_string + ' with a random number seed of %d'%self.seed
          print retry_string
        elif attempts <= 10:
          retry_string = 'Retrying'
          attempts_ts = attempts-5
          new_integrator_script = integrator_script.replace('timestep','timestep %.4f ;#'%(original['timestep']/attempts_ts))
          new_output_script = output_script.replace('dcdfreq','dcdfreq %d ;#'%(original['dcdfreq']*attempts_ts)).replace('outputEnergies','outputEnergies %d ;#'%(original['dcdfreq']*attempts_ts))
          new_execution_script = execution_script.replace('run','run %d ;#'%(original['run']*attempts_ts))
          if not (self.seed==None):
            self.seed = self.seed + 1
            retry_string = retry_string + ' with a random number seed of %d'%self.seed
          self._writeConfiguration(outputname, temperature,
            new_integrator_script, new_output_script, new_execution_script)
          retry_string = retry_string + ' and a time step of %.4f'%(original['timestep']/attempts_ts)
          print retry_string
        else:
          print 'Too many attempts!'
          noRunError = False

    # Clean up
    if not keepScript:
      self._removeFile(outputname+'.namd')
    # Even if keepOutput=False, keep the output if there are run errors
    if (not keepOutput) or (not noRunError):
      self._removeFile(outputname+'.out')
    if not keepCoor:
      self._removeFile(outputname+'.coor')
      self._removeFile(outputname+'.vel')
    self._removeFile(outputname+'.xsc')
    self._removeFile(outputname+'.fep')
    self._removeFile(outputname+'.colvars.traj')
    self._removeFile(outputname+'.colvars.state')
    self._removeFile(outputname+'.colvars.state.old')
    # Return to original directory
    if execution_dir != '':
      os.chdir(original_dir)
    if noRunError:
      if writeEnergyDatGZ:
        self._writeEnergyDatGZ(outputname+'.dat.gz',energies)
      return energies
    else:
      return None

  def mintherm(self, outputname, keepScript=False):
    """
    Minimizes a configuration and thermalizes it to 300 K.
    """
    
    integrator_script = '''
# Integrator parameters
rigidBonds           none
timestep             1.0
nonbondedFreq        2
fullElectFrequency   4
stepspercycle        20
'''

    output_script = '''
outputEnergies       500 ;# 1 ps'''

    execution_script = '''
minimize 1000
for { set curTemp 10 } { $curTemp <= $temperature } { incr curTemp 10 } {
  reinitvels           $curTemp
  langevinTemp         $curTemp
  run                  100
}
'''

    energies = self._execute(outputname,300.0,
        integrator_script, output_script,execution_script,
        keepScript=keepScript)
    return energies

  def simulate(self, outputname, temperature=300.0, steps=10000000, \
    energyFields=[12], \
    keepScript=False, keepCoor=False, writeEnergyDatGZ=False):
    """
    Runs an MD simulation.
    """
    
    integrator_script = '''
# Integrator parameters
rigidBonds           none
timestep             1.0
nonbondedFreq        2
fullElectFrequency   4
stepspercycle        20
'''

    output_script = '''
dcdfreq             1000 ;# 1 ps
outputEnergies      1000 ;# 1 ps
'''

    execution_script = '''
reinitvels           $temperature
run %d
'''%steps

    energies = self._execute(outputname, temperature,
        integrator_script, output_script, execution_script,
        energyFields=energyFields,
        keepScript=keepScript, keepCoor=keepCoor,
        writeEnergyDatGZ=writeEnergyDatGZ)
    return energies

  def _energy_scripts(self,dcdname,stride=1,test=False):
    """
    Scripts for calculating energies
    """
    
    integrator_script = '''
# Integrator parameters
rigidBonds           none
nonbondedFreq        1
fullElectFrequency   1
stepspercycle        1 
'''

    output_script = '''
outputEnergies       1
'''
    if not test:
      execution_script = '''
set ts 0
coorfile open dcd {'''+dcdname+'''}
while { ![coorfile read] } {
  if { [expr $ts %'''+'''%d == 0] } { 
    firstTimestep $ts
    run 0
  }
  incr ts 1
}
coorfile close
'''%stride
    else:
      execution_script = '''run 0'''

    return (integrator_script,output_script,execution_script)

  def energies_PE(self, outputname, dcdname=None, energyFields=[12], \
      stride=1, keepScript=False, writeEnergyDatGZ=True, test=False):
    """
    Calculates potential energies in a dcd file.
    
    outputname - the prefix for the resulting energy file
    dcdname - the dcd file file read [Default - outputname.dcd]
    energyFields - the NAMD energy fields to keep [Default 12, total potential energy]
    """
    
    if dcdname==None:
      dcdname = outputname+'.dcd'
    
    (integrator_script,output_script,execution_script) = self._energy_scripts(dcdname,stride=stride,test=test)
  
    if (os.path.exists('%s.dat.gz'%outputname)):
      energies = self._readEnergyDatGZ('%s.dat.gz'%outputname)
    else:
      energies = self._execute(outputname, 0.0,
        integrator_script, output_script, execution_script,  
        energyFields=energyFields,
        writeEnergyDatGZ=writeEnergyDatGZ,
        keepScript=keepScript, retry=False)
    return energies

  def energies_LJ_ELE_INT(self, outputname, dcdname=None, keepScript=False):
    """
    Calculates Lennard-Jones and electrostatic interaction energies with a grid, and ligand internal energy
    """
    
    if (self.grid==None):
      raise Exception('Function requires grid')
    
    if dcdname==None:
      dcdname = outputname+'.dcd'
    
    (integrator_script,output_script,execution_script) = self._energy_scripts(dcdname)

    grid_script_LJ = self.grid.script_LJ()
    grid_script_ELE = self.grid.script_ELE()

    if (os.path.exists('%s.LJ.dat.gz'%outputname)):
      energies = self._readEnergyDatGZ('%s.LJ.dat.gz'%outputname)
    else:
      energies = self._execute(outputname+'.LJ', 0.0,
        integrator_script, output_script, execution_script, grid_script_LJ,
        energyFields=[8,12], # MISC and POTENTIAL energy fields
        writeEnergyDatGZ=True, keepScript=keepScript, retry=False)
    E_LJ = [energy[0] for energy in energies]
    E_INT = [energy[1]-energy[0] for energy in energies]

    if (os.path.exists('%s.ELE.dat.gz'%outputname)):
      energies = self._readEnergyDatGZ('%s.ELE.dat.gz'%outputname)
    else:
      energies = self._execute(outputname+'.ELE', 0.0,
        integrator_script, output_script, execution_script, grid_script_ELE,
        energyFields=[8], # MISC field has grid energy
        writeEnergyDatGZ=True, keepScript=keepScript, retry=False)
    E_ELE = [energy[0] for energy in energies]
    
    return [E_LJ,E_ELE,E_INT]
