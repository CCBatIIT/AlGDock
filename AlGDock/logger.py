import os

class Logger:
  """Handles log file, locks, and timers

  Attributes
  ----------
  args : argument_parser.SimulationArguments object
    Arguments for the simulation
  log : file
    Log file
  max_time : float
    Maximum simulation time, in minutes
  start_times : dict of float
    Times when calculations started
  timings : dict of float
    Amount of time spend in each type of calculation


  Methods
  -------
  set_lock
  clear_lock
  tee
  """
  def __init__(self, args, max_time = None, run_type = None):
    """

    Parameters
    ----------
    args : argument_parser.SimulationArguments object
      Arguments for the simulation
    max_time : float
      Maximum simulation time, in minutes
    run_type : str
      Type of simulation
    """
    self.args = args
    self.max_time = max_time
    self.run_type = run_type
    self.start_times = {}
    self.timings = {}

  def set_lock(self, process):
    """Creates a lock and opens a log file

    Parameters
    ----------
    process : str
      Process, either 'BC' or 'CD'
    """
    if not os.path.isdir(self.args.dir[process]):
      os.system('mkdir -p ' + self.args.dir[process])
    if process == 'CD' and self.args.params['CD']['pose'] > -1:
      lockFN = os.path.join(self.args.dir[process], \
        '.lock_pose%03d'%self.args.params['CD']['pose'])
    else:
      lockFN = os.path.join(self.args.dir[process], '.lock')
    if os.path.isfile(lockFN):
      raise Exception(process + ' is locked')
    else:
      lockF = open(lockFN, 'w')
      lockF.close()
    if process == 'CD' and self.args.params['CD']['pose'] > -1:
      logFN = os.path.join(self.args.dir[process],'%s_pose%03d_log.txt'%(\
        process, self.args.params['CD']['pose']))
    else:
      logFN = os.path.join(self.args.dir[process], process + '_log.txt')
    self.log = open(logFN, 'a')

  def clear_lock(self, process):
    """Clears a lock and closes a log file

    Parameters
    ----------
    process : str
      Process, either 'BC' or 'CD'
    """
    if process == 'CD' and self.args.params['CD']['pose'] > -1:
      lockFN = os.path.join(self.args.dir[process], \
        '.lock_pose%03d'%self.args.params['CD']['pose'])
    else:
      lockFN = os.path.join(self.args.dir[process], '.lock')
    if os.path.isfile(lockFN):
      os.remove(lockFN)
    if hasattr(self, 'log'):
      self.log.close()
      del self.log

  def tee(self, var, process=None):
    """Sends output to both the standard output and the log file, if applicable

    Parameters
    ----------
    var : str
      The object to write
    process : str
      Process, either 'BC' or 'CD' or None
    """
    print var
    if hasattr(self, 'log'):
      if isinstance(var, str):
        self.log.write(var + '\n')
      else:
        self.log.write(repr(var) + '\n')
      self.log.flush()
    elif process is not None:
      self.set_lock(process)
      if isinstance(var, str):
        self.log.write(var + '\n')
      else:
        self.log.write(repr(var) + '\n')
      self.log.flush()
      self.clear_lock(process)
