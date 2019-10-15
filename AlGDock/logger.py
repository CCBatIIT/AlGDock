import os
import time


class NullDevice():
  """A device to suppress output
  """
  def write(self, s):
    pass

  def flush(self):
    pass


class Logger:
  """Handles log file, locks, and timers

  ...

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
  """
  def __init__(self, args, max_time=None, run_type=None):
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

  def recordStart(self, event_key):
    """Records the time at the beginning of an event

    Parameters
    ----------
    event_key : str
      Event description
    """
    self.start_times[event_key] = time.time()

  def timeSince(self, event_key):
    """Amount of time since an event

    Parameters
    ----------
    event_key : str
      Event description

    Returns
    -------
    float
      The amount of time since the event, in seconds
    """
    if event_key in self.start_times.keys():
      return time.time() - self.start_times[event_key]
    else:
      raise Exception(event_key + ' has not started!')

  def isTimeRemaining(self):
    """Checks whether there is any time remaining

    Returns
    -------
    bool
      If True, there is time remaining
    """
    if self.run_type.startswith('timed'):
      time_since_start = (time.time() - self.start_times['run'])
      remaining_time = self.max_time * 60 - time_since_start
      if remaining_time < 0:
        return False
    else:
      return True

  def isTimeForTask(self, task_times):
    """Projects whether there is enough time remaining to complete a task

    Parameters
    ----------
    task_times : list of float
      The amount of time required for similar previous tasks

    Returns
    -------
    bool
      If True, there is expected to be sufficient time remaining to complete the task
    """
    if self.run_type.startswith('timed'):
      time_since_start = (time.time() - self.start_times['run'])
      remaining_time = self.max_time * 60 - time_since_start
      mean_task_time = np.mean(task_times)
      self.tee("  projected task time: %s, remaining time: %s"%(\
        HMStime(mean_task_time), HMStime(remaining_time)), process=process)
      if mean_task_time > remaining_time:
        return False
    else:
      return True
