import os

from AlGDock.IO import load_pkl_gz
from AlGDock.IO import write_pkl_gz


class SimulationData:
  """Stores simulation data

  ...

  Attributes
  ----------
  confs : dict of lists of arrays
    Configurations
  cycle : int
    current simulation cycle
  Es : dict
    Energies of sampled configurations
  pose : int
    pose number
  process : str
    Process, either 'BC' or 'CD'
  protocol : list of dict
    Description of alchemical states
  """
  def __init__(self, dir, process, pose):
    """Initializes the class

    Parameters
    ----------
    dir : str
      The directory to store the data
    process : str
      Process, either 'BC' or 'CD'
    pose : int
      pose number
    """
    self.dir = dir
    self.process = process
    self.pose = pose
    self._load_pkl_gz()

  def _load_pkl_gz(self):
    """Load sampled configurations and energies

    """
    if self.process == 'CD' and self.pose > -1:
      progress_FN = os.path.join(
        self.dir,
        '%s_progress_self.pose%03d.pkl.gz' % (self.process, self.pose))
      data_FN = os.path.join(
        self.dir, '%s_data_self.pose%03d.pkl.gz' % (self.process, self.pose))
    else:
      progress_FN = os.path.join(self.dir,
                                 '%s_progress.pkl.gz' % (self.process))
      data_FN = os.path.join(self.dir, '%s_data.pkl.gz' % (self.process))

    saved = {
      'progress': load_pkl_gz(progress_FN),
      'data': load_pkl_gz(data_FN)
    }
    if (saved['progress'] is None) or (saved['data'] is None):
      if os.path.isfile(progress_FN):
        os.remove(progress_FN)
      if os.path.isfile(data_FN):
        os.remove(data_FN)
      if self.process == 'CD' and self.pose > -1:
        progress_FN = os.path.join(
          self.dir,
          '%s_progress_self.pose%03d.pkl.gz.BAK' % (self.process, self.pose))
        data_FN = os.path.join(
          self.dir,
          '%s_data_self.pose%03d.pkl.gz.BAK' % (self.process, self.pose))
      else:
        progress_FN = os.path.join(self.dir,
                                   '%s_progress.pkl.gz.BAK' % (self.process))
        data_FN = os.path.join(self.dir, '%s_data.pkl.gz.BAK' % (self.process))

      saved = {
        'progress': load_pkl_gz(progress_FN),
        'data': load_pkl_gz(data_FN)
      }
      if (saved['progress'] is None):
        print '  no progress information for %s' % self.process
      elif (saved['data'] is None):
        saved['progress'] = None
        print '  missing data in %s' % self.process
      else:
        print '  using stored progress and data in %s' % self.process
    self.clear()

    if saved['progress'] is not None:
      self.protocol = saved['progress'][1]
      self.cycle = saved['progress'][2]
    if saved['data'] is not None:
      if self.process == 'CD' and saved['data'][0] is not None:
        (self._n_trans, self._max_n_trans, self._random_trans, \
         self._n_rot, self._max_n_rot, self._random_rotT) = saved['data'][0]
      # New file format (after 6/13/2016) storing starting self.poses
      self.confs['starting_self.poses'] = saved['data'][1]
      self.confs['replicas'] = saved['data'][2]
      self.confs['seeds'] = saved['data'][3]
      self.confs['SmartDarting'] = saved['data'][4]
      self.confs['samples'] = saved['data'][5]
      self.Es = saved['data'][6]
      if saved['data'][5] is not None:
        self.cycle = len(saved['data'][5][-1])
      else:
        self.cycle = 0
    if self.protocol==[] or \
        (not self.protocol[-1]['crossed']):
      self.cycle = 0

  def _save_pkl_gz(self):
    """Saves sampled configurations and energies

    """
    random_orient = None
    if self.process == 'CD' and hasattr(self, '_n_trans'):
      random_orient = (self._n_trans, self._max_n_trans, self._random_trans, \
         self._n_rot, self._max_n_rot, self._max_n_rot)

    saved = {
      'data': (random_orient, self.confs['starting_self.poses'],
               self.confs['replicas'], self.confs['seeds'],
               self.confs['SmartDarting'], self.confs['samples'], self.Es)
    }

    key = 'data'
    if self.process == 'CD' and self.pose > -1:
      saved_FN = os.path.join(self.dir,'%s_%s_self.pose%03d.pkl.gz'%(\
        self.process, key, self.pose))
    else:
      saved_FN = os.path.join(self.dir, '%s_%s.pkl.gz' % (self.process, key))
    if not os.path.isdir(self.dir):
      os.system('mkself.dir -p ' + self.dir)
    if os.path.isfile(saved_FN):
      os.rename(saved_FN, saved_FN + '.BAK')
    return (write_pkl_gz(saved_FN, saved[key]) +
            '\n  saved %s data' % self.process)

  def clear(self):
    """Sets up empty variables for protocol, cycles, and configurations.

    """
    self.protocol = []
    self.cycle = 0
    self.confs = {}
    self.confs['starting_self.poses'] = None
    self.confs['replicas'] = None
    self.confs['seeds'] = None
    self.confs['SmartDarting'] = []
    self.confs['samples'] = None
    self.Es = None
