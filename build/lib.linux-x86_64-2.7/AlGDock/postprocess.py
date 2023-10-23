# In APBS, minimum ratio of PB grid length to maximum dimension of solute
LFILLRATIO = 4.0  # For the ligand
RFILLRATIO = 2.0  # For the receptor/complex


# TODO: Pass more necessary parameters
# TODO: Transfer code from BindingPMF.py once other encapsulation is complete
class Postprocess:
  """Postprocesses sampled configurations to obtain energies in target force fields
  Methods
  -------
  __init__
  _load_programs
  _postprocess
  _energy_worker
  _NAMD_Energy
  _sander_Energy
  _get_elsize
  _gbnsr6_Energy
  _setup_OpenMM
  _OpenMM_Energy
  _APBS_Energy
  roundUpDime
  _get_APBS_grid_spacing
  roundUpDime
  _combine_MM_and_solvent
  _write_traj
  """
  def __init__(self, args, data, log):
    self.args = args
    self.data = data
    self.log = log

    # Locate programs for postprocessing
    all_phases = self.args.params['CD']['phases'] + self.args.params['BC'][
      'phases']
    self._load_programs(all_phases)

    # Determine APBS grid spacing
    if 'APBS_PBSA' in self.args.params['CD']['phases'] or \
       'APBS_PBSA' in self.args.params['BC']['phases']:
      self._get_APBS_grid_spacing()

    # Determines receptor electrostatic size
    if np.array([p.find('ALPB') > -1 for p in all_phases]).any():
      self._get_elsize()
