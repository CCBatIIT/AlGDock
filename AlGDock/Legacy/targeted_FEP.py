"""Legacy code

Ideas that were tested and appeared unsuccessful.
The code is not being maintained but is kept in case they are resurrected.
"""

  def targeted_FEP(self):
    start_string = "\n>>> Targeted FEP calculations, starting at " + \
      time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()) + "\n"
    self.log.start_times['targeted_FEP'] = time.time()

    # Cluster poses
    if self.args.params['CD']['pose'] == -1:
      (self.stats_RL['pose_inds'], self.stats_RL['scores']) = \
        self._get_pose_prediction()
    confs = [self.data['CD'].confs['samples'][-1][cycle][n] \
      for (cycle,n) in self.stats_RL['pose_inds']]

    # Set up OpenMM
    process = 'CD'
    moiety = 'RL'
    phase = 'OpenMM_OBC2'
    key = moiety + phase
    self._setup_OpenMM(moiety, phase)

    # Minimize each pose
    import simtk.unit

    receptor_0 = self.data['CD'].confs['receptor'][:self._ligand_first_atom, :]
    receptor_1 = self.data['CD'].confs['receptor'][self._ligand_first_atom:, :]

    confs_min = []
    for conf in confs:
      full_conf = np.vstack(
          (receptor_0, conf[self.molecule.prmtop_atom_order, :], receptor_1))
      self._OpenMM_sims[key].context.setPositions(full_conf)
      s = self._OpenMM_sims[key].context.getState(getEnergy=True)
      Eo = s.getPotentialEnergy() / simtk.unit.kilojoule * simtk.unit.mole
      self._OpenMM_sims[key].minimizeEnergy()
      s = self._OpenMM_sims[key].context.getState(getPositions=True,
                                                  getEnergy=True)
      confs_min.append(s.getPositions(asNumpy=True) \
                       [self._ligand_first_atom:self._ligand_first_atom+self._ligand_natoms,:] \
                       [self.molecule.inv_prmtop_atom_order,:]/simtk.unit.nanometer)

    # Define mappings to minimized poses
    import AlGDock.BAT
    BAT_converter = AlGDock.BAT.converter(self.universe, self.molecule)
    self.stats_RL['mappings'] = []
    for (conf_o, conf_min) in zip(confs, confs_min):
      BAT_o = BAT_converter.BAT(conf_o, extended=True)
      BAT_min = BAT_converter.BAT(conf_min, extended=True)
      mapping = BAT_min - BAT_o
      # Don't change bond length or bond angles
      mapping[6::3] = 0.
      mapping[7::3] = 0.
      self.stats_RL['mappings'].append(mapping)

    # If relevant, store the rmsd of the mapped cluster poses
    if self.args.params['CD']['rmsd']:
      self.stats_RL['mapped_rmsd'] = []
      for mapping in self.stats_RL['mappings']:
        # RMSD of mapped configurations
        confs_map = []
        for conf_o in confs:
          BAT_o = BAT_converter.BAT(conf_o, extended=True)
          confs_map.append(BAT_converter.Cartesian(BAT_o + mapping))
        self.stats_RL['mapped_rmsd'].append(self.get_rmsds(confs_map))

    # Gather all snapshots in milestone D
    for k in range(self.stats_RL['equilibrated_cycle'][-1], \
          self.data[process].cycle):
      if not isinstance(self.data[process].confs['samples'][-1][k], list):
        self.data[process].confs['samples'][-1][k] = [
            self.data[process].confs['samples'][-1][k]
        ]
    import itertools
    confs = np.array([conf for conf in itertools.chain.from_iterable(\
      [self.data[process].confs['samples'][-1][c] \
        for c in range(self.stats_RL['equilibrated_cycle'][-1], self.data[process].cycle)])])

    extractCycles = range(self.stats_RL['equilibrated_cycle'][-1],
                          self.data[process].cycle)
    u_sampled = np.concatenate([\
              self.stats_RL['u_K_sampled'][c] for c in extractCycles])
    f_R_solv = self.args.original_Es[0][0]['R' + phase][:, -1] / self.RT_TARGET

    # For each mapping, map all samples and evalulate energies
    self.stats_RL['u_K_mapped_OpenMM_OBC2'] = []
    for mapping in self.stats_RL['mappings']:
      # Energy of mapped configurations
      umap = []
      for conf_o in confs:
        BAT_o = BAT_converter.BAT(conf_o, extended=True)
        conf_map = BAT_converter.Cartesian(BAT_o + mapping)
        full_conf = np.vstack((receptor_0, \
          conf_map[self.molecule.prmtop_atom_order,:], receptor_1))
        self._OpenMM_sims[key].context.setPositions(full_conf)
        s = self._OpenMM_sims[key].context.getState(getEnergy=True)
        umap.append(s.getPotentialEnergy() / simtk.unit.kilojoule *
                    simtk.unit.mole / self.RT_TARGET)
      umap = np.array(umap)
      self.stats_RL['u_K_mapped_OpenMM_OBC2'].append(umap)

    # For each mapping, calculate the free energy
    self.f_RL['mapped_OpenMM_OBC2_solv'] = []
    for umap in self.stats_RL['u_K_mapped_OpenMM_OBC2']:
      # Targeted FEP
      du = umap - u_sampled
      min_du = min(du)
      weights = np.exp(-du + min_du)

      # Filter outliers
      if self.args.params['CD']['pose'] > -1:
        toKeep = du > (np.mean(du) - 3 * np.std(du))
        du = du[toKeep]
        weights[~toKeep] = 0.

      weights = weights / sum(weights)

      # Exponential average
      f_RL_solv = -np.log(np.exp(-du + min_du).mean()) + min_du - f_R_solv
      self.f_RL['mapped_OpenMM_OBC2_solv'].append(f_RL_solv)

    # Select mapping that leads to the lowest free energy
    selected_mapping = np.argmin(self.f_RL['mapped_OpenMM_OBC2_solv'])

    # Calculate BPMF based on targeted FEP with the selected mapping
    self.B['mapped_OpenMM_OBC2_MBAR'] = self.B['OpenMM_OBC2_MBAR'][-1] \
      - self.f_RL['OpenMM_OBC2_solv'][-1] \
      + self.f_RL['mapped_OpenMM_OBC2_solv'][selected_mapping]

    # Score poses based on selected mapping
    if self.args.params['CD']['rmsd']:
      self.stats_RL['scores']['mapped_rmsd'] = self.stats_RL['mapped_rmsd'][
          selected_mapping]

    un = self.stats_RL['u_K_mapped_OpenMM_OBC2'][selected_mapping]
    du = un - u_sampled
    min_du = min(du)
    weights = np.exp(-du + min_du)

    rmsd_matrix = self._get_rmsd_matrix()
    assignments = self._cluster_samples(rmsd_matrix)
    cluster_counts = np.histogram(assignments, \
      bins=np.arange(len(set(assignments))+1)-0.5,
      weights=weights)[0]

    cluster_fe = -np.log(cluster_counts)
    cluster_fe -= np.min(cluster_fe)
    self.stats_RL['scores']['mapped_' + phase + '_fe_u'] = cluster_fe
    # by minimum and mean energy
    self.stats_RL['scores']['mapped_' + phase + '_min_u'] = []
    self.stats_RL['scores']['mapped_' + phase + '_mean_u'] = []
    for n in range(max(assignments) + 1):
      un_n = [un[i] for i in range(len(assignments)) if assignments[i] == n]
      self.stats_RL['scores']['mapped_' + phase + '_min_u'].append(
          np.min(un_n))
      self.stats_RL['scores']['mapped_' + phase + '_mean_u'].append(
          np.mean(un_n))

    # Store data
    if self.args.params['CD']['pose'] == -1:
      f_RL_FN = os.path.join(self.args.dir['CD'], 'f_RL.pkl.gz')
    else:
      f_RL_FN = os.path.join(self.args.dir['CD'], \
        'f_RL_pose%03d.pkl.gz'%self.args.params['CD']['pose'])
    self.log.tee(write_pkl_gz(f_RL_FN, (self.f_L, self.stats_RL, self.f_RL, self.B)))
    self.log.tee("\nElapsed time for targeted FEP: " + \
      HMStime(time.time()-self.log.start_times['targeted_FEP']))
