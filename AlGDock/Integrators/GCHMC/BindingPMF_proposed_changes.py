# export LD_LIBRARY_PATH=/home/lspirido/AstexDiv_xtal/scripts/GCHMC/:$LD_LIBRARY_PATH
  def _mixed_sampler(self, steps, steps_per_trial, T, delta_t, random_seed, normalize=False, adapt=False):
    ntrials = steps/steps_per_trial
    if(ntrials > 5):
      Csteps = steps/5
      Csteps_per_trial = steps_per_trial / 5

      (mixed_confs0, mixed_potEs0, mixed_accs0, mixed_ntrials0, mixed_dts0) = \
        self.mixed_samplers[0].Call(Csteps, Csteps_per_trial, T, delta_t, (random_seed))

      (mixed_confs1, mixed_potEs1, mixed_accs1, mixed_ntrials1, mixed_dts1) = \
        self.mixed_samplers[1](steps=Csteps, steps_per_trial=Csteps_per_trial, T=T, delta_t=delta_t, \
        normalize=normalize, adapt=adapt, random_seed=random_seed)

      (mixed_confs2, mixed_potEs2, mixed_accs2, mixed_ntrials2, mixed_dts2) = \
        self.mixed_samplers[0].Call(Csteps, Csteps_per_trial, T, delta_t, (random_seed))

      (mixed_confs3, mixed_potEs3, mixed_accs3, mixed_ntrials3, mixed_dts3) = \
        self.mixed_samplers[1](steps=Csteps, steps_per_trial=Csteps_per_trial, T=T, delta_t=delta_t, \
        normalize=normalize, adapt=adapt, random_seed=random_seed)

      (mixed_confs4, mixed_potEs4, mixed_accs4, mixed_ntrials4, mixed_dts4) = \
        self.mixed_samplers[0].Call(Csteps, Csteps_per_trial, T, delta_t, (random_seed))
    else:
      (mixed_confs3, mixed_potEs3, mixed_accs3, mixed_ntrials3, mixed_dts3) = \
        self.mixed_samplers[1](steps=steps, steps_per_trial=steps_per_trial, T=T, delta_t=delta_t, \
        normalize=normalize, adapt=adapt, random_seed=random_seed)

      (mixed_confs4, mixed_potEs4, mixed_accs4, mixed_ntrials4, mixed_dts4) = \
        self.mixed_samplers[0].Call(steps, steps_per_trial, T, delta_t, (random_seed))
      return (mixed_confs4, mixed_potEs4, mixed_accs4, mixed_ntrials4, mixed_dts4)

    mixed_confs = mixed_confs0 + mixed_confs1 + mixed_confs2 + mixed_confs3 + mixed_confs4
    mixed_potEs = mixed_potEs0 + mixed_potEs1 + mixed_potEs2 + mixed_potEs3 + mixed_potEs4
    mixed_accs = (mixed_accs0 + mixed_accs1 + mixed_accs2 + mixed_accs3 + mixed_accs4) / 5
    return (mixed_confs, mixed_potEs, mixed_accs, ntrials, delta_t)
  #
## [...]

    for p in ['cool', 'dock']:
      if self.params[p]['sampler'] == 'TDHMC':
        self.mixed_samplers = []
        from GCHMC.GCHMC import GCHMCIntegrator
        GCHMCintegrator = GCHMCIntegrator(self.universe, os.path.dirname(self._FNs['ligand_database']), os.path.dirname(self._FNs['forcefield']))
        self.mixed_samplers.append(GCHMCintegrator)
        from AlGDock.Integrators.HamiltonianMonteCarlo.HamiltonianMonteCarlo import HamiltonianMonteCarloIntegrator
        self.mixed_samplers.append(HamiltonianMonteCarloIntegrator(self.universe))
        self.sampler[p] = self._mixed_sampler
## [...]
      for T in T_SERIES:
        self.sampler['cool'](steps = 500, steps_per_trial = 100, T=T,\
                             delta_t=self.delta_t, random_seed=(random_seed))
## [...]
          for T in T_SERIES:
            self.sampler['dock'](steps = 500, steps_per_trial = 100, T=T,\
                                 delta_t=self.delta_t, random_seed=random_seed)
## [...]

  def _sim_one_state(self, seed, process, lambda_k, \
      initialize=False, reference=0):
    # [...]
    dat = sampler(\
      steps=steps, steps_per_trial=steps_per_trial, \
      T=lambda_k['T'], delta_t=delta_t, \
      normalize=(process=='cool'), adapt=initialize, random_seed=random_seed)
    results['acc_Sampler'] = dat[2]
    results['att_Sampler'] = dat[3]
    results['delta_t'] = dat[4]
    results['time_Sampler'] = (time.time() - time_start_Sampler)

##

