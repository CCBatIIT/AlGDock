#!/usr/bin/env python

from AlGDock.BindingPMF import *

class BPMF_plots(BPMF):
  def plot_energies(self, process='cool', firstCycle=0, toCycle=None):
    """
    Plots timeseries and histograms of the energies for each state.
    Requires matplotlib extension for python.
    """
    try:
      import matplotlib.pyplot as plt  # @UnresolvedImport
    except:
      print 'plot_energies requires matplotlib'
      return
    
    K = len(getattr(self,process+'_protocol'))
    
    if toCycle is None:
      toCycle = getattr(self,'_%s_cycle'%process)
    Es = [getattr(self,process+'_Es')[0][firstCycle:toCycle]]
    for Es_state in getattr(self,process+'_Es'):
      Es.append(Es_state[firstCycle:toCycle])
    (u_kln,N_k) = self._u_kln(Es, getattr(self,process+'_protocol'), noBeta=True)
    
    plt.figure(0)
    for k in range(K-1):
      plt.plot(u_kln[k,k,:N_k[k]],'.-')
    
    plt.figure(1)
    for k in range(K-1):
      (p_k,x_k) = np.histogram(u_kln[k,k,:N_k[k]], bins=20)
      x_c = x_k[:-1]+(x_k[1]-x_k[0])/2.
      plt.plot(x_c,p_k,'.-')

    plt.figure(2)
    for k in range(K-1):
      i1 = k
      i2 = k+1
      min_N_k = min(N_k[i1],N_k[i2])
      e1 = u_kln[i1,i1,:min_N_k]
      e2 = u_kln[i2,i2,:min_N_k]
      Erange = (max(min(e1),min(e2)),min(max(e1),max(e2)))

      (p_h,x_h) = np.histogram(e1, bins=20, range=Erange)
      (p_l,x_l) = np.histogram(e2, bins=20, range=Erange)
      x_c = x_h[:-1]+(x_h[1]-x_h[0])/2.
      plt.plot(x_c,np.log(np.array(p_h,dtype=float)/np.array(p_l,dtype=float)),'.-')

  def plot_energy_ratio(self):
    try:
      import matplotlib.pyplot as plt  # @UnresolvedImport
    except:
      print 'plot_energy_ratio requires matplotlib'
      return

    K = len(self.dock_protocol)
    e_ratio = []
    for k in range(K):
      e_ratio_k = np.array([])
      for c in range(len(self.dock_Es[k])):
        e_ratio_k = np.hstack((e_ratio_k,np.abs(self.dock_Es[k][c]['ELE']/self.dock_Es[k][c]['sLJr'])))
      e_ratio.append(e_ratio_k)
    plt.plot(np.transpose(e_ratio))

  def show_replicas(self, process='dock'):
    import AlGDock.IO
    IO_dcd = AlGDock.IO.dcd(self.molecule,
      ligand_atom_order = self.molecule.prmtop_atom_order, \
      receptorConf = self.confs['receptor'], \
      ligand_first_atom = self._ligand_first_atom)
      
    dcd_FN = 'replicas.dcd'
    confs = self.confs['dock']['replicas']
    IO_dcd.write(dcd_FN, confs, includeLigand=True, includeReceptor=False)
    
    script = ''
    
    # Show the original ligand position
    original_ligand_dcd_FN = 'L.dcd'
    IO_dcd.write(original_ligand_dcd_FN, self.confs['ligand'], \
      includeLigand=True, includeReceptor=False)
    script += 'set original_ligand [mol new '+self._FNs['prmtop']['L']+']\n'
    script += 'mol addfile '+original_ligand_dcd_FN+' type dcd waitfor all\n'
    script += 'mol modstyle 0 $original_ligand ' + \
              'Licorice 0.300000 10.000000 10.000000\n'
    
    # Show receptor coordinates
    receptor_dcd_FN = 'R.dcd'
    IO_dcd.write(receptor_dcd_FN, self.confs['receptor'], \
      includeLigand=False, includeReceptor=True)
    script += 'set receptor [mol new '+self._FNs['prmtop']['R']+']\n'
    script += 'mol addfile '+receptor_dcd_FN+' type dcd waitfor all\n'
    script += 'mol modstyle 0 $receptor NewCartoon 0.300000 10.000000 4.100000 0\n'
    script += 'mol modmaterial 0 $receptor Transparent\n'

    # Show samples
    script +=  'set ligand [mol new '+self._FNs['prmtop']['L']+']\n'
    script += 'mol addfile '+dcd_FN+' type dcd waitfor all\n'

    script_FN = 'show_samples.vmd'
    script_F = open('show_samples.vmd','w')
    script_F.write(script)
    script_F.close()

    import subprocess
    subprocess.call([self._FNs['vmd'], '-e', script_FN, '-size', '800', '800'])
    for FN in [ligand_dcd_FN, original_ligand_dcd_FN, \
               receptor_dcd_FN, 'show_samples.vmd']:
      if os.path.isfile(FN):
        os.remove(FN)

  def show_samples(self, process='dock', state=-1, \
      show_original_ligand=True, show_receptor=False, \
      save_image=False, scale=True, execute=True, quit=False):
    if state==-1:
      state = len(self.confs[process]['samples'])-1
    import AlGDock.IO
    IO_dcd = AlGDock.IO.dcd(self.molecule,
      ligand_atom_order = self.molecule.prmtop_atom_order, \
      receptorConf = self.confs['receptor'], \
      ligand_first_atom = self._ligand_first_atom)
      
    # Gather and write ligand coordinates
    confs = []
    for c in range(len(self.confs[process]['samples'][state])):
      if len(self.confs[process]['samples'][state][c])>0:
        if not isinstance(self.confs[process]['samples'][state][c],list):
          self.confs[process]['samples'][state][c] = \
            [self.confs[process]['samples'][state][c]]
        confs += self.confs[process]['samples'][state][c]
    rep = 0
    while os.path.isfile('%s-%05d-%d.dcd'%(process,state,rep)):
      rep += 1
    ligand_dcd_FN = '%s-%05d-%d.dcd'%(process,state,rep)
    IO_dcd.write(ligand_dcd_FN, confs, includeLigand=True, includeReceptor=False)
    
    script = ''
    
    # Show the original ligand position
    original_ligand_dcd_FN = 'L.dcd'
    if show_original_ligand:
      IO_dcd.write(original_ligand_dcd_FN, self.confs['ligand'], \
        includeLigand=True, includeReceptor=False)
      script += 'set original_ligand [mol new '+self._FNs['prmtop']['L']+']\n'
      script += 'mol addfile '+original_ligand_dcd_FN+' type dcd waitfor all\n'
      script += 'mol modstyle 0 $original_ligand ' + \
                'Licorice 0.300000 10.000000 10.000000\n'
    
    # For docking, write receptor coordinates
    receptor_dcd_FN = 'R.dcd'
    if show_receptor:
      IO_dcd.write(receptor_dcd_FN, self.confs['receptor'], \
        includeLigand=False, includeReceptor=True)
      script += 'set receptor [mol new '+self._FNs['prmtop']['R']+']\n'
      script += 'mol addfile '+receptor_dcd_FN+' type dcd waitfor all\n'
      script += 'mol modstyle 0 $receptor NewCartoon 0.300000 10.000000 4.100000 0\n'
      script += 'mol modmaterial 0 $receptor Transparent\n'

    # Show samples
    script +=  'set ligand [mol new '+self._FNs['prmtop']['L']+']\n'
    script += 'mol addfile '+ligand_dcd_FN+' type dcd waitfor all\n'
    script += 'mol drawframes $ligand 0 {0:%d}\n'%len(confs)

    if scale:
      script += 'scale to 0.07\n'
      script += 'axes location Off\n'

    if save_image:
      script += 'render snapshot %s-%05d.tga\n'%(process,state)
    if quit:
      script += 'quit\n'
    if execute:
      script_FN = 'show_samples.vmd'
      script_F = open('show_samples.vmd','w')
      script_F.write(script)
      script_F.close()

      import subprocess
      subprocess.call([self._FNs['vmd'], '-e', script_FN, '-size', '800', '800'])
      for FN in [ligand_dcd_FN, original_ligand_dcd_FN, \
                 receptor_dcd_FN, 'show_samples.vmd']:
        if os.path.isfile(FN):
          os.remove(FN)

    return script

  def show_all_samples(self, process='dock', stride=1):
    for state_ind in range(0,len(self.confs[process]['samples']),stride):
      self.show_samples(process, state=state_ind, quit=True)
    os.system('convert -delay 10 dock*.tga dock-movie.gif')

