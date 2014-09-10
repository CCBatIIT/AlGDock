from AlGDock.HREX import *

class BindingPMF_plots(BindingPMF):
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

  def show_samples(self, process='dock', state=-1,
      show_original_ligand=True, quit=False):
    if state==-1:
      state = len(self.confs[process]['samples'])-1
    # Gather and write ligand coordinates
    confs = []
    for c in range(len(self.confs[process]['samples'][state])):
      if len(self.confs[process]['samples'][state][c])>0:
        if not isinstance(self.confs[process]['samples'][state][c],list):
          self.confs[process]['samples'][state][c] = \
            [self.confs[process]['samples'][state][c]]
        confs += self.confs[process]['samples'][state][c]
    ligand_dcd_FN = '%s-%05d.dcd'%(process,state)
    self._write_dcd(ligand_dcd_FN, confs, includeLigand=True)
    
    script  =  'set ligand [mol new '+self._FNs['prmtop']['L']+']\n'
    script += 'mol addfile '+ligand_dcd_FN+' type dcd waitfor all\n'
    script += 'mol drawframes $ligand 0 {0:%d}\n'%len(confs)

    # Show the original ligand position
    if show_original_ligand:
      original_ligand_dcd_FN = 'L.dcd'
      self._write_dcd(original_ligand_dcd_FN, self.confs['ligand'], \
        includeLigand=True, includeReceptor=False)
      script += 'set original_ligand [mol new '+self._FNs['prmtop']['L']+']\n'
      script += 'mol addfile '+original_ligand_dcd_FN+' type dcd waitfor all\n'
      script += 'mol modstyle 0 $original_ligand ' + \
                'Licorice 0.300000 10.000000 10.000000\n'
    
    # For docking, write receptor coordinates
    if process=='dock':
      receptor_dcd_FN = 'R.dcd'
      self._write_dcd(receptor_dcd_FN, self.confs['receptor'], \
        includeLigand=False, includeReceptor=True)
      script += 'set receptor [mol new '+self._FNs['prmtop']['R']+']\n'
      script += 'mol addfile '+receptor_dcd_FN+' type dcd waitfor all\n'
      script += 'mol modstyle 0 $receptor NewCartoon 0.300000 10.000000 4.100000 0\n'
      script += 'mol modmaterial 0 $receptor Transparent\n'
    
    center_matrix = '[[[1 0 0 {0[0]}] [0 1 0 {0[1]}] [0 1 0 {0[2]}] [0 0 0 1]]]'
    center_matrix = center_matrix.format(-10*self._forceFields['site'].center)
    center_matrix = center_matrix.replace('[','{').replace(']','}')
    
    script += 'molinfo 0 set center_matrix '+center_matrix+'\n'
    script += 'molinfo 1 set center_matrix '+center_matrix+'\n'
    
    script += 'scale to 0.07\n'
    
#    script += 'rotate y by -25.000000\n'
#    script += 'scale by 2.0\n'
    script += 'axes location Off\n'
    script += 'render snapshot %s-%05d.tga\n'%(process,state)
    if quit:
      script += 'quit\n'
    
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

  def show_all_samples(self, process='dock', stride=1):
    for state_ind in range(0,len(self.confs[process]['samples']),stride):
      self.show_samples(process, state=state_ind, quit=True)
    os.system('convert -delay 10 dock*.tga dock-movie.gif')

