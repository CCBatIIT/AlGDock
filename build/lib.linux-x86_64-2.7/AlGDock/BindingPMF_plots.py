#!/usr/bin/env python

from AlGDock.BindingPMF import *

# vmd procedure to align the principal axes of the ligand to y and x.
vmd_principal_axes_procedure = '''
package require Orient
namespace import Orient::orient

set sel [atomselect 0 "all"]
set I [Orient::calc_principalaxes $sel]
set M1 [orient $sel [lindex $I 2] {0 1 0}]
$sel move $M1
set I [Orient::calc_principalaxes $sel]
set M2 [orient $sel [lindex $I 1] {1 0 0}]
$sel move $M2

for {set molid 1} {$molid < %d} {incr molid} {
  set nframes [molinfo $molid get numframes]
  for {set frameid 0} {$frameid < $nframes} {incr frameid} {
    set sel [atomselect $molid "all" frame $frameid]
    $sel move $M1
    $sel move $M2
  }
}
'''

# If the protein center of mass is positive (in front of the ligand),
# flip it around
vmd_unblock_procedure = '''
proc center_of_mass {selection} {
  # some error checking
  if {[$selection num] < 0} {
          error "center_of_mass: needs a selection with atoms"
  }
  # set the center of mass to 0
  set com [veczero]
  # set the total mass to 0
  set mass 0
  # [$selection get {x y z}] returns the coordinates {x y z} 
  # [$selection get {mass}] returns the masses
  # so the following says "for each pair of {coordinates} and masses,
  #  do the computation ..."
  foreach coord [$selection get {x y z}] m [$selection get mass] {
     # sum of the masses
     set mass [expr $mass + $m]
     # sum up the product of mass and coordinate
     set com [vecadd $com [vecscale $m $coord]]
  }
  # and scale by the inverse of the number of atoms
  if {$mass == 0} {
          error "center_of_mass: total mass is zero"
  }
  # The "1.0" can't be "1", since otherwise integer division is done
  return [vecscale [expr 1.0/$mass] $com]
}

set sel [atomselect %d "all"]
if {[expr [lindex [center_of_mass $sel] 2]>0]} {
  echo Flipped to unblock view
  rotate x by 180
}

'''


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

    K = len(getattr(self, process + '_protocol'))

    if toCycle is None:
      toCycle = getattr(self, '_%s_cycle' % process)
    Es = [getattr(self,process+'_Es')[k][firstCycle:toCycle] \
      for k in range(len(getattr(self,process+'_Es')))]
    (u_kln, N_k) = self._u_kln(Es,
                               getattr(self, process + '_protocol'),
                               noBeta=True)

    plt.figure(0)
    for k in range(K - 1):
      plt.plot(u_kln[k, k, :N_k[k]], '.-')

    plt.figure(1)
    for k in range(K - 1):
      (p_k, x_k) = np.histogram(u_kln[k, k, :N_k[k]], bins=20)
      x_c = x_k[:-1] + (x_k[1] - x_k[0]) / 2.
      plt.plot(x_c, p_k, '.-')

    plt.figure(2)
    for k in range(K - 1):
      i1 = k
      i2 = k + 1
      min_N_k = min(N_k[i1], N_k[i2])
      e1 = u_kln[i1, i1, :min_N_k]
      e2 = u_kln[i2, i2, :min_N_k]
      Erange = (max(min(e1), min(e2)), min(max(e1), max(e2)))

      (p_h, x_h) = np.histogram(e1, bins=20, range=Erange)
      (p_l, x_l) = np.histogram(e2, bins=20, range=Erange)
      x_c = x_h[:-1] + (x_h[1] - x_h[0]) / 2.
      plt.plot(x_c,
               np.log(np.array(p_h, dtype=float) / np.array(p_l, dtype=float)),
               '.-')

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
        e_ratio_k = np.hstack(
          (e_ratio_k,
           np.abs(self.dock_Es[k][c]['ELE'] / self.dock_Es[k][c]['sLJr'])))
      e_ratio.append(e_ratio_k)
    plt.plot(np.transpose(e_ratio))

  def show_replicas(self, process='dock', \
        show_ref_ligand=True, show_starting_pose=True, show_receptor=False, \
        save_image=False, image_labels=None, execute=True, \
        principal_axes_alignment=False, clear_files=True, quit=False, \
        view_args={}):
    """
    Show replicas from replica exchange
    """
    confs = self.confs['dock']['replicas']
    prefix = 'replicas'
    return self.show_poses(confs, prefix, \
      thickness=None,
      show_ref_ligand=show_ref_ligand, \
      show_starting_pose=show_starting_pose,
      show_receptor=show_receptor, \
      save_image=save_image, \
      image_labels=image_labels,
      execute=execute, \
      principal_axes_alignment=principal_axes_alignment, \
      clear_files=clear_files, \
      quit=quit, \
      view_args=view_args)

  def show_samples(self, prefix=None, process='dock', state=-1, \
        show_ref_ligand=True, show_starting_pose=True, show_receptor=False, \
        save_image=False, image_labels=None, execute=True, \
        principal_axes_alignment=False, clear_files=True, quit=False, \
        view_args={}):
    """
    Show samples from replica exchange
    """
    if state == -1:
      state = len(self.confs[process]['samples']) - 1
      if prefix is None:
        prefix = '%s-last' % (process)
    else:
      if prefix is None:
        prefix = '%s-%05d' % (process, state)

    if process == 'dock':
      first_cycle = self.stats_RL['equilibrated_cycle'][-1]
    else:
      first_cycle = self.stats_L['equilibrated_cycle'][-1]

    confs = []
    for c in range(first_cycle, len(self.confs[process]['samples'][state])):
      if len(self.confs[process]['samples'][state][c]) > 0:
        if not isinstance(self.confs[process]['samples'][state][c], list):
          self.confs[process]['samples'][state][c] = \
            [self.confs[process]['samples'][state][c]]
        confs += self.confs[process]['samples'][state][c]
    return self.show_poses(confs, prefix, \
      thickness=None,
      show_ref_ligand=show_ref_ligand, \
      show_starting_pose=show_starting_pose,
      show_receptor=show_receptor, \
      save_image=save_image, \
      image_labels=image_labels,
      execute=execute, \
      principal_axes_alignment=principal_axes_alignment, \
      clear_files=clear_files, \
      quit=quit, \
      view_args=view_args)

  def show_pose_prediction(self, score='OpenMM_OBC2_fe_u', \
        show_ref_ligand=True, show_starting_pose=False, show_receptor=True, \
        save_image=False, image_labels=None, execute=True, \
        principal_axes_alignment=False, clear_files=True, quit=False, \
        view_args={}):
    ws = np.exp(-self.stats_RL['scores'][score] / self.RT_TARGET)
    ws = ws / sum(ws)
    toShow = np.arange(len(ws))[ws > 0.001]

    confs = [self.confs['dock']['samples'][-1][cycle][n] \
      for (cycle,n) in self.stats_RL['pose_inds']]
    confs = [confs[n] for n in np.arange(len(ws))[toShow]]
    return self.show_poses(confs, 'prediction-'+score,
      thickness=ws[toShow], \
      show_ref_ligand=show_ref_ligand, \
      show_starting_pose=show_starting_pose,
      show_receptor=show_receptor, \
      save_image=save_image, \
      image_labels=image_labels,
      execute=execute, \
      principal_axes_alignment=principal_axes_alignment, \
      clear_files=clear_files, \
      quit=quit, \
      view_args=view_args)

  def show_poses(self, confs, prefix, \
        thickness=None,
        show_ref_ligand=True, show_starting_pose=True, show_receptor=False, \
        save_image=False, image_labels=None, execute=True, \
        principal_axes_alignment=False, clear_files=True, quit=False, \
        view_args={}):
    """
    Show poses in the context of the protein
    """
    # Write ligand coordinates
    import AlGDock.IO
    IO_dcd = AlGDock.IO.dcd(self.molecule,
      ligand_atom_order = self.molecule.prmtop_atom_order, \
      receptorConf = self.confs['receptor'], \
      ligand_first_atom = self._ligand_first_atom)

    rep = 0
    import os
    while os.path.isfile('%s-%d.dcd' % (prefix, rep)):
      rep += 1
    ligand_dcd_FN = os.path.join(self.dir[process],
                                 '%s-%d.dcd' % (prefix, rep))
    IO_dcd.write(ligand_dcd_FN,
                 confs,
                 includeLigand=True,
                 includeReceptor=False)

    script = ''
    molids = []

    # Show the reference ligand position
    ref_ligand_dcd_FN = os.path.join(self.dir[process], 'L.dcd')
    if show_ref_ligand:
      IO_dcd.write(ref_ligand_dcd_FN, self.confs['ligand'], \
        includeLigand=True, includeReceptor=False)
      script += 'set ref_ligand [mol new ' + self._FNs['prmtop']['L'] + ']\n'
      script += 'mol addfile ' + ref_ligand_dcd_FN + ' type dcd waitfor all\n'
      script += 'mol modstyle 0 $ref_ligand ' + \
                'Licorice 0.250000 10.000000 10.000000\n'
      script += 'mol rename $ref_ligand {Reference Pose}\n'
      # Change atom type to change coloring
      script += 'set sel [atomselect $ref_ligand carbon]\n'
      script += '$sel set type J\n'
      script += 'color Type J purple\n'
      script += 'mol modcolor 0 $ref_ligand Type\n'
      molids.append('Reference Pose')

    # Show the starting pose
    start_ligand_dcd_FN = os.path.join(self.dir[process], 'L_start.dcd')
    if show_starting_pose:
      if self.confs['dock']['starting_poses'] is not None:
        IO_dcd.write(start_ligand_dcd_FN, \
          self.confs['dock']['starting_poses'][-1], \
          includeLigand=True, includeReceptor=False)
        script += 'set start_ligand [mol new ' + self._FNs['prmtop'][
          'L'] + ']\n'
        script += 'mol addfile ' + start_ligand_dcd_FN + ' type dcd waitfor all\n'
        script += 'mol modstyle 0 $start_ligand ' + \
                  'Licorice 0.250000 10.000000 10.000000\n'
        # Change atom type to change coloring
        script += 'set sel [atomselect $start_ligand carbon]\n'
        script += '$sel set type R\n'
        script += 'color Type R green2\n'
        script += 'mol modcolor 0 $start_ligand Type\n'  # green
        script += 'mol rename $start_ligand {Starting Pose}\n'
        molids.append('Starting Pose')

    # For docking, write complex coordinates
    complex_dcd_FN = os.path.join(self.dir[process], 'RL.dcd')
    if show_receptor:
      IO_dcd.write(complex_dcd_FN, self.confs['ligand'], \
        includeLigand=True, includeReceptor=True)
      script += 'set complex [mol new ' + self._FNs['prmtop']['RL'] + ']\n'
      script += 'mol addfile ' + complex_dcd_FN + ' type dcd waitfor all\n'
      script += 'mol modstyle 0 $complex Ribbons 0.200000 25.000000 2.000000\n'
      script += 'mol modcolor 0 $complex ColorID 9\n'  # Pink backbone
      script += 'mol modmaterial 0 $complex MetallicPastel\n'
      script += 'mol rename $complex {Complex}\n'
      molids.append('Receptor')
      molid_receptor = len(molids) - 1
    else:
      molid_receptor = None

    # Show samples
    if thickness is None:
      script += 'set ligand [mol new ' + self._FNs['prmtop']['L'] + ']\n'
      script += 'mol addfile ' + ligand_dcd_FN + ' type dcd waitfor all\n'
      script += 'mol drawframes $ligand 0 {0:%d}\n' % len(confs)
      if len(self.dir[process].split('/')) > 3:
        label = '/'.join(self.dir[process].split('/')[-3:])
        script += 'mol rename $ligand {%s}\n' % (label)
    else:
      for n in range(len(confs)):
        script += 'set ligand [mol new ' + self._FNs['prmtop']['L'] + ']\n'
        script += 'mol addfile ' + ligand_dcd_FN + \
          ' type dcd first %d last %d waitfor all\n'%(n,n)
        script += 'mol modstyle 0 $ligand CPK ' + \
          '%.6f %.6f 10.000000 10.000000\n'%(thickness[n]+0.2, thickness[n]+0.1)
        script += 'mol rename $ligand {%s}\n' % (n)
    molids.append('Samples')

    # Use the reference ligand structure as a basis for setting the view
    if show_ref_ligand:
      script += 'mol top $ref_ligand\n'
    elif show_starting_pose:
      script += 'mol top $start_ligand\n'
    script += 'display resetview\n'

    script += self.parse_view_args(view_args,\
      principal_axes_alignment=principal_axes_alignment, \
      nmolecules=len(molids), molid_receptor=molid_receptor)

    if save_image:
      image_path = os.path.join(self.dir['dock'], prefix + '.tga')
      if 'render' in view_args.keys():
        render = view_args['render']
      else:
        render = 'snapshot'
      script += 'render %s %s\n' % (render, image_path)
    if quit:
      script += 'quit\n'
    if execute:
      script_FN = os.path.join(self.dir[process], 'show_samples.vmd')
      script_F = open(script_FN, 'w')
      script_F.write(script)
      script_F.close()

      import os
      os.environ['VMDNOCUDA'] = "True"

      import subprocess
      if not 'vmd' in self._FNs.keys():
        self._FNs['vmd'] = a.findPaths(['vmd'])['vmd']

      vmd_args = [self._FNs['vmd']]
      if quit:
        vmd_args.extend(['-dispdev', 'text'])
      if 'size' in view_args.keys():
        vmd_args.extend(['-size',\
          '%d'%view_args['size'][0],'%d'%view_args['size'][1]])
      vmd_args.extend(['-e', script_FN])
      print vmd_args
      p = subprocess.Popen(vmd_args, \
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      (vmd_stdout, vmd_stderr) = p.communicate()
      p.wait()

      print 'script:'
      print script
      print 'stdout:'
      print vmd_stdout
      print 'stderr:'
      print vmd_stderr

      if clear_files:
        for FN in [ligand_dcd_FN, ref_ligand_dcd_FN, start_ligand_dcd_FN, \
                   complex_dcd_FN, script_FN]:
          if os.path.isfile(FN):
            os.remove(FN)

    if save_image and os.path.isfile(image_path) and (image_labels is
                                                      not None):
      import PIL.Image
      im = PIL.Image.open(image_path)
      import PIL.ImageDraw
      draw = PIL.ImageDraw.Draw(im)
      import PIL.ImageFont
      self._FNs['font'] = a.findPaths(['font'])['font']
      if self._FNs['font'] is not None:
        font = PIL.ImageFont.truetype(self._FNs['font'], size=84)
        textsize = draw.textsize(image_labels[0], font=font)
        ytop = im.size[1] - textsize[1] - 20
        draw.rectangle((5, ytop-5, 15+textsize[0], im.size[1] - 5), \
          fill='Black', outline='White')
        draw.text((10, ytop), image_labels[0], font=font)
        textsize = draw.textsize(image_labels[1], font=font)
        draw.rectangle((im.size[0]/2.+5, ytop-5, \
          im.size[0]/2.+15+textsize[0], im.size[1] - 5), \
          fill='Black', outline='White')
        draw.text((im.size[0] / 2. + 10, ytop), image_labels[1], font=font)
      else:
        textsize = draw.textsize(image_labels[0])
        ytop = im.size[1] - textsize[1] - 20
        draw.rectangle((5, ytop-5, 15+textsize[0], im.size[1] - 5), \
          fill='Black', outline='White')
        draw.text((10, ytop), image_labels[0])
        textsize = draw.textsize(image_labels[1])
        draw.rectangle((im.size[0]/2.+5, ytop-5, \
          im.size[0]/2.+15+textsize[0], im.size[1] - 5), \
          fill='Black', outline='White')
        draw.text((im.size[0] / 2. + 10, ytop), image_labels[0])
      im.save(image_path)

    return script

  def parse_view_args(self, view_args, principal_axes_alignment=False, \
      nmolecules=None, molid_receptor=None):
    script = 'display projection Orthographic\n'

    if ('rotate_matrix' in view_args.keys()) or principal_axes_alignment:
      if ('rotate_matrix' in view_args.keys()):
        script += 'set view_rotate_matrix %s\n' % view_args['rotate_matrix']
      elif principal_axes_alignment and (nmolecules is not None):
        script += vmd_principal_axes_procedure % nmolecules
        if molid_receptor is not None:
          script += vmd_unblock_procedure % molid_receptor
    elif 'rotate_by' in view_args.keys():
      script += 'rotate x by %d\n' % (view_args['rotate_by'][0])
      script += 'rotate y by %d\n' % (view_args['rotate_by'][1])
      script += 'rotate z by %d\n' % (view_args['rotate_by'][2])

    if 'scale_by' in view_args.keys():
      script += 'scale by %f\n' % view_args['scale_by']
    if 'scale_to' in view_args.keys():
      script += 'scale to %f\n' % view_args['scale_to']
    if 'nearclip' in view_args.keys():
      script += 'display nearclip set %f\n' % view_args['nearclip']
    if ('axes_off' in view_args.keys()) and (view_args['axes_off'] is True):
      script += 'axes location Off\n'
    return script

  def render_intermediates(self, process='dock', movie_name=None, \
        stride=1, nframes=None, view_args={}):
    # Determine state indices by the number of frames of by the stride
    if nframes is not None:
      state_inds = [int(s) \
        for s in np.linspace(0,len(self.dock_protocol)-1,nframes)]
    else:
      state_inds = range(0, len(self.confs[process]['samples']), stride)

    # Confirm that intermediate conformations exist
    for state_ind in state_inds:
      if self.confs['dock']['samples'][state_ind] == []:
        raise Exception('No snapshots in state %d' % state_ind)

    # Generate each snapshot
    view_args[
      'render'] = 'TachyonInternal'  # This works without the vmd display
    for s in range(len(state_inds)):
      # Format label for snapshot
      lambda_s = getattr(self, '%s_protocol' % process)[state_inds[s]]
      if process == 'dock':
        labels = [u'\u03B1 = %-8.2g' % lambda_s['a']]
      else:
        labels = []
      labels.append('T = %4.1f K' % lambda_s['T'])
      # Generate the snapshot
      self.show_samples(prefix=None, process=process, state=state_inds[s], \
        show_ref_ligand=True, show_starting_pose=False, show_receptor=True, \
        save_image=True, image_labels=labels, execute=True, quit=True, \
        view_args=view_args)

    # Make a movie
    if movie_name is not None:
      self._FNs['convert'] = a.findPaths(['convert'])['convert']
      if self._FNs['convert'] is not None:
        import subprocess
        subprocess.call([self._FNs['convert'],'-delay','15',\
          os.path.join(self.dir[process],'dock*.tga'), movie_name])
      else:
        raise Exception('ImageMagick convert required to make movies')

      import glob
      tga_FNs = glob.glob(os.path.join(self.dir[process], 'dock*.tga'))
      for tga_FN in tga_FNs:
        os.remove(tga_FN)
