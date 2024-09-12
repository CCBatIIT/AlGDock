import numpy as np
from munkres import Munkres


class hRMSD():
  """
  A class to compute the Hungarian symmetry-corrected heavy-atom
  root mean square deviation between a set of configurations
  and a reference configuration.

  See http://dock.compbio.ucsf.edu/DOCK_6/dock6_manual.htm#LigandRMSD
  """
  def __init__(self, prmtop_FN, inv_prmtop_atom_order, ref_conf=None):
    """
    prmtop_FN is an AMBER prmtop file
    inv_prmtop_atom_order is a list of indices that maps
      from the prmtop order to the MMTK order
    ref_conf is the default reference configuration,
      an Nx3 array in MMTK order
    """
    import AlGDock.IO
    IO_prmtop = AlGDock.IO.prmtop()
    ligand_prmtop = IO_prmtop.read(prmtop_FN,\
      ['ATOMIC_NUMBER','ATOM_TYPE_INDEX'])

    atom_type_indices = ligand_prmtop['ATOM_TYPE_INDEX'][inv_prmtop_atom_order]
    atom_numbers = ligand_prmtop['ATOMIC_NUMBER'][inv_prmtop_atom_order]
    atom_indices = np.array(range(len(atom_type_indices)))

    self.atom_sets_to_compare = []
    for j in set(atom_type_indices):
      indices_j = (atom_type_indices == j)
      if (atom_numbers[indices_j] != 1).any():
        self.atom_sets_to_compare.append((sum(indices_j),\
                                          atom_indices[indices_j]))
    self.nheavy_atoms = np.sum(atom_numbers != 1)
    self.ref_conf = ref_conf
    self.munkres = Munkres()

  def __call__(self, confs, ref_conf=None):
    """
    confs is a list of Nx3 arrays in MMTK order.
    ref_conf is a Nx3 array in MMTK order.
    """
    rmsds = []
    if ref_conf is None:
      ref_conf = self.ref_conf
    for conf in confs:
      ssd = 0.
      for (nelements, atom_set) in self.atom_sets_to_compare:
        if nelements == 1:
          j = atom_set[0]
          ssd += np.sum(np.square(conf[j, :] - ref_conf[j, :]))
        else:
          cost_matrix = np.array([[\
            np.sum(np.square(conf[atom_set[j],:]-ref_conf[atom_set[k],:])) \
              for j in range(nelements)] \
                for k in range(nelements)])
          path = self.munkres.compute(cost_matrix)
          ssd += np.sum([np.sum(np.square(\
            conf[atom_set[j],:]-ref_conf[atom_set[k],:])) for (j,k) in path])
      rmsds.append(np.sqrt(ssd / self.nheavy_atoms))
    return rmsds

  def set_ref_configuration(self, ref_conf):
    """
    Sets a new reference configuration
    """
    self.ref_conf = ref_conf
