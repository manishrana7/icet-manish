import numpy as np
from ase.db import connect
from ase.build import bulk
from icetdev.cluster_space import ClusterSpace
from ase.visualize import view
from ase.io import write
from spglib import get_spacegroup
from spglib import niggli_reduce
from icetdev.permutation_map import __get_primitive_structure
from icetdev.tools.geometry import transform_cell_to_cell

from icetdev.tools.map_sites import map_configuration_to_reference

prim = bulk('Au', a=4.0, crystalstructure='hcp')
subelements = ['Au', 'Pd']
cutoffs = [8.0, 6.0]
cs = ClusterSpace(prim, cutoffs, subelements)

atoms_prim = cs.get_primitive_structure().to_atoms()
# view(atoms_prim)
db = connect('hcp-equivalent.db')
for entry in db.select():
    atoms = entry.toatoms()
    view(atoms)
    atoms = transform_cell_to_cell(atoms, atoms_prim)
    # atoms, dr_max, dr_avg = map_configuration_to_reference(atoms, atoms_prim, tolerance_mapping = 1.0)
    # atoms.cell = atoms_prim.cell
    # atoms.wrap()
    view(atoms)