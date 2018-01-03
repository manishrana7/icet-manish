'''
This example demonstrates how to count the number of clusters for a structure.
'''
# Import modules
import numpy as np

from ase.build import bulk

from icetdev.cluster_counts import ClusterCounts
from icetdev.structure import Structure
from icetdev.orbit_list import create_orbit_list
from icetdev.permutation_map import __get_primitive_structure

# Create a titanium, single-layered, sheet and randomly populate some of the
# sites with W atoms. Also set the cutoffs for pairs to 4 Å.
atoms = bulk("Ti", "bcc", a=3.43).repeat([2, 2, 1])
atoms.pbc = [True, True, False]

atoms.set_chemical_symbols([['Ti', 'W'][n] for n in
                            np.round(np.random.random((len(atoms),
                                                       ))).astype(int)])
cutoffs = [4]

# Create the orbit list that corresponds to the primitive structure.
prim_atoms = __get_primitive_structure(atoms)
prim_structure = Structure.from_atoms(prim_atoms)
prim_orbitlist = create_orbit_list(prim_structure, cutoffs)

# Use the primitive orbitlist to count the number of clusters.
cluster_counts = ClusterCounts()
structure = Structure.from_atoms(atoms)
cluster_counts.count_clusters(structure, prim_orbitlist)

# Print all of the clusters that were found.
print("number of atoms {0}".format(len(atoms)))
print("Found {} clusters".format(len(cluster_counts)))
cluster_counts.print()
