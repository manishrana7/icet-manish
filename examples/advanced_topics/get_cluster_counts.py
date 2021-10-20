"""
This example demonstrates how to count the number of clusters for a structure.
"""

# Start import
from ase.build import bulk
from icet.core.orbit_list import OrbitList
# End import

# Create a titanium, single-layered, sheet and randomly populate some of the
# sites with W atoms.
# Start setup
prim_structure = bulk('Ti', 'sc', a=3.0)
structure = prim_structure.repeat([2, 1, 1])
structure.set_chemical_symbols(['Ti', 'W'])
cutoffs = [5.0]
# End setup

# Determine the orbit list for the corresponding primitive structure for all
# pair clusters within the cutoff distance
symprec = 1e-5  # tolerance used by spglib
position_tolerance = 1e-5  # tolerance used when comparing positions
fractional_position_tolerance = position_tolerance / 3  # ... in fractional coordinates
prim_orbitlist = OrbitList(prim_structure, cutoffs, symprec,
                           position_tolerance, fractional_position_tolerance)
# Use the primitive orbit list to count the number of clusters.
cluster_counts = prim_orbitlist.get_cluster_counts(structure, fractional_position_tolerance)
# Print all of the clusters that were found.
print('Number of atoms: {0}'.format(len(structure)))
print('Found {} orbits'.format(len(cluster_counts)))
orbit_names = {1: '(singlet)', 2: '(pair)', 3: '(triplet)', 4: '(quadruplet)'}
for orbit_index, counts in cluster_counts.items():
    orbit = prim_orbitlist.get_orbit(orbit_index)
    print()
    print(f'Orbit index: {orbit_index}')
    print(f'Number of atoms in cluster: {orbit.order} {orbit_names.get(orbit.order, "")}')
    if orbit.order > 1:
        print(f'Cluster radius: {orbit.radius:.4f}')
        print('Distances between atoms in the clusters: ', end='')
        print(' '.join(f'{d:.4f}' for d in orbit.representative_cluster.distances))
    for symbols, count in counts.items():
        print(' '.join(f'{sym:3s}' for sym in symbols), f'{count:5d}')
