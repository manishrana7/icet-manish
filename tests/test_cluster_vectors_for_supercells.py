from ase.build import bulk, cut

import itertools

from icetdev import ClusterSpace


# initialize cluster space and get the internal primitive atoms object
prim = bulk('Au', a=4.0, crystalstructure='hcp')
subelements = ['Au', 'Pd']
cutoffs = [0.0]
cs = ClusterSpace(prim, cutoffs, subelements)

atoms_prim = cs.get_primitive_structure().to_atoms()


# Create a supercell using permutation matrix
p_trial = [[1, 2, 1], [3, 1, 20], [0, 0, 2]]
atoms = cut(atoms_prim, p_trial[0], p_trial[1], p_trial[2])

# Setup cartesian input to generate a random population
cartesian_product_input = []
for i in range(len(atoms)):
    cartesian_product_input.append(['Pd', 'Au'])

# Loop over element combinations and assert expected singlet value
for subset in itertools.product(*cartesian_product_input):
    for atom, element in zip(atoms, subset):
        atom.symbol = element
    cv = cs.get_cluster_vector(atoms)
    expected_singlet = (-subset.count("Pd") + subset.count("Au")) / len(atoms)
    assert abs(
        cv[1] - expected_singlet) < 1e-3, "wrong singlet clustervector element"
