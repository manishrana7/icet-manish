from ase.build import bulk
from ase.build import cut
from ase.build import niggli_reduce

import itertools

from icetdev import ClusterSpace
from icetdev.tools.geometry import get_primitive_structure
import unittest

import numpy as np

tc = unittest.TestCase('__init__')

# initialize cluster space and get the internal primitive atoms object
prim = bulk('Au', a=4.0, crystalstructure='hcp')
subelements = ['Au', 'Pd']
cutoffs = [0.0]
cs = ClusterSpace(prim, cutoffs, subelements)

atoms_prim = cs.get_primitive_structure().to_atoms()


# Create a supercell using permutation matrix
p_trial = [[1, 0, 0], [0, 1, 5], [0, 0, 2]]
atoms2 = cut(atoms_prim, p_trial[0], p_trial[1], p_trial[2])

print("det:",np.linalg.det(p_trial))
# exit(1)
# Setup cartesian input to generate a random population
cartesian_product_input = []
for i in range(len(atoms2)):
    cartesian_product_input.append(['Pd', 'Au'])

# Loop over element combinations and assert expected singlet value
for subset in itertools.product(*cartesian_product_input):    
    for atom, element in zip(atoms2, subset):
        atom.symbol = element
    # atoms = get_primitive_structure(atoms2)
    # niggli_reduce(atoms)
    atoms=atoms2
    cv = cs.get_cluster_vector(atoms)
    expected_singlet = (-atoms.get_chemical_symbols().count("Pd") + atoms.get_chemical_symbols().count("Au")) / len(atoms)
    # print("testRes: ",cv[1], expected_singlet, len(atoms),  atoms.get_chemical_symbols().count("Pd"), atoms.get_chemical_symbols().count("Au"))
    # print(atoms.cell)
    # print("cell det",np.linalg.det(atoms.cell))
    tc.assertAlmostEqual(cv[1], expected_singlet)
    # assert abs(
        # cv[1] - expected_singlet) < 1e-3, "wrong singlet clustervector element"
