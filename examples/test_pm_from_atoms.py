from ase import Atoms
from ase.build import bulk
from icetdev.permutationMap import PermutationMap, permutation_maps_from_atoms
import numpy as np

atoms = bulk("Al", "fcc", a=2.0).repeat(1)

atoms = bulk("Ti", "bcc", a=3.3).repeat(5)

atoms.set_chemical_symbols([['Ti', 'W'][n] for n in np.round(
    np.random.random((len(atoms),))).astype(int)])


cutoffs = [30, 10, 10, 5, 1]

permutation_maps_from_atoms(atoms, cutoffs, verbosity=0)
