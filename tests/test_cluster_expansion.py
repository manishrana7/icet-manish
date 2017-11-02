"""
This scripts checks that cluster expansion model can be initialized with any structure
in the test database and can predict a property

"""

import numpy as np

from ase.db import connect
from icetdev.clusterspace import create_clusterspace
from icetdev.cluster_expansion import ClusterExpansion
print(__doc__)


def test_clusterexpansion_model(atoms, cutoffs, subelements):
    """
    Test clusterexpansion init and prediction
        :param atoms: ASE Atoms object
        :param subelements: list of atomic symbols
    """
    cs = create_clusterspace(atoms, cutoffs, subelements)
    params_len = cs.get_clusterspace_size() + 1  # plus one for singlet
    params = np.random.rand(params_len)

    ce = ClusterExpansion(cs, params)
    predicted_val = ce.predict(atoms)

    assert isinstance(predicted_val, float)


db = connect("structures_for_testing.db")
subelements = ["H", "He", "Pb"]

for row in db.select():

    atoms_row = row.toatoms()
    cutoffs = [1.4] * 3

    print("Testing cluster expansion for structure: {} with cutoffs {}".format(
        row.tag, cutoffs))

    for at in atoms_row:
        at.symbol = subelements[0]
    test_clusterexpansion_model(atoms_row, cutoffs, subelements)
