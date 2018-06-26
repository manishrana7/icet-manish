from icet import ClusterSpace
from ase.build import bulk
from ase.db import connect
import numpy as np
import os

prim = bulk('Au', a=4.0, crystalstructure='hcp')
cutoffs = [7.0, 7.0, 7.0]
subelements = ['Au', 'Pd']
cs = ClusterSpace(prim, cutoffs, subelements)

db = connect(os.path.dirname(__file__) + '/equivalent_structure_pairs.db')

# Loop over all pairs
for structure in db.select():
    # Do not check the pair that was just checked
    if structure.equivalent_structure < structure.id:
        continue

    cv_1 = cs.get_cluster_vector(structure.toatoms())
    cv_2 = cs.get_cluster_vector(
        db.get(structure.equivalent_structure).toatoms())

    assert np.all(np.abs(cv_2 - cv_1) < 1e-6)
