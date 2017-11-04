'''
This scripts checks the computation of cluster vectors for three body centerd
cubic based structures.
'''

from icetdev.clusterspace import create_clusterspace
from icetdev.structure import structure_from_atoms
from ase.build import bulk, make_supercell
import numpy as np

cutoffs = [8.0, 7.0]
subelements = ['W', 'Ti']

print('')
prototype = bulk('W')
cs = create_clusterspace(prototype, cutoffs, subelements)

# structure #1
print(' structure #1')
conf = structure_from_atoms(prototype.copy())
cv = cs.get_clustervector(conf)
cv_target = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                      1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                      1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                      1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                      1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                      1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
assert np.all(np.abs(cv_target - cv) < 1e-6)

# structure #2
print(' structure #2')
conf = make_supercell(prototype, [[2, 0, 1],
                                  [0, 1, 0],
                                  [0, 1, 2]])
conf[0].symbol = 'Ti'
conf[1].symbol = 'Ti'
conf = structure_from_atoms(conf)
cv = cs.get_clustervector(conf)
cv_target = np.array([1.0, 0.0, 0.0, -0.3333333333333333,
                      -0.3333333333333333, 0.0, 1.0, 1.0,
                      0.0, -0.3333333333333333,
                      -0.3333333333333333, 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
assert np.all(np.abs(cv_target - cv) < 1e-6)

# structure #3
print(' structure #3')
conf = make_supercell(prototype, [[1,  0, 1],
                                  [0,  1, 1],
                                  [0, -1, 3]])
conf[0].symbol = 'Ti'
conf[1].symbol = 'Ti'
conf[2].symbol = 'Ti'
conf = structure_from_atoms(conf)
cv = cs.get_clustervector(conf)
cv_target = np.array([1.0, -0.5, 0.0, 0.6666666666666666,
                      0.3333333333333333, 0.0, 0.0, 1.0,
                      0.0, 0.6666666666666666, 0.3333333333333333,
                      -0.16666666666666666, 0.5, 0.16666666666666666,
                      -0.5, -0.16666666666666666, -0.5,
                      0.16666666666666666, -0.5, -0.5,
                      0.16666666666666666, 0.5, -0.16666666666666666,
                      -0.5, 0.16666666666666666, 0.16666666666666666,
                      -0.5, 0.5, -0.16666666666666666, 0.5,
                      0.16666666666666666, -0.16666666666666666,
                      -0.5, 0.5, 0.16666666666666666, -0.5])
assert np.all(np.abs(cv_target - cv) < 1e-6)
