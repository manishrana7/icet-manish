"""
This scripts checks the computation of cluster vectors for three BCC-based
structures.
"""

from icetdev import clusterspace
from icetdev.clusterspace import create_clusterspace
from icetdev.structure import structure_from_atoms
from ase.build import bulk, make_supercell
import numpy as np

print(__doc__)

cutoffs = [8, 7]
subelements = ['W', 'Ti']

prototype = bulk('W')


clusterspace = create_clusterspace(prototype, cutoffs, subelements)

conf = structure_from_atoms(prototype.copy())

print('Structure no. 1 (nat= {}):'.format(len(conf)))
cv = clusterspace.get_clustervector(conf)
cv_target = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                     1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                     1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                     1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                     1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                     1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
print(cv)
assert np.all(np.abs(cv_target - cv) < 1e-6)

conf = make_supercell(prototype, [[2, 0, 1],
                                  [0, 1, 0],
                                  [0, 1, 2]])

conf[0].symbol = 'Ti'
conf[1].symbol = 'Ti'
conf = structure_from_atoms(conf)
print('Structure no. 2 (nat= {}):'.format(len(conf)))
cv = clusterspace.get_clustervector(conf)

cv_target = np.array([1.0, 0.0, 0.0, -0.3333333333333333,
                     -0.3333333333333333, 0.0, 1.0, 1.0,
                     0.0, -0.3333333333333333, 
                     -0.3333333333333333, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
print(cv)
assert np.all(np.abs(cv_target - cv) < 1e-6)

conf = make_supercell(prototype, [[1,  0, 1],
                                  [0,  1, 1],
                                  [0, -1, 3]])

conf[0].symbol = 'Ti'
conf[1].symbol = 'Ti'
conf[2].symbol = 'Ti'
print('Structure no. 3 (nat= {}):'.format(len(conf)))
conf = structure_from_atoms(conf)
cv = clusterspace.get_clustervector(conf)
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
print(cv)
assert np.all(np.abs(cv_target - cv) < 1e-6)
