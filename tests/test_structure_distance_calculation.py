import numpy as np
from ase.db import connect
import icetdev as icet
from ase.neighborlist import NeighborList

'''
Testing the calculation of distances with offsets
'''


''' Tolerance for comparing distances '''
DISTTOL = 1e-8

''' Fetch structures from database '''
db = connect('structures_for_testing.db')

print('')
for row in db.select():
    print(' structure: {}'.format(row.tag))
    atoms_row = row.toatoms()
    structure = icet.structure_from_atoms(atoms_row)
    nl = NeighborList(len(atoms_row)*[2.6])
    nl.update(atoms_row)
    for index in range(len(atoms_row)):
        indices, offsets = nl.get_neighbors(index)
        for i, offset in zip(indices, offsets):
            dvec = atoms_row.positions[index] - atoms_row.positions[i]
            dvec -= np.dot(offset, atoms_row.get_cell())
            dist_ase = np.linalg.norm(dvec)
            dist_struct = structure.get_distance2(index, [0, 0, 0], i, offset)
            msg = 'Testing distance calculator failed'
            msg += ' for structure {}'.format(row.tag)
            assert dist_ase - dist_struct < DISTTOL, msg
