import numpy as np
from ase.neighborlist import NeighborList

import spglib
from icet.tools.geometry import (get_fractional_positions_from_ase_neighbor_list,
                                 get_primitive_structure)


class PermutationMatrix(object):
    '''
    Permutation matrix object
    '''

    def __init__(self, atoms, cutoff, find_prim=True, verbosity=0):
        atoms = atoms.copy()
        self.cutoff = cutoff
        self.verbosity = verbosity
        # set each element to the same since we only care about geometry when
        # taking primitive

        atoms.set_chemical_symbols(len(atoms) * [atoms[0].symbol])

        if find_prim:
            atoms = get_primitive_structure(atoms)

        self.primitive_structure = atoms

        if self.verbosity >= 3:
            print('size of atoms {}'.format(len(self.primitive_structure)))

        # Get symmetry information and load into a permutation map object
        symmetry = spglib.get_symmetry(self.primitive_structure)
        self.translations = symmetry['translations']
        self.rotations = symmetry['rotations']
        self.build()

    def build(self):
        """
        Builds the permutation matrix
        """
        # Create neighbor_lists from the different cutoffs
        neighbor_list = NeighborList(len(
            self.primitive_structure) * [self.cutoff / 2.0], skin=1e-8,
            bothways=True, self_interaction=False)
        neighbor_list.update(self.primitive_structure)

        # get fractional positions for neighbor_list
        frac_positions = get_fractional_positions_from_ase_neighbor_list(
            self.primitive_structure, neighbor_list)
        
        if self.verbosity >= 3:
            print('number of positions: {}'.format(len(frac_positions)))

        

        permutation_matrix = []
        for frac_pos in frac_positions:
            permutation_row = [frac_pos]
            # permutation_row = []
            for rotation, translation in zip(self.rotations, self.translations):
                permutated_position = translation + np.dot(frac_pos, rotation)
                permutation_row.append(permutated_position)
            permutation_matrix.append(permutation_row)
        
        self.permutation_matrix = permutation_matrix