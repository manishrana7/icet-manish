import spglib
from icet.tools.geometry import (get_primitive_structure)
from ase.neighborlist import NeighborList


class PermutationMatrix(object):
    '''
    Permutation matrix object
    '''

    def __init__(self, atoms, cutoff, find_prim=True, verbosity=0):
        atoms = atoms.copy()
        # set each element to the same since we only care about geometry when
        # taking primitive

        atoms.set_chemical_symbols(len(atoms) * [atoms[0].symbol])

        if find_prim:
            atoms = get_primitive_structure(atoms)

        self.primitive_structure = atoms

        if verbosity >= 3:
            print('size of atoms {}'.format(len(self.primitive_structure)))

        # Get symmetry information and load into a permutation map object
        symmetry = spglib.get_symmetry(self.primitive_structure)
        self.translations = symmetry['translations']
        self.rotations = symmetry['rotations']

        # Create neighbor_lists from the different cutoffs
        neighbor_list = NeighborList(len(atoms) * [cutoff / 2.0], skin=1e-8,
                                     bothways=True, self_interaction=False)
        neighbor_list.update(self.primitive_structure)

        # get fractional positions for neighbor_list
        # frac_positions = get_fractional_positions_from_neighbor_list(
        # self.primitive_structure, neighbor_list)
        # if verbosity >= 3:
        #     print('number of positions: {}'.format(len(frac_positions)))
        # if len(frac_positions) > 0:
        #     self.build(frac_positions)

    def build(self, frac_positions):
        """
        Builds the permutation matrix
        """
        pass
