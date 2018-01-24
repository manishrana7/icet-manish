import numpy as np
from ase.neighborlist import NeighborList

import spglib
from icet.tools.geometry import (get_fractional_positions_from_ase_neighbor_list,
                                 find_lattice_site_by_position,
                                 get_primitive_structure,
                                 fractional_to_cartesian)


class PermutationMatrix(object):
    '''
    Permutation matrix object

    1: get frac position version of PM
    2: construct the lattice site version
    3: prune the lattice site version by only
       having unique lattice sites in column 1
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
        self.build_lattice_site_permutation_matrix()
        # self.prune_permutation_matrix()

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

        self.permutaded_matrix_frac = permutation_matrix

    def build_lattice_site_permutation_matrix(self):
        """
        Builds the lattice site version of permutation matrix.

        It will loop through all positions, pos_ij in the matrix,
        convert to lattice sites and construct the matrix,
        lattice_site_ij
        """
        pm_lattice_sites = []
        for row in self.permutaded_matrix_frac:
            positions = fractional_to_cartesian(self.primitive_structure, row)
            sites = []
            if np.all(self.primitive_structure.pbc):
                sites = [find_lattice_site_by_position(
                    self.primitive_structure, position) for position in positions]
            else:
                for position in positions:
                    try:
                        lattice_site = find_lattice_site_by_position(
                            self.primitive_structure, position)
                    except:
                        pass
            if len(sites) > 0:
                pm_lattice_sites.append(sites)
            else:
                raise Exception("Lattice sites are empty")


        self.pm_lattice_sites = pm_lattice_sites
            
