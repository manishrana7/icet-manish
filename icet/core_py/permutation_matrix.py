import numpy as np
from ase.neighborlist import NeighborList

import spglib
from icet.tools.geometry import (
    get_fractional_positions_from_ase_neighbor_list,
    find_lattice_site_by_position,
    get_primitive_structure,
    fractional_to_cartesian)


class PermutationMatrix(object):
    '''
    Permutation matrix object.

    parameters
    ----------
    atoms : ASE Atoms object
    cutoff : float
        this sets the radius for the the sites
        to include in the permutation matrix.
    find_prim : bool (default True)
        if True the incoming atoms object will
        be used to construct a primitive structure
        via spglib. This primitive structure will
        then be used to construct the neighbor list,
        be what the lattice sites refer to etc.
    verbosity : int
        sets the verbosity of this class
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
        self._build_lattice_site_permutation_matrix()
        self._prune_permutation_matrix()
        self._sort()

    def build(self):
        """
        Builds the permutation matrix.
        """
        # Create neighbor_lists from the different cutoffs
        neighbor_list = NeighborList(len(
            self.primitive_structure) * [self.cutoff / 2.0], skin=1e-8,
            bothways=True, self_interaction=False)
        neighbor_list.update(self.primitive_structure)

        # get fractional positions for neighbor_list
        frac_positions = get_fractional_positions_from_ase_neighbor_list(
            self.primitive_structure, neighbor_list)

        # frac_positions.sort()
        if self.verbosity >= 3:
            print('number of positions: {}'.format(len(frac_positions)))

        permutation_matrix = []
        for frac_pos in frac_positions:
            # permutation_row = [frac_pos]
            permutation_row = []
            for rotation, translation in zip(
                    self.rotations, self.translations):
                permutated_position = translation + \
                    np.dot(frac_pos, rotation.T)
                permutation_row.append(permutated_position)
            permutation_matrix.append(permutation_row)

        self.permutaded_matrix_frac = permutation_matrix

    def _build_lattice_site_permutation_matrix(self):
        """
        Builds the lattice site version of permutation matrix.

        It will loop through all positions, pos_ij in the matrix,
        convert to lattice sites and construct the matrix,
        lattice_site_ij

        This is a helper method for the build method.
        """
        pm_lattice_sites = []
        for row in self.permutaded_matrix_frac:
            positions = fractional_to_cartesian(self.primitive_structure, row)
            sites = []
            if np.all(self.primitive_structure.pbc):
                sites = [find_lattice_site_by_position(
                    self.primitive_structure, position)
                    for position in positions]
            else:
                for position in positions:
                    try:
                        lattice_site = find_lattice_site_by_position(
                            self.primitive_structure, position)
                        sites.append(lattice_site)
                    except:  # NOQA
                        pass
            if len(sites) > 0:
                pm_lattice_sites.append(sites)
            else:
                raise Exception("Lattice sites are empty")

        self.pm_lattice_sites = pm_lattice_sites

    def _prune_permutation_matrix(self):
        """
        Prunes columns in the permutation matrix.
        This is a helper method for the build method.
        """
        for i in range(len(self.pm_lattice_sites)):
            for j in reversed(range(len(self.pm_lattice_sites))):
                if j <= i:
                    continue
                if self.pm_lattice_sites[i][0] == self.pm_lattice_sites[j][0]:
                    self.pm_lattice_sites.pop(j)
                    if self.verbosity > 2:
                        msg = ['Removing duplicate in permutation matrix']
                        msg += ['i: {} j: {}'.format(i, j)]
                        print(' '.join(msg))

    def _sort(self):
        """
        Sorts the Lattice Site permutation matrix.
        This is a helper method for the build method.
        """
        self.pm_lattice_sites.sort()
