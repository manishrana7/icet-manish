import spglib
from _icet import PermutationMap
from .neighbor_list import NeighborList
from .structure import Structure
from ..tools.geometry import (get_primitive_structure,
                              get_fractional_positions_from_neighbor_list)

from icet.io.logging import logger
logger = logger.getChild('permutation_map')


def permutation_matrix_from_atoms(atoms, cutoff, find_prim=True):
    '''Set up a list of permutation maps from an atoms object.

    Parameters
    ----------
    atoms : ASE Atoms object
        input structure
    cutoff : float
        cutoff radius
    find_primitive : boolean
        if True the symmetries of the primitive structure will be employed

    Returns
    -------
    matrix, icet Structure object, NeighborList object
        the tuple comprises the permutation matrix, the primitive structure,
        and the neighbor list
    '''

    atoms = atoms.copy()
    # set each element to the same since we only care about geometry when
    # taking primitive
    if len(atoms) > 0:
        atoms.set_chemical_symbols(len(atoms) * [atoms[0].symbol])
    else:
        raise Exception('Length of atoms is {}'.format(len(atoms)))

    atoms_prim = atoms
    if find_prim:
        atoms_prim = get_primitive_structure(atoms)

    logger.debug('Size of primitive structure: {}'.format(len(atoms_prim)))

    # get symmetry information
    symmetry = spglib.get_symmetry(atoms_prim)
    translations = symmetry['translations']
    rotations = symmetry['rotations']

    # set up a permutation map object
    permutation_matrix = PermutationMap(translations, rotations)

    # create neighbor_lists from the different cutoffs
    prim_structure = Structure.from_atoms(atoms_prim)
    neighbor_list = NeighborList(cutoff)
    neighbor_list.build(prim_structure)

    # get fractional positions for neighbor_list
    frac_positions = get_fractional_positions_from_neighbor_list(
        prim_structure, neighbor_list)
    # frac_positions.sort()

    logger.debug('Number of fractional positions: {}'.format(len(frac_positions)))
    if len(frac_positions) > 0:
        permutation_matrix.build(frac_positions)

    return permutation_matrix, prim_structure, neighbor_list
