'''This module has the purpose of enumerating structures. Given a lattice
(possibly with as basis) and a number of elements, the code generates all
the derivative superstructures having a certain size defined by the user.
'''

from spglib import get_symmetry
import numpy as np
from ase import Atoms
from icetdev.enumeration.hermite_normal_form import get_reduced_hnfs


def dehash_labelkey(labelkey, natoms, nelements):
    '''
    Calculate labeling from hashkey.

    Parameters
    ----------
    labelkey : int
        Hash key to a labeling.
    natoms : int
        Number of atoms in the labeling.
    nelements : int
        Number of elements in the enumeration.

    Returns
    -------
    tuple of ints
        Labeling corresponding to labelkey.
    '''
    labeling = ()
    for _ in range(natoms):
        labeling += (labelkey % nelements,)
        labelkey = labelkey // nelements
    return labeling


def hash_labeling(labeling, natoms, nelements):
    '''
    Calculate hash.

    Parameters
    ----------
    labelkey : int
        Hash key to a labeling.
    natoms : int
        Number of atoms in the labeling.
    nelements : int
        Number of elements in the enumeration.

    Returns
    -------
    tuple of ints
        Labeling corresponding to labelkey.
    '''
    labelkey = 0
    for i in range(natoms):
        labelkey += labeling[i] * nelements**i
    return labelkey


def get_unique_snfs(hnfs):
    '''
    For a list of Hermite Normal Forms, obtain the set of unique Smith Normal
    Forms.

    Paramters
    ---------
    hnfs : list of HermiteNormalForm objects

    Returns
    -------
    list of SmithNormalForm objects
        The unique Smith Normal Form matrices.
    '''
    snfs = []
    for hnf in hnfs:
        # Check whether the snf is new or already encountered
        snf = hnf.snf
        snf_is_new = True
        for snf_comp in snfs:
            if snf_comp.S == snf.S:
                snf_is_new = False
                snf_comp.add_hnf(hnf)
                break
        if snf_is_new:
            snf.group_order = get_group_order(snf.S)
            snf.add_hnf(hnf)
            snfs.append(snf)
    return snfs


def translation_permutations(labeling, snf, nsites, nelements, include_self=False):
    '''
    Yield labelings that are equivalent to original labeling
    under translations as dictated by snf.

    Paramters
    ---------
    original : tuple
        labelings to be translated
    snf : SmithNormalForm object

    include_self : bool
        inlcude original labeling or not

    Yields
    ------
    tuple
        labeling that is equivalent to original
    '''

    # Compute size of each block within which translations occur
    size_1 = nsites * snf.blocks[0]
    size_2 = nsites * snf.blocks[1]
    size_3 = nsites * snf.blocks[2]

    # Loop over each possible translation
    for trans_1 in range(snf.S[0]):
        for trans_2 in range(snf.S[1]):
            for trans_3 in range(snf.S[2]):
                if not include_self and trans_1 + trans_2 + trans_3 == 0:
                    continue

                # Calculate the positions to which the atoms will be moved
                new_positions_1 = [(i + trans_1) % snf.S[0]
                                   for i in range(snf.S[0])]
                new_positions_2 = [(i + trans_2) % snf.S[1]
                                   for i in range(snf.S[1])]
                new_positions_3 = [(i + trans_3) % snf.S[2]
                                   for i in range(snf.S[2])]

                labelkey = 0
                count = 0

                for i in new_positions_1:
                    block_1 = labeling[size_1 * i:size_1 * (i + 1)]
                    for j in new_positions_2:
                        block_2 = block_1[size_2 * j:size_2 * j + size_2]
                        for k in new_positions_3:
                            for label in block_2[size_3 * k:size_3 * k + size_3]:
                                labelkey += label * nelements**count
                                count += 1
                yield labelkey


def get_group_order(snf):
    '''
    Get group represantation of an SNF matrix (the G matrix in HarFor08).

    Paramters
    ---------
    snf : tuple
        Diagonal of matrix on Smith Normal Form

    Returns
    -------
    ndarray
        Group representation, shape 3xN where N is product of elements in snf.
    '''
    group_order = []
    for i in range(snf[0]):
        for j in range(snf[1]):
            for k in range(snf[2]):
                group_order.append([i, j, k])
    return np.array(group_order)


def get_labelkeys(snf, nelements, nsites):
    '''
    Get all labelings corresponding to a Smith Normal Form. Superperiodic
    labelings as well as labelings that are equivalent for this particular SNF
    will not be included. However, labelings that are equivalent by rotations
    that leave the cell (but not the labeling) unchanged will still be
    included, since these have to be removed for each HNF.

    Paramters
    ---------
    snf : tuple
        The three diagonal elements of an SNF matrix.
    sites : int
        Number of sites per primitive cell.

    Returns
    -------
    list of tuples
        Symmetrically inequivalent labelings
    '''
    ncells = snf.N
    natoms = ncells * nsites
    labelkey_tracker = [False] * nelements**natoms
    labelkeys = []
    for labelkey in range(nelements**natoms):
        labeling = dehash_labelkey(labelkey, natoms, nelements)
        unique = True

        for labelkey_trans in translation_permutations(labeling, snf, nsites,
                                                       nelements,
                                                       include_self=False):
            if labelkey == labelkey_trans:
                unique = False
                break

            # Check with previous labelings,
            # if labeling can be translated into a previously
            # added labeling, then it is not unique
            if labelkey_tracker[labelkey_trans]:
                unique = False
                break
            if not unique:
                break
        if unique:
            labelkeys.append(labelkey)
            labelkey_tracker[labelkey] = True
    return labelkeys


def permute_labeling(labeling, snf, transformation, nsites, nelements):
    '''
    Rotate labeling based on group representation defined by Gp.

    Paramters
    ---------
    labeling : tuple
        Labeling to be rotated
    snf : tuple
        Diagonal of matrix on Smith Normal Form
    Gp : ndarray
        New group representation of labeling in SNF

    Returns
    -------
    tuple
        Labeling rotated based on Gp.
    '''

    # Calculate transformation imposed by LRL multiplication
    new_group_order = np.dot(snf.group_order, transformation[0].T)

    # Loop over every atom to find its new position
    labelkey = 0
    for member_index, member in enumerate(new_group_order):

        # Transform according to Gp,
        # but each site also transforms in its own way
        for basis in range(nsites):
            new_cell = member + transformation[1][basis]

            # Calculate new index, first by finding the right block,
            # then the basis index in that block
            new_index = 0
            for i in range(3):
                new_index += (new_cell[i] % snf.S[i]) * \
                    snf.blocks[i] * nsites
            new_index += transformation[2][basis]

            # Add the contribution to the hash key
            element = labeling[member_index * nsites + basis]
            labelkey += element * nelements**new_index
    return labelkey


def yield_unique_labelings(labelings, snf, hnf, nsites, nelements):
    '''
    Yield labelings that are unique in every imaginable sense.

    Parameters
    ----------
    labelings : list of tuples
        List of labelings that may still contain labelings that are equivalent
        under rotations that leaves the supercell shape unchanged.
    snf : tuple
        Diagonal elements of matrix on Smith Normal Form.
    hnf_rots : list of ndarrays
        Rotations that leaves the supercell unchanged but not necessarily the
        labeling.
    G : ndarray
        Group representation corresponding to snf.
    sites : int
        Number of sites in the primitive cell.

    Yields
    ------
    tuple
        Labeling, each and every one unique.
    '''
    natoms = snf.N * nsites
    labelkey_tracker = [False] * nelements**natoms
    for labelkey in labelings:
        labeling = dehash_labelkey(labelkey, natoms, nelements)

        # Check whether labeling is just a rotated version of a previous
        # labeling. Apply transformation that is specific to the hnf
        # and check all translations of the transformed labeling.
        unique = True
        for transformation in hnf.transformations:

            labelkey_rot = permute_labeling(labeling, snf, transformation,
                                            nsites, nelements)

            # Commonly, the transformation leaves the labeling
            # unchanged, so check that first as a special case
            # (yields a quite significant speedup)
            if labelkey_rot == labelkey:
                continue

            labeling_rot = dehash_labelkey(labelkey_rot, natoms, nelements)

            # Translate in all possible ways
            for labelkey_rot_trans in \
                    translation_permutations(labeling_rot, snf, nsites,
                                             nelements, include_self=True):
                if labelkey_tracker[labelkey_rot_trans]:
                    # Then we have rotated and translated the labeling
                    # into one that was already yielded
                    unique = False
                    break
            if not unique:
                break
        if unique:
            # Then we have finally found a unique structure
            # defined by an HNF matrix and a labeling
            yield labeling
            labelkey_tracker[labelkey] = True


def get_symmetry_operations(atoms, basis):
    '''
    Use spglib to calculate the symmetry operations of atoms and return their
    inverse matrices. basis_shifts correspond to d_N,d, rotation_translations t_N,d

    Parameters
    ----------
    atoms : ASE Atoms
        Structure for which the symmetry operations are sought.

    Returns
    -------
    list of ndarrays
        Inverse of matrices for rotation operatios of atoms.
    '''

    symmetries = get_symmetry(atoms)

    rotations = symmetries['rotations']
    translations = symmetries['translations']

    # Calculate how atoms within the primitive cell are
    # shifted upon operation with rotation matrix
    basis_shifts = np.zeros((len(rotations), len(basis)), dtype='int64')
    rotation_translations = []
    for i, rotation in enumerate(rotations):
        translation = translations[i]
        rotation_translation = []
        for j, basis_element in enumerate(basis):
            rot_trans = np.dot(rotation, basis_element) + translation
            translation_basis = [0, 0, 0]
            for index in range(3):
                while rot_trans[index] < -1e-3:
                    rot_trans[index] += 1
                    translation_basis[index] -= 1
                while rot_trans[index] > 1 - 1e-3:
                    rot_trans[index] -= 1
                    translation_basis[index] += 1
            rotation_translation.append(translation_basis)
            found = False
            for basis_index, basis_element_comp in enumerate(basis):
                if (abs(rot_trans - basis_element_comp) < 1e-3).all():
                    assert not found
                    basis_shifts[i, j] = basis_index
                    found = True
            assert found
        rotation_translations.append(np.array(rotation_translation))

    symmetries['translations'] = rotation_translations
    symmetries['basis_shifts'] = basis_shifts
    return symmetries


def get_atoms_from_labeling(labeling, cell, hnf, subelements, basis):
    '''
    Get ASE Atoms object from labeling, HNF matrix and parent lattice.

    Parameters
    ---------
    labeling : tuple
        Permutation of index of elements.
    A : ndarray
        Parent lattice (basis vectors listed column wise).
    hnf : ndarray
        HNF matrix defining the supercell.
    subelements : list of str
        List of elements, e.g. ['Au', 'Ag']

    Returns
    -------
    ASE Atoms
        Atoms object corresponding to the given labeling.
    '''
    symbols = []
    positions = []
    count = 0
    for i in range(hnf.H[0, 0]):
        coord = i * hnf.H[1, 0]
        offset10 = coord // hnf.H[0, 0] + coord % hnf.H[0, 0]
        coord = i * hnf.H[2, 0]
        offset20 = coord // hnf.H[0, 0] + coord % hnf.H[0, 0]
        for j in range(hnf.H[1, 1]):
            coord = j * hnf.H[2, 1]
            offset21 = coord // hnf.H[1, 1] + coord % hnf.H[1, 1]
            for k in range(hnf.H[2, 2]):
                for basis_vector in basis:

                    positions.append(i * cell[0] +
                                     (j + offset10) * cell[1] +
                                     (k + offset20 + offset21) * cell[2] +
                                     np.dot(cell, basis_vector))
                    symbols.append(subelements[labeling[count]])
                    count += 1
    return Atoms(symbols, positions, cell=np.dot(hnf.H.T, cell),
                 pbc=(True, True, True))


def enumerate_structures(atoms, sizes, subelements):
    '''
    Generate enumerated structures, i.e. all inequivalent structures up to a
    certain size

    Paramters
    ---------
    atoms : ASE Atoms
        Primitive structure from which derivative superstructures should be
        generated.
    size : list of ints
        Maximum number of atoms in the returned structures.
    subelements : list of str
        Elements to decorate the structure, e.g. ['Au', 'Ag']

    Yields
    ------

    '''
    nelements = len(subelements)
    assert nelements > 1

    nsites = len(atoms)
    basis = atoms.get_scaled_positions()
    symmetries = get_symmetry_operations(atoms, basis)

    # Loop over each cell size
    for ncells in sizes:
        count = 0

        hnfs = get_reduced_hnfs(ncells, symmetries)
        snfs = get_unique_snfs(hnfs)

        for snf in snfs:
            labelkeys = get_labelkeys(snf, nelements, nsites)

            for hnf in snf.hnfs:
                for labeling in yield_unique_labelings(labelkeys, snf, hnf,
                                                       nsites, nelements):
                    yield get_atoms_from_labeling(labeling, atoms.cell, hnf,
                                                  subelements, basis)
                    count += 1
