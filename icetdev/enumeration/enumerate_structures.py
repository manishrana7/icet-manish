import itertools
from spglib import get_symmetry
from icetdev.enumeration.smith_normal_form import get_smith_normal_form
import numpy as np
from ase import Atoms
from icetdev.enumeration.hermite_normal_form \
    import yield_hermite_normal_forms, get_reduced_hermite_normal_forms


def get_snfs_and_dangerous_rotations(hnfs, A, Ainv, rotations_inv):
    '''
    For a list of HNF matrices, calculate their corresponding SNFs (Smith
    Normal Form matrices). Also calculate the rotations of the parent lattice
    that leaves the supercell corresponding to each HNF unchanged (these
    rotation may still change the labeling, so we need them later on).

    Paramters
    ---------
    hnfs : list of ndarrays
        HNF matrices, typically calculated with get_reduced_hermite_normal_forms()
    A : ndarray
        Parent lattice (basis vectors listed column wise).
    Ainv : ndarray
        Inverse of parent lattice (i.e., inverse of A).
    rotations_inv : list of ndarrays
        Inverse rotation matrices of the parent lattice.

    Returns
    -------
    list of tuples
        SNF matrices. Each tuple contains the diagonal element of each
        occuring SNF matrix.
    list of lists
        Map from SNF to HNF. Each sublist correponds to an SNF, and the
        elements of the SNF are indices to HNF in the input list of hnfs.
    list of list ndarrays
        Rotations that leave the supercell unchanged. Each sublist corresponds
        to an HNF matrix and contains all the rotations that leaves the cell
        defined by that HNF unchanged.
    '''
    snfs = []
    snf_to_hnf_map = []
    hnf_rots = []

    for hnf_index, hnf in enumerate(hnfs):
        # Get SNF for HNF matrix
        snf_matrix, L, _ = get_smith_normal_form(hnf)
        snf = tuple([snf_matrix[i, i] for i in range(3)])

        # Check whether the snf is new or already encountered
        snf_is_new = True
        for snf_index, snf_comp in enumerate(snfs):
            if snf_comp == snf:
                snf_is_new = False
                snf_to_hnf_map[snf_index].append(hnf_index)
                break
        if snf_is_new:
            snfs.append(snf)
            snf_to_hnf_map.append([hnf_index])

        # Save transformations (based on rotations) that turns the
        # supercell into an equivalent supercell
        LA = np.dot(L, Ainv)
        hnf_rots_single = []
        B_i = np.dot(A, hnf)
        B_i_inv = np.linalg.inv(B_i)
        for R_inv in rotations_inv:
            check = np.dot(np.dot(B_i_inv, R_inv), B_i)
            check = check - np.round(check)
            if (abs(check) < 1e-3).all():
                LARLA = np.dot(LA, np.linalg.inv(np.dot(LA, R_inv)))
                assert (abs(LARLA - np.round(LARLA)) < 1e-3).all()
                LARLA = np.round(LARLA).astype(np.int64)
                hnf_rots_single.append(LARLA)
        hnf_rots.append(hnf_rots_single)

    return snfs, snf_to_hnf_map, hnf_rots


def yield_symmetry_equivalent(original, snf, include_self=False):
    '''
    Yield labelings that are equivalent to original labeling
    under translations as dictated by snf.

    Paramters
    ---------
    original : tuple
        labelings to be translated
    snf : tuple
        diagonal of matrix on Smith Normal Form
    include_self : bool
        inlcude original labeling or not

    Yields
    ------
    tuple
        labeling that is equivalent to original
    '''

    N = snf[0] * snf[1] * snf[2]
    assert(len(original) == N)

    # Compute size of each block within which translations occur
    size_1 = N // snf[0]
    size_2 = size_1 // snf[1]
    size_3 = size_2 // snf[2]

    # Loop over each possible translation
    for trans_1 in range(snf[0]):
        for trans_2 in range(snf[1]):
            for trans_3 in range(snf[2]):
                if not include_self and trans_1 + trans_2 + trans_3 == 0:
                    continue
                # Calculate the positions to which the atoms will be moved
                new_positions_1 = [(i + trans_1) % snf[0]
                                   for i in range(snf[0])]
                new_positions_2 = [(i + trans_2) % snf[1]
                                   for i in range(snf[1])]
                new_positions_3 = [(i + trans_3) % snf[2]
                                   for i in range(snf[2])]

                # Construct new labeling by building up tuple step by step
                labeling_translated = ()
                for i in new_positions_1:
                    block_1 = original[size_1 * i:size_1 * i + size_1]
                    for j in new_positions_2:
                        block_2 = block_1[size_2 * j:size_2 * j + size_2]
                        for k in new_positions_3:
                            labeling_translated += block_2[
                                size_3 * k:size_3 * k + size_3]
                yield labeling_translated


def is_superperiodic(labeling, N):
    '''
    Check whether labeling is superperiodic, i.e. ABCABC or AABBAABB. These
    should be generated with a smaller cell, i.e. we should already have ABC
    and AABB among our labelings.

    Parameters
    ----------
    labeling : tuple
        labeling to be checked

    Returns
    -------
    bool
        True if labeling is superperiodic, False otherwise
    '''
    for translate in range(1, N):
        labeling_translated = tuple(
            [labeling[(i + translate) % N] for i in range(N)])
        if labeling_translated == labeling:
            return True
    return False


def get_group_representation(snf):
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
    G = []
    for i in range(snf[0]):
        for j in range(snf[1]):
            for k in range(snf[2]):
                G.append([i, j, k])
    G = np.array(G).T
    return G


def get_labelings(snf, nbr_of_elements):
    '''
    Get all labelings corresponding to a Smith Normal Form. Superperiodic
    labelings as well as labelings that are equivalent for this particular SNF
    will not be included. However, labelings that are equivalent by rotations
    that leave the cell (but not the labeling) unchanged will still be
    included, since these have to be removed for each HNF.

    Paramters
    ---------
    snf : tuple
        The three diagonal elements of an SNF matrix

    Returns
    -------
    list of tuples
        Symmetrically inequivalent labelings
    '''
    N = snf[0] * snf[1] * snf[2]
    labelings = []
    for labeling in itertools.product(list(range(nbr_of_elements)), repeat=N):
        unique = True
        if is_superperiodic(labeling, N):
            continue
        for labeling_symm in yield_symmetry_equivalent(labeling, snf):
            if labeling_symm == labeling:
                unique = False
                break
            for labeling_comp in labelings:
                if labeling_symm == labeling_comp:
                    unique = False
                    break
            if not unique:
                break
        if unique:
            labelings.append(labeling)
    return labelings


def get_rotated_labeling(labeling, snf, Gp):
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
    N = snf[0] * snf[1] * snf[2]
    assert N == len(labeling)
    assert N == len(Gp[0])

    # Gp should only contain integers
    assert((abs(Gp - np.round(Gp)) < 1e-3).all())
    Gp = np.round(Gp).astype(np.int64).T

    blocks = [N // snf[0]]
    for i in range(1, 3):
        blocks.append(blocks[-1] // snf[i])

    labeling_rotated = ()
    for member in Gp:
        index = 0
        for i in range(3):
            index += (member[i] % snf[i]) * blocks[i]
        labeling_rotated += (labeling[index],)
    return labeling_rotated


def yield_unique_labelings(labelings, snf, hnf_rots, G):
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
        labeling
    G : ndarray
        Group representation corresponding to snf.

    Yields
    ------
    tuple
        Labeling, each and every one unique.
    '''
    labelings_yielded = []
    for labeling in labelings:

        # Check whether labeling is just a rotated version of a previous
        # labeling
        unique = True
        for transformation in hnf_rots:
            labeling_rot = get_rotated_labeling(
                labeling, snf, np.dot(transformation, G))

            # Commonly, the transformation leaves the labeling
            # unchanged, so check that first as a special case
            # (yields a quite significant speedup)
            if labeling_rot == labeling:
                continue
            for labeling_rot_symm in \
                    yield_symmetry_equivalent(labeling_rot, snf,
                                              include_self=True):
                for labeling_previous in labelings_yielded:
                    if labeling_previous == labeling_rot_symm:
                        unique = False
                        break
                # if not unique: # This is not necessarily wrong
                #    break       # but gives no speedup
            if not unique:
                break
        if unique:
            # Then we have finally found a unique structure
            # defined by an HNF matrix and a labeling
            yield labeling

            # Still needs to save the labeling so perhaps this should not be a
            # generator...
            labelings_yielded.append(labeling)


def get_inverse_rotation_matrices(atoms):
    '''
    Use spglib to calculate the symmetry operations of atoms and return their
    inverse matrices.

    Parameters
    ----------
    atoms : ASE Atoms
        Structure for which the symmetry operations are sought.

    Returns
    -------
    list of ndarrays
        Inverse of matrices for rotation operatios of atoms.
    '''

    # The spglib documentations says atoms.get_cell().T but that seems to be
    # wrong...?
    A = atoms.get_cell()
    lattice = np.array(A, dtype='double', order='C')
    positions = np.array(atoms.get_scaled_positions(), dtype='double',
                         order='C')
    numbers = np.array(atoms.get_atomic_numbers(), dtype='intc')
    cell = (lattice, positions, numbers)
    symmetries = get_symmetry(cell, symprec=1e-5)
    rotations = symmetries['rotations']

    # Save rotations in Cartesian coordinates.
    # We actually only need the inverse matrices.
    rotations_inv = []
    for r in rotations:
        R = np.dot(np.dot(A.T, r), np.linalg.inv(A.T))
        rotations_inv.append(np.linalg.inv(R))
    return rotations_inv


def get_atoms_from_labeling(labeling, A, hnf, subelements):
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
    for i in range(hnf[0, 0]):
        a = i * hnf[1, 0]
        offset_10 = a // hnf[0, 0] + a % hnf[0, 0]
        a = i * hnf[2, 0]
        offset_20 = a // hnf[0, 0] + a % hnf[0, 0]
        for j in range(hnf[1, 1]):
            a = j * hnf[2, 1]
            offset_21 = a // hnf[1, 1] + a % hnf[1, 1]
            for k in range(hnf[2, 2]):
                positions.append(i * A[:, 0] + (j + offset_10) * A[:, 1] + \
                    (k + offset_20 + offset_21) * A[:, 2])
                symbols.append(subelements[labeling[count]])
                count += 1
    return Atoms(symbols, positions, cell=np.dot(A, hnf).T, pbc=(True, True, True)) 


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
    nbr_of_elements = len(subelements)
    assert nbr_of_elements > 1

    rotations_inv = get_inverse_rotation_matrices(atoms)
    A = atoms.cell.T
        
    # Loop over each cell size
    for N in sizes:
        hnfs = get_reduced_hermite_normal_forms(N, A, rotations_inv)
        snfs, snf_to_hnf_map, hnf_rots = get_snfs_and_dangerous_rotations(
            hnfs, A, np.linalg.inv(A), rotations_inv)

        for snf_index, snf in enumerate(snfs):

            labelings = get_labelings(snf, nbr_of_elements)
            G = get_group_representation(snf)

            for hnf_index in snf_to_hnf_map[snf_index]:
                hnf = hnfs[hnf_index]
                for labeling in yield_unique_labelings(labelings,
                                                       snf,
                                                       hnf_rots[hnf_index], G):
                    yield get_atoms_from_labeling(labeling, A, hnf,
                                                  subelements)
                    # yield labeling
                    if False:
                        '''
                        print('{:7}  '.format(count_structures), end='')
                        for i in snf:
                            print('{:3}'.format(i), end='')
                        print('  ', end='')
                        hnf = hnfs[hnf_index]
                        for i in range(3):
                            for j in range(i + 1):
                                print('{:3}'.format(hnf[i, j]), end='')
                        print('    ', end='')
                        for i in labeling:
                            print(i, end='')
                        print()
                        '''


if __name__ == '__main__':
    from ase.build import bulk
    import time
    import argparse
    from icetdev.clusterspace import ClusterSpace, Structure

    parser = argparse.ArgumentParser(description='Enumerate structures')
    parser.add_argument('size', help='Maximum number of atoms in enumerated'
                        ' structures')
    args = parser.parse_args()

    size = int(args.size)

    atoms = bulk('Au', a=4.0, crystalstructure='fcc')
    cell = atoms.cell
    cell[0] = cell[0]
    atoms.set_cell(cell)

    subelements = ['Au', 'Ag']
    cutoffs = [11, 8, 5]

    #cs = ClusterSpace(atoms, cutoffs, subelements)

    start = time.time()

    from ase.db import connect
    #db = connect(
    #    '~/repos/structure-enumeration/databases/fcc_binary_12atoms.db')
    db = connect('test.db')
    
    #from ase.visualize import view
    count_structures = 0
    cvs_my = []
    cvs_en = []
    a = 0
    for structure in get_enumerated_structures(atoms, size, subelements):

        # if len(structure) < size:
        #    continue
        #print(np.linalg.det(structure.get_cell()))
        count_structures += 1
        #print(count_structures)
        #cv = cs.get_clustervector(structure)

        #try:
        #    cv = cs.get_clustervector(structure)
        #except:
            #view(structure)
            #enumlib_structure = db.get_atoms(count_structures)
            #view(enumlib_structure)
            #exit(0)

        #print(count_structures)
        #db.write(structure)
        #cvs_my.append(cv)
        
        '''
        for c in cv:
            print('{:8.4f}'.format(c), end='')
        print()
        '''

        #enumlib_structure = db.get_atoms(count_structures)
        #cv = cs.get_clustervector(enumlib_structure)
        #cvs_en.append(cv)

        '''
        for c in cv:
            print('{:8.4f}'.format(c), end='')
        print()
        print()
        '''

    end = time.time()
    print('Time consumed: {}'.format(end - start))
    print('Number of structures: {}'.format(count_structures))

    '''
    used = []
    for cv_my in cvs_my:
        match_found = False
        for j, cv_en in enumerate(cvs_en):
            if np.allclose(cv_my, cv_en):
                match_found = True
                if j in used:
                    print('Duplicate found')
                else:
                    used.append(j)
        if not match_found:
            print('MATCH LACKING')
    '''
