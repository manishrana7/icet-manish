import itertools
from spglib import get_symmetry
#from icetdev.enumeration.smith_normal_form import get_smith_normal_form
import numpy as np
from ase import Atoms
from icetdev.enumeration.hermite_normal_form \
    import yield_hermite_normal_forms, get_reduced_hermite_normal_forms


def dehash_labelkey(labelkey, nbr_of_atoms, nbr_of_elements):
    labeling = ()
    for i in range(nbr_of_atoms):
        labeling += (labelkey % nbr_of_elements,)
        labelkey = labelkey // nbr_of_elements
    return labeling


def hash_labeling(labeling, nbr_of_atoms, nbr_of_elements):
    labelkey = 0
    for i in range(nbr_of_atoms):
        labelkey += labeling[i] * nbr_of_elements**i
    return labelkey


def get_snfs_and_dangerous_rotations(hnfs, rotations, basis_shifts):
    '''
    For a list of HNF matrices, calculate their corresponding SNFs (Smith
    Normal Form matrices). Also calculate the rotations of the parent lattice
    that leaves the supercell corresponding to each HNF unchanged (these
    rotation may still change the labeling, so we need them later on).

    Paramters
    ---------
    hnfs : list of ndarrays
        HNF matrices, typically calculated with get_reduced_hermite_normal_forms()
    rotations_inv : list of ndarrays
        Inverse rotation matrices of the parent lattice.

    Returns
    -------
    list of tuples
        SNF matrices. Each tuple contains the diagonal element of each
        occuring SNF matrix.
    list of ndarrays
        Matrices multiplying SNF from the left.
    list of lists
        Map from SNF to HNF. Each sublist correponds to an SNF, and the
        elements of the SNF are indices to HNF in the input list of hnfs.
    list of list ndarrays
        Rotations that leave the supercell unchanged. Each sublist corresponds
        to an HNF matrix and contains all the rotations that leaves the cell
        defined by that HNF unchanged.
    '''
    snfs = []
    
    for hnf_index, hnf in enumerate(hnfs):
        # Get SNF for HNF matrix
        #snf_matrix, L, _ = get_smith_normal_form(hnf)
        #snf = tuple([snf_matrix[i, i] for i in range(3)])

        # Check whether the snf is new or already encountered
        snf = hnf.snf
        snf_is_new = True
        for snf_index, snf_comp in enumerate(snfs):
            if snf_comp.S == snf.S:
                snf_is_new = False
                snf_comp.add_hnf(hnf)
                break
        if snf_is_new:
            snf.G = get_group_representation(snf.S)
            snf.add_hnf(hnf)
            snfs.append(snf)

        # Save transformations (based on rotations) that turns the
        # supercell into an equivalent supercell
        # Should be moved to HNF
        hnf_rots_single = []
        for R, basis_shift in zip(rotations, basis_shifts):
            check = np.dot(np.dot(np.linalg.inv(hnf.H), R), hnf.H)
            check = check - np.round(check)
            if (abs(check) < 1e-3).all():
                LRL = np.dot(hnf.snf.L, np.dot(R, np.linalg.inv(hnf.snf.L)))

                # Should be an integer matrix
                assert (abs(LRL - np.round(LRL)) < 1e-3).all()
                LRL = np.round(LRL).astype(np.int64)
            
                hnf.add_transformation([LRL, basis_shift])
    return snfs


def yield_translation_permutations(original, snf, nbr_of_sites, nbr_of_elements, include_self=False):
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

    # Compute size of each block within which translations occur
    size_1 = nbr_of_sites * snf.blocks[0]
    size_2 = nbr_of_sites * snf.blocks[1]
    size_3 = nbr_of_sites * snf.blocks[2]
    
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
                    block_1 = original[size_1 * i:size_1 * (i + 1)]
                    for j in new_positions_2:
                        block_2 = block_1[size_2 * j:size_2 * j + size_2]
                        for k in new_positions_3:
                            for label in block_2[size_3 * k:size_3 * k + size_3]: 
                                labelkey += label*nbr_of_elements**count
                                count += 1
                yield labelkey
    

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
    return np.array(G)


def get_labelings(snf, nbr_of_elements, nbr_of_sites):
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
    N = snf.N
    nbr_of_atoms = N*nbr_of_sites
    labelkey_tracker = [False] * nbr_of_elements**nbr_of_atoms
    labelings = []
    for labeling in itertools.product(list(range(nbr_of_elements)), repeat=nbr_of_atoms):
        labelkey = hash_labeling(labeling, nbr_of_atoms, nbr_of_elements)
        
        unique = True
        #if is_superperiodic(labeling, N):
        #    continue
        for labelkey_symm in yield_translation_permutations(labeling, snf, nbr_of_sites, nbr_of_elements, include_self=False):
            if labelkey == labelkey_symm:
                unique = False
                break

            # Check with previous labelings,
            # if labeling can be translated into a previously
            # added labeling, then it is not unique
            if labelkey_tracker[labelkey_symm]:
                unique = False
                break
            if not unique:
                break
        if unique:
            labelings.append(labeling)
            labelkey_tracker[labelkey] = True
    return labelings


def permute_labeling(labeling, snf, nbr_of_sites, nbr_of_elements, rotation, basis_shift):
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


    Gp = np.dot(snf.G, rotation.T)
   
    labelkey = 0
    count = 0
    for member in Gp:
        cell_index = 0
        for i in range(3):
            cell_index += (member[i] % snf.S[i]) * snf.blocks[i]
        index = cell_index
        for basis_index in basis_shift:
            index = nbr_of_sites*cell_index + basis_index
            labelkey += labeling[index]*nbr_of_elements**count
            count += 1
    return labelkey


def yield_unique_labelings(labelings, snf, transformations, nbr_of_sites, nbr_of_elements):
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
    nbr_of_atoms = snf.N * nbr_of_sites
    labelkey_tracker = [False]*nbr_of_elements**nbr_of_atoms
    labelings_yielded = []
    for labeling in labelings:
        labelkey = hash_labeling(labeling, nbr_of_atoms, nbr_of_elements)
        
        # Check whether labeling is just a rotated version of a previous
        # labeling
        unique = True
        for transformation in transformations:
            labelkey_rot = permute_labeling(labeling, snf, nbr_of_sites, nbr_of_elements, 
                                            rotation=transformation[0],
                                            basis_shift=transformation[1])
    
            # Commonly, the transformation leaves the labeling
            # unchanged, so check that first as a special case
            # (yields a quite significant speedup)
            if labelkey_rot == labelkey:
                continue
            labeling_rot = dehash_labelkey(labelkey_rot, nbr_of_atoms, nbr_of_elements)
            
            for labelkey_rot_trans in \
                    yield_translation_permutations(labeling_rot, snf, nbr_of_sites, nbr_of_elements,
                                              include_self=True):
                                
                if labelkey_tracker[labelkey_rot_trans]:
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
    
    symmetries = get_symmetry(atoms)
    rotations = symmetries['rotations']
    translations = symmetries['translations']

    # Calculate how atoms within the primitive cell are
    # shifted upon operation with rotation matrix
    basis_shifts = np.zeros((len(rotations), len(basis)), dtype='int64')
    for i in range(len(rotations)):
        rotation = rotations[i]
        translation = translations[i]
        for j, basis_element in enumerate(basis):
            Rd = np.dot(rotation, basis_element) + translation
            for index in range(3):
                while Rd[index] < -1e-3:
                    Rd[index] += 1
                while Rd[index] > 1-1e-3:
                    Rd[index] -= 1
            found = False
            for basis_index, basis_element_comp in enumerate(basis):
                if (abs(Rd - basis_element_comp) < 1e-3).all():
                    assert not found
                    basis_shifts[i, j] = basis_index
                    found = True
            assert found

    return rotations, basis_shifts


def get_atoms_from_labeling(labeling, A, hnf, subelements, basis):
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
                for basis_vector in basis:
                    positions.append(i * A[:, 0] + (j + offset_10) * A[:, 1] +
                                    (k + offset_20 + offset_21) * A[:, 2] +
                                    np.dot(A, basis_vector))
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

    nbr_of_sites = len(atoms)
    basis = atoms.get_scaled_positions()
    rotations, basis_shifts = get_symmetry_operations(atoms, basis)
    A = atoms.cell.T

    #exit(0)    

    # Loop over each cell size
    for N in sizes:
        count = 0

        nbr_of_atoms = N*nbr_of_sites
        hnfs = get_reduced_hermite_normal_forms(N, rotations)

        snfs = get_snfs_and_dangerous_rotations(hnfs, rotations, basis_shifts)
        for snf_index, snf in enumerate(snfs):
            labelings = get_labelings(snf, nbr_of_elements, nbr_of_sites)
            
            for hnf in snf.hnfs:
                for labeling in yield_unique_labelings(labelings,
                                                       snf,
                                                       hnf.transformations,
                                                       nbr_of_sites, nbr_of_elements):
                    yield get_atoms_from_labeling(labeling, A, hnf.H,
                                                  subelements, basis)
                    count += 1
                    if False:
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


if __name__ == '__main__':
    from ase.build import bulk
    from ase.build import fcc111
    import time
    import argparse

    parser = argparse.ArgumentParser(description='Enumerate structures')
    parser.add_argument('size', help='Maximum number of atoms in enumerated'
                        ' structures')
    args = parser.parse_args()

    size = range(1, int(args.size)+1)
    atoms = bulk('Au', a=4.0, crystalstructure='fcc')
    #atoms = fcc111('Au', a=4.0, size=(1,1,5), vacuum=7.0)
    
    #print(atoms)
    cell = atoms.cell
    cell[0] = 1.0*cell[0]
    atoms.set_cell(cell)
    
    subelements = ['Au', 'Ag']
    #cutoffs = [11, 8, 5]

    #cs = ClusterSpace(atoms, cutoffs, subelements)

    start = time.time()

    #from ase.db import connect
    # db = connect(
    #    '~/repos/structure-enumeration/databases/fcc_binary_12atoms.db')
    #db = connect('test.db')

    #from ase.visualize import view
    count_structures = 0
    cvs_my = []
    cvs_en = []
    a = 0
    for structure in enumerate_structures(atoms, size, subelements):
        # if len(structure) < size:
        #    continue
        # print(np.linalg.det(structure.get_cell()))
        count_structures += 1
        # print(count_structures)
        #cv = cs.get_clustervector(structure)

        # try:
        #    cv = cs.get_clustervector(structure)
        # except:
        # view(structure)
        #enumlib_structure = db.get_atoms(count_structures)
        # view(enumlib_structure)
        # exit(0)

        # print(count_structures)
        # db.write(structure)
        # cvs_my.append(cv)

        '''
        for c in cv:
            print('{:8.4f}'.format(c), end='')
        print()
        '''

        #enumlib_structure = db.get_atoms(count_structures)
        #cv = cs.get_clustervector(enumlib_structure)
        # cvs_en.append(cv)

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
