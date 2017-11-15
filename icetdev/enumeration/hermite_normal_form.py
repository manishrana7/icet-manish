import numpy as np


def yield_hermite_normal_forms(det):
    '''
    Yield all Hermite Normal Form matrices with determinant det.

    Paramters
    ---------
    det : int
        Target determinant of HNFs

    Yields
    ------
    ndarray
        3x3 HNF matrix
    '''
    for a in range(1, det + 1):
        if det % a == 0:
            for c in range(1, det // a + 1):
                if det // a % c == 0:
                    f = det // (a * c)
                    for b in range(0, c):
                        for d in range(0, f):
                            for e in range(0, f):
                                hnf = [[a, 0, 0],
                                       [b, c, 0],
                                       [d, e, f]]
                                yield np.array(hnf)


def get_reduced_hermite_normal_forms(N, A, rotations_inv):
    '''
    For a fixed determinant N (i.e., a number of atoms N), yield all
    Hermite Normal Forms (HNF) that are inequivalent under symmetry
    operations of the parent lattice.'

    Paramters
    ---------
    N : int
        Determinant (or, equivalently, the number of atoms) of the HNF.
    A : ndarray
        Parent lattice (basis vectors listed column-wise).
    rotations_inv : list of ndarrays
        Inverse rotation matrices of the parent lattice.

    Returns
    ------
    list of ndarrays
        Symmetrically inequivalent HNFs with determinant N.
    '''
    hnfs = []
    Bs_previous = []
    for hnf in yield_hermite_normal_forms(N):

        B_i = np.dot(A, hnf)
        B_i_inv = np.linalg.inv(B_i)

        # Throw away HNF:s that yield equivalent supercells
        duplicate = False
        for R_inv in rotations_inv:
            BR = np.dot(B_i_inv, R_inv)
            for B_j in Bs_previous:
                check = np.dot(BR, B_j)
                check = check - np.round(check)
                if (abs(check) < 1e-3).all():
                    duplicate = True
                    break
            if duplicate:
                break
        if duplicate:
            continue

        # If it's not a duplicate, save the hnf
        # and the supercell so that it can be compared to
        hnfs.append(hnf)
        Bs_previous.append(B_i)
    return hnfs
