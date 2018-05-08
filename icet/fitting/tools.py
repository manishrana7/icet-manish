"""
Collection of tools for managing and analysisng linear models

Todo
-----
Consider what functionality we actually want here

"""

import numpy as np
from collections import namedtuple


ScatterData = namedtuple('ScatterData', ['target', 'predicted'])


def compute_correlation_matrix(A):
    """
    Computes correlation matrix for rows in input matrix.

    Naive implementation.

    Parameters
    ----------
    A : NumPy array
        fit matrix

    Returns
    -------
    NumPy array
        correlation matrix
    """
    N = A.shape[0]
    C = np.zeros((N, N))
    for i in range(N):
        for j in range(i+1, N):
            norm = np.linalg.norm(A[i, :]) * np.linalg.norm(A[j, :])
            c_ij = np.dot(A[i, :], A[j, :]) / norm
            C[i, j] = c_ij
            C[j, i] = c_ij
    return C
