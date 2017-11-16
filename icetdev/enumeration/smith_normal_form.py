import numpy as np

class SmithNormalForm:
    def __init__(self, H):
        '''
        Get Smith Normal Form for 3x3 matrix.
        
        Parameters
        ----------
        H : ndarray
            3x3 matrix
        
        Returns
        -------
        ndarray
            Smith Normal form matrix S for H
        ndarray
            Matrix L that multiplies H from the left.
        ndarray
            Matrix R that multuplies R from the right such that L*H*R=S
        '''
        A = H.copy()
        L = np.eye(3, dtype=int)
        R = np.eye(3, dtype=int)
        while True:
            # Clear upper row and leftmost column in such
            # a way that greatest common denominator ends
            # up in A[0, 0], in a standard Smith Normal Form way
            # (Euclidean algorithm for finding greatest common divisor)
            while(sorted(A[0])[1] != 0 or sorted(A[:, 0])[1] != 0):
                A, R = _clear_row(A, R, 0)
                A, L = _clear_column(A, L, 0)

            # Do the same thing for lower 2x2 matrix 
            while(sorted(A[1, 1:])[0] != 0 or sorted(A[1:, 1])[0] != 0):
                A, R = _clear_row(A, R, 1)
                A, L = _clear_column(A, L, 1)

            # If last diagonal entry is negative,
            # make it positive
            if A[2, 2] < 0:
                A[2, 2] = -A[2, 2]
                L[2] = -L[2]
                
            # Check that the diagonal entry i,i divides
            # diagonal entry i+1, i+1. Otherwise,
            # add row i+1 to i and start over.
            if A[2, 2] % A[1, 1] != 0:
                A[1] = A[1] + A[2]
                L[1] = L[1] + L[2]
            elif A[1, 1] % A[0, 0] != 0:
                A[0] = A[0] + A[1]
                L[0] = L[0] + L[1]
            else:
                break
        assert (abs(np.dot(np.dot(L, H), R) - A) < 1e-3).all()

        self.S_matrix = A
        self.S = tuple([A[i, i] for i in range(3)])
        self.L = L
        self.R = R
        self.G = None
        self.hnfs = []


    def add_hnf(self, hnf):
        self.hnfs.append(hnf)



def _switch_rows(A, i, j):
    '''
    Switch rows in matrix.
    
    Parameters
    ---------
    A : ndarray
        Matrix in which rows will be swapped.
    i : int
        Index of row 1 to be swapped.
    j : int
        Index of row 2 to be swapped.

    Returns
    -------
    ndarray
        Matrix with swapped rows.
    '''
    row = A[j].copy()
    A[j] = A[i]
    A[i] = row
    return A


def _switch_columns(A, i, j):
    '''
    Switch columns in matrix.
    
    Parameters
    ---------
    A : ndarray
        Matrix in which columns will be swapped.
    i : int
        Index of column 1 to be swapped.
    j : int
        Index of column 2 to be swapped.

    Returns
    -------
    ndarray
        Matrix with swapped columns.
    '''
    col = A[:, j].copy()
    A[:, j] = A[:, i]
    A[:, i] = col
    return A


def _clear_row(A, R, i):
    '''
    Use column operations to make A[i, i] the greatest common
    denominator of the elements in the row and the other elements
    zero.

    Parameters
    ----------
    A : ndarray
        Matrix whose row is to be cleared.
    R : ndarray
        Matrix that should be subject to the same operations.
    i : int
        Index of row to be treated.
    
    Returns
    -------
    ndarray
        Treated matrix A.
    ndarray
        Matrix that has been subject to the same operations.
    '''
    for j in range(i, 3):
        if A[i, j] < 0:
            A[:, j] = -1*A[:, j]
            R[:, j] = -1*R[:, j]
    while(np.sort(A[i, i:])[1-i] > 0):
        max_index = np.argmax(A[i, i:]) + i
        min_index = np.argmin(A[i, i:]) + i
        if max_index == min_index:
            max_index += 1
        if A[i, min_index] == 0 and i == 0:
            if np.sort(A[i])[1] > 0:
                min_index += 1
                min_index = min_index % 3
                if min_index == max_index:
                    min_index += 1
                    min_index = min_index % 3
            if A[i, min_index] == A[i, max_index]:
                tmp = min_index
                min_index = max_index
                max_index = tmp
        A[:, max_index] = A[:, max_index] - A[:, min_index]
        R[:, max_index] = R[:, max_index] - R[:, min_index]
    max_index = np.argmax(A[i])
    A = _switch_columns(A, i, max_index)
    R = _switch_columns(R, i, max_index)
    return A, R

def _clear_column(A, L, j):
    '''
    Use row operations to make A[i, i] the greatest common
    denominator of the elements in the column and the other elements
    zero.

    Parameters
    ----------
    A : ndarray
        Matrix whose column is to be cleared.
    R : ndarray
        Matrix that should be subject to the same operations.
    i : int
        Index of column to be treated.
    
    Returns
    -------
    ndarray
        Treated matrix A.
    ndarray
        Matrix that has been subject to the same operations.
    '''
    for i in range(j, 3):
        if A[i, j] < 0:
            A[i] = -1*A[i]
            L[i] = -1*L[i]
    while(np.sort(A[j:, j])[1-j] > 0):
        max_index = np.argmax(A[j:, j]) + j
        min_index = np.argmin(A[j:, j]) + j
        if max_index == min_index:
            max_index += 1
        if A[min_index, j] == 0 and j == 0:
            if np.sort(A[:, j])[1] > 0:
                min_index += 1
                min_index = min_index % 3
                if min_index == max_index:
                    min_index += 1
                    min_index = min_index % 3
            if A[min_index, j] == A[max_index, j]:
                tmp = min_index
                min_index = max_index
                max_index = tmp
        A[max_index] = A[max_index] - A[min_index]
        L[max_index] = L[max_index] - L[min_index]
    max_index = np.argmax(A[:, j])
    A = _switch_rows(A, j, max_index)
    L = _switch_rows(L, j, max_index)
    return A, L

if __name__=='__main__':
    A = np.array([[2, 4, 4], [-6, 6, 12], [10, -4, -16]])
    S, L, R = get_smith_normal_form(A.copy())
    print(S)
    print(L)
