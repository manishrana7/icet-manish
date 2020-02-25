import numpy as np
from scipy.linalg import null_space


class Constraints:
    """ Class for handling linear constraints with right hand side equal to zero.

    Parameters
    ----------
    n_params : int
        number of parameters in model
    """

    def __init__(self, n_params):
        self.M = np.empty((0, n_params))
        self.constraint_vectors = np.eye(n_params)

    def transform(self, A):
        """ Transform array to constrained parameter space

        Parameters
        ----------
        A : numpy array
            array to be transformed
         """
        return A.dot(self.constraint_vectors)

    def inverse_transform(self, A):
        """ Inverse transform array from constrained parameter space
        to unconstrained space

        Parameters
        ----------
        A : numpy array
            array to be inversed transformed
        """
        return self.constraint_vectors.dot(A)

    def add_constraint(self, M):
        """ Add a constraint matrix and resolve for the constraint space

        Parameters
        ----------
        M : numpy array
            Constraint matrix with each constraint as a row
        """
        self.M = np.vstack((self.M, M))
        self.constraint_vectors = null_space(self.M)
