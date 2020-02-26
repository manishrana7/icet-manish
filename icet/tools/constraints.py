import numpy as np
from scipy.linalg import null_space
from typing import List

class Constraints:
    """ Class for handling linear constraints with right hand side equal to zero.

    Parameters
    ----------
    n_params
        number of parameters in model

    Example
    -------
    The following example demonstrates fitting of a cluster expansion under the
    constraint that ECI 2 and ECI 4 should be equal.

        >>> from icet.tools import Constraints
        >>> import numpy as np
        
        >>> # Set up random sensing matrix and target "energies"
        >>> n_params = 10
        >>> A = np.random((10, n_params))
        >>> y = np.random(10)

        >>> # Define constraints
        >>> c = Constraints(n_params=n_params)
        >>> c.add_constraint([2, 4])

        >>> # Do the actual fit and finally extract parameters
        >>> A_constrained = c.transform(A)
        >>> opt = Optimizer((A_constrained, y), fit_method='ridge')
        >>> opt.train()
        >>> parameters = c.inverse_transform(opt.parameters)

    """

    def __init__(self, n_params: int):
        self.M = np.empty((0, n_params))
        self.constraint_vectors = np.eye(n_params)

    def transform(self, A: np.ndarray) -> np.ndarray:
        """ Transform array to constrained parameter space

        Parameters
        ----------
        A
            array to be transformed
         """
        return A.dot(self.constraint_vectors)

    def inverse_transform(self, A: np.ndarray) -> np.ndarray:
        """ Inverse transform array from constrained parameter space
        to unconstrained space

        Parameters
        ----------
        A
            array to be inversed transformed
        """
        return self.constraint_vectors.dot(A)

    def add_constraint(self, parameter_indices: List[int]) -> None:
        """ Add a constraint and resolve for the constraint space

        Parameters
        ----------
        parameter_indices
            Indices of parameters the sum of which should be zero
        """
        M = np.zeros((1, self.M.shape[1]))
        M[0, parameter_indices] = 1
        self.M = np.vstack((self.M, M))
        self.constraint_vectors = null_space(self.M)
