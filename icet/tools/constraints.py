import numpy as np
from scipy.linalg import null_space
from .structure_enumeration import enumerate_structures


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
        >>> from icet.fitting import Optimizer
        >>> import numpy as np

        >>> # Set up random sensing matrix and target "energies"
        >>> n_params = 10
        >>> A = np.random.random((10, n_params))
        >>> y = np.random.random(10)

        >>> # Define constraints
        >>> c = Constraints(n_params=n_params)
        >>> M = np.zeros((1, n_params))
        >>> M[0, [2, 4]] = 1
        >>> c.add_constraint(M)

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

    def add_constraint(self, M: np.ndarray) -> None:
        """ Add a constraint matrix and resolve for the constraint space

        Parameters
        ----------
        M
            Constraint matrix with each constraint as a row. Can (but need not be)
            cluster vectors.
        """
        M = np.array(M)
        self.M = np.vstack((self.M, M))
        self.constraint_vectors = null_space(self.M)


def get_mixing_energy_constraints(cluster_space) -> Constraints:
    """
    A cluster expansion of *mixing energy* should ideally predict zero energy
    for concentration 0 and 1. This function constructs a `Constraint` object
    that enforces that condition during fitting.

    Parameters
    ----------
    cluster_space : ClusterSpace
        Cluster space corresponding to cluster expansion for which constraints
        should be imposed
    """
    M = []
    for structure in enumerate_structures(structure=cluster_space.primitive_structure,
                                          sizes=[1],
                                          chemical_symbols=cluster_space.chemical_symbols):
        M.append(cluster_space.get_cluster_vector(structure))
    c = Constraints(n_params=len(cluster_space))
    c.add_constraint(M)
    return c
