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


if __name__ == '__main__':
    from hiphive.fitting import Optimizer
    np.random.seed(42)
    n, m = 20, 10

    A = np.random.random((n, m))
    y = np.random.random(n)

    # constraint sum eci[inds] = 0
    inds1 = [1, 3, 4, 5]
    inds2 = [2, 6, 7, 8]
    M = np.zeros((2, m))
    M[0, inds1] = 1
    M[1, inds2] = 1

    c = Constraints(m)
    c.add_constraint(M)

    Ac = c.transform(A)
    opt = Optimizer((Ac, y), fit_method='ridge')
    opt.train()

    parameters = c.inverse_transform(opt.parameters)
    print('constraints 1, ', parameters[inds1].sum())
    print('constraints 2, ', parameters[inds2].sum())
