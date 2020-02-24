import numpy as np
from scipy.linalg import null_space


class Constraints:

    def __init__(self, n_params):
        self._n_params = n_params
        self.M = np.zeros((1, n_params))
        self.rhs = [0]
        self.constraint_vectors = np.eye(self._n_params)

    def transform(self, A):
        return A.dot(self.constraint_vectors)

    def inverse_transform(self, A):
        return self.constraint_vectors.dot(A)

    def add_constraint(self, M, rhs):
        if not np.allclose(rhs, 0):
            raise ValueError('Cant handle rhs != 0')

        # add constraints
        self.M = np.vstack((self.M, M))
        self.rhs.extend(rhs)

        # find new constraint vectors
        self.constraint_vectors = null_space(self.M)


if __name__ == '__main__':
    from hiphive.fitting import Optimizer
    np.random.seed(42)
    n, m = 20, 5

    A = np.random.random((n, m))
    y = np.random.random(n)

    # constraint eci 2 and 4 to sum to 0
    M = np.zeros((1, m))
    M[0, 2], M[0, 4] = 1, 1

    c = Constraints(m)
    c.add_constraint(M, [0])

    Ac = c.transform(A)
    opt = Optimizer((Ac, y), fit_method='ridge')
    opt.train()

    parameters = c.inverse_transform(opt.parameters)
    print('eci-2 {:12.8f}'.format(parameters[2]))
    print('eci-4 {:12.8f}'.format(parameters[4]))
