
"""
BaseOptimizer serves as base for all optimizers.
"""

import numpy as np
from .fit_methods import available_fit_methods


class BaseOptimizer:
    """
    BaseOptimizer class.

    Serves as base class for all Optimizers solving `Ax = y`.

    Parameters
    ----------
    fit_data : tuple of NumPy (N, M) array and NumPy (N) array
        the first element of the tuple represents the fit matrix `A`
        whereas the second element represents the vector of target
        values `y`; here `N` (=rows of `A`, elements of `y`) equals the number
        of target values and `M` (=columns of `A`) equals the number of
        parameters
    fit_method : string
        method to be used for training; possible choice are
        "least-squares", "lasso", "elasticnet", "bayesian-ridge", "ardr"
    standardize : bool
        whether or not to standardize the fit matrix before fitting
    seed : int
        seed for pseudo random number generator
    """

    def __init__(self, fit_data, fit_method, standardize=True, seed=42):
        """
        Attributes
        ----------
        _A : NumPy (N, M) array
            fit matrix
        _y : NumPy (N) array
            target values
        """

        if fit_method not in available_fit_methods:
            raise ValueError('Unknown fit_method: {}'.format(fit_method))

        if fit_data[0].shape[0] != fit_data[1].shape[0]:
            raise ValueError('Invalid fit data; shapes of fit matrix'
                             ' and target vector do not match')

        if len(fit_data[0].shape) != 2:
            raise ValueError('Invalid fit matrix; must have two dimensions')

        self._A, self._y = fit_data
        self._n_rows = self._A.shape[0]
        self._n_cols = self._A.shape[1]
        self._fit_method = fit_method
        self._standarize = standardize
        self._seed = seed
        self._fit_results = {'parameters': None}

    def compute_rmse(self, A, y):
        """
        Compute the root mean square error using the `A`, `y`, and the
        vector of fitted parameters `x` corresponding to `||Ax-y||_2`.

        Parameters
        ----------
        A : NumPy (N, M) array
            fit matrix where `N` (=rows of `A`, elements of `y`) equals the
            number of target values and `M` (=columns of `A`) equals the number
            of parameters (=elements of `x`)
        y : NumPy (N) array
            vector of target values

        Returns
        -------
        float
            root mean squared error
        """
        y_predicted = self.predict(A)
        delta_y = y_predicted - y
        rmse = np.sqrt(np.mean(delta_y**2))
        return rmse

    def predict(self, A):
        """
        Predict data given an input matrix `A`, i.e., `Ax`, where `x` is
        the vector of the fitted parameters.

        Parameters
        ----------
        A : NumPy (N, M) array or NumPy (M, )
            fit matrix where `N` (=rows of `A`, elements of `y`) equals the
            number of target values and `M` (=columns of `A`) equals the number
            of parameters

        Returns
        -------
        NumPy (N) array, or float if single row is inputed
            vector of predicted values
        """
        return np.dot(A, self.parameters)

    def get_contributions(self, A):
        """
        Compute the average contribution to the predicted values from each
        element of the parameter vector.

        Parameters
        ----------
        A : NumPy (N, M) array
            fit matrix where `N` (=rows of `A`, elements of `y`) equals the
            number of target values and `M` (=columns of `A`) equals the number
            of parameters

        Returns
        -------
        NumPy (N, M) array
            average contribution for each row of `A` from each parameter
        """
        return np.mean(np.abs(np.multiply(A, self.parameters)), axis=0)

    @property
    def summary(self):
        """ dict : Comprehensive information about the optimizer. """
        info = dict()
        info['fit_method'] = self.fit_method
        info['standardize'] = self.standardize
        info['number_of_target_values'] = self.number_of_target_values
        info['number_of_parameters'] = self.number_of_parameters
        info['number_of_nonzero_parameters'] = \
            self.number_of_nonzero_parameters
        return {**info, **self._fit_results}

    def __str__(self):
        width = 54
        s = []
        s.append(' {} '.format(self.__class__.__name__).center(width, '='))
        for key, value in self.summary.items():
            if isinstance(value, (str, int)):
                s.append('{:30} : {}'.format(key, value))
            elif isinstance(value, (float)):
                s.append('{:30} : {:.7g}'.format(key, value))
        s.append(''.center(width, '='))
        return '\n'.join(s)

    def __repr__(self):
        return 'BaseOptimizer((A, y), {}, {}'.format(
            self.fit_method, self.seed)

    @property
    def fit_method(self):
        """ string : fit method. """
        return self._fit_method

    @property
    def parameters(self):
        """ NumPy array : copy of parameter vector. """
        if self._fit_results['parameters'] is None:
            return None
        else:
            return self._fit_results['parameters'].copy()

    @property
    def number_of_nonzero_parameters(self):
        """ int : number of non-zero parameters """
        if self.parameters is None:
            return None
        else:
            return np.count_nonzero(self.parameters)

    @property
    def number_of_target_values(self):
        """ int : number of target values (=rows in `A` matrix). """
        return self._n_rows

    @property
    def number_of_parameters(self):
        """ int : number of parameters (=columns in `A` matrix). """
        return self._n_cols

    @property
    def standardize(self):
        """ bool : whether or not to standardize the fit matrix before
                   fitting.
        """
        return self._standarize

    @property
    def seed(self):
        """ int : seed used to initialize pseudo random number generator."""
        return self._seed
