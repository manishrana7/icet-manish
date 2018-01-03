'''
Ensemble Optimizer

https://en.wikipedia.org/wiki/Bootstrap_aggregating
http://scikit-learn.org/stable/modules/generated/sklearn.ensemble.BaggingRegressor.html  # NOQA
'''

import numpy as np
from .base_optimizer import BaseOptimizer
from .optimizer import Optimizer


class EnsembleOptimizer(BaseOptimizer):
    '''Ensemble optimizer that carries out a series of single optimization runs
    using the :class:`Optimizer` class and then provides access to various
    ensemble averaged quantities including e.g., errors and parameters.

    Parameters
    ----------
    fit_data : tuple of (N, M) numpy.ndarray and (N) numpy.ndarray
        the first element of the tuple represents the fit matrix `A`
        whereas the second element represents the vector of target
        values `y`; here `N` (=rows of `A`, elements of `y`) equals the number
        of target values and `M` (=columns of `A`) equals the number of
        parameters
    fit_method : str
        method to be used for training; possible choice are
        "least-squares", "lasso", "bayesian-ridge", "ardr"
    n_splits : int
        number of fits in the ensemble
    train_fraction : float
        fraction of input data (=rows) to be used for training
    test_fraction : float
        fraction of input data (=rows) to be used for testing
    bootstrap : boolean
        if True sampling will be carried out with replacement
    seed : int
        seed for pseudo random number generator
    '''

    def __init__(self, fit_data, fit_method='least-squares', n_splits=50,
                 train_fraction=0.7, test_fraction=None,
                 bootstrap=True, seed=42, **kwargs):

        BaseOptimizer.__init__(self, fit_data, fit_method, seed)
        self._n_splits = n_splits
        self._training_set_size = int(np.round(train_fraction * self._Nrows))
        if test_fraction is not None:
            self._testing_set_size = int(np.round(test_fraction * self._Nrows))
            assert (self._training_set_size +
                    self._testing_set_size <= self._Nrows), \
                'Size of training and testing sets exceeds available data'
        else:
            self._testing_set_size = self._Nrows - self._training_set_size

        self._bootstrap = bootstrap
        self._kwargs = kwargs

    def train(self):
        ''' Carry out ensemble training. '''
        self._run_ensemble()
        self._construct_final_model()

    def _run_ensemble(self):
        ''' Carry out training. '''

        np.random.seed(self.seed)

        parameters_list = []
        rmse_training_set_list, rmse_testing_set_list = [], []
        training_set_list, testing_set_list = [], []
        for _ in range(self.n_splits):
            # construct training and testing sets
            rows = np.random.choice(
                self._Nrows, self.training_set_size + self.testing_set_size,
                replace=self.bootstrap)
            train_rows = np.random.choice(rows, self.training_set_size,
                                          replace=self.bootstrap)
            test_rows = np.setdiff1d(rows, range(self._Nrows))

            # train
            opt = Optimizer((self._A, self._y), self.fit_method,
                            train_rows=train_rows, test_rows=test_rows,
                            **self._kwargs)
            opt.train()

            # collect results
            parameters_list.append(opt.parameters)
            rmse_training_set_list.append(opt.rmse_training_set)
            rmse_testing_set_list.append(opt.rmse_testing_set)
            training_set_list.append(train_rows)
            testing_set_list.append(test_rows)

        self._parameters_set = np.array(parameters_list)
        self._training_set_list = training_set_list
        self._testing_set_list = testing_set_list
        self._average_rmse_training_set = np.average(rmse_training_set_list)
        self._average_rmse_testing_set = np.average(rmse_testing_set_list)

    def _construct_final_model(self):
        ''' Construct final model. '''
        self._fit_results['parameters'] = np.mean(self.parameters_set, axis=0)

    def get_errors(self):
        ''' Get the errors for each fit and each target value.

        Returns
        -------
        (N,M) numpy.ndarray
            matrix of fit errors where `N` is the number of target values and
            `M` is the number of fits (i.e., the size of the ensemble)
        '''
        error_matrix = np.zeros((self._Nrows, self.n_splits))
        for i, parameters in enumerate(self.parameters_set):
            error_matrix[:, i] = np.dot(self._A, parameters) - self._y
        return error_matrix

    def get_parameters_avg(self):
        '''Get average of each parameter over the ensemble.

        Returns
        -------
        numpy.ndarray
            vector of average values
        '''
        return np.average(self.parameters_set, axis=0)

    def get_parameters_std(self):
        '''Get standard deviation of each parameter over the ensemble.

        Returns
        -------
        numpy.ndarray
            vector of standard deviations
        '''
        return np.std(self.parameters_set, axis=0)

    @property
    def parameter_vectors(self):
        ''' list : all parameter vectors in the ensemble '''
        return self._parameters_set

    @property
    def ensemble_size(self):
        ''' int : number of rounds of training '''
        return self._n_splits

    @property
    def rmse_training_set(self):
        ''' float : ensemble average of root mean squared error over training
        set '''
        return self._average_rmse_training_set

    @property
    def rmse_testing_set(self):
        ''' float : ensemble average of root mean squared error over testing
        set '''
        return self._average_rmse_testing_set

    @property
    def training_set_size(self):
        ''' int : number of rows included in training sets '''
        return self._training_set_size

    @property
    def testing_set_size(self):
        ''' int : number of rows included in testing sets '''
        return self._testing_set_size

    @property
    def fractional_training_set_size(self):
        ''' float : fraction of input data used for training; this value can
                    differ slightly from the value set during initialization
                    due to rounding '''
        return float(self.training_set_size) / self._Nrows

    @property
    def bootstrap(self):
        ''' bool : True if sampling is carried out with replacement '''
        return self._bootstrap
