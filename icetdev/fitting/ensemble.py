'''
Ensemble Optimizer

https://en.wikipedia.org/wiki/Bootstrap_aggregating
http://scikit-learn.org/stable/modules/generated/sklearn.ensemble.BaggingRegressor.html
'''

import numpy as np
from icetdev.fitting.base_optimizer import BaseOptimizer
from icetdev.fitting import Optimizer


class EnsembleOptimizer(BaseOptimizer):
    ''' Ensemble Optimizer a.k.a Bagging

    This optimizers carries out an ensemble of fits using a fraction of the
    available fit data.

    Todo
    ----
    Potentially an ensemble estimate of the rmse_train and rmse_test could be
    added to this class. (ScatterData might be too much for this optimizer)
    '''

    def __init__(self, fit_data, fit_method='least-squares', n_splits=100,
                 train_fraction=0.7, bootstrap=True, seed=42, **kwargs):

        BaseOptimizer.__init__(self, fit_data, fit_method, seed)
        self._n_splits = n_splits
        self._train_size = int(np.round(train_fraction * self.Nrows))

        self._bootstrap = bootstrap
        self._kwargs = kwargs

    def train(self):
        ''' Run ensemble fitting '''
        self._run_ensemble()
        self._construct_final_model()

    def _run_ensemble(self):
        ''' run all fits '''

        np.random.seed(self.seed)

        parameters_list = []
        train_rows_list, test_rows_list = [], []
        for _ in range(self.n_splits):
            # select train and test rows
            train_rows = np.random.choice(
                np.arange(self.Nrows), self.train_size, replace=self.bootstrap)
            test_rows = np.array(list(set(np.arange(self.Nrows)) -
                                      set(train_rows)))

            # train
            opt = Optimizer((self._A, self._y), self.fit_method,
                            train_rows=train_rows, test_rows=test_rows,
                            **self._kwargs)
            opt.train()
            parameters_list.append(opt.parameters)
            train_rows_list.append(train_rows)
            test_rows_list.append(test_rows)

        self._parameters_set = np.array(parameters_list)
        self._train_rows_list = train_rows_list
        self._test_rows_list = test_rows_list

    def _construct_final_model(self):
        ''' Construct final model '''
        self._fit_results['parameters'] = np.mean(self.parameters_set, axis=0)

    def get_error_matrix(self):
        ''' Get error matrix

        error_matrix[i][j] is the error for row i for fit j
        '''
        error_matrix = np.zeros((self.Nrows, self.n_splits))
        for i, parameters in enumerate(self.parameters_set):
            error_matrix[:, i] = np.dot(self._A, parameters) - self._y
        return error_matrix

    def get_parameters_stds(self):
        ''' Get standard deviation of parameters '''
        return np.std(self.parameters_set, axis=0)

    @property
    def parameters_set(self):
        ''' list : list of all obtained parameters '''
        return self._parameters_set

    @property
    def n_splits(self):
        ''' int : Number of fits (size of ensemble) '''
        return self._n_splits

    @property
    def bootstrap(self):
        ''' bool : Boolean dictating if sampling with replacement or not '''
        return self._bootstrap

    @property
    def train_size(self):
        ''' int : Number of columns included in each fit '''
        return self._train_size

    @property
    def train_fraction(self):
        ''' float : Training fraction (note this might differ slightly from the
                    value set in the __init__) '''
        return self.train_size/self.Nrows
