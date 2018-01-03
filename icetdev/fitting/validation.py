'''
Optimizer with cross validation score
'''

import numpy as np
from sklearn.model_selection import KFold, ShuffleSplit
from .optimizer import Optimizer
from .base_optimizer import BaseOptimizer
from .tools import ScatterData


validation_methods = {
    'k-fold': KFold,
    'shuffle-split': ShuffleSplit,
}


class CrossValidationEstimator(BaseOptimizer):
    ''' Optimizer with cross validation.

    This optimizer carries out cross validation and computes the final model
    with the full data set.

    Attributes
    ----------
    train_scatter_data : ScatterData object (namedtuple)
        contains target and predicted value for the final training set
    validation_scatter_data : ScatterData object (namedtuple)
        contains target and predicted value for each row in each validation set
    '''

    def __init__(self, fit_data, fit_method='least-squares',
                 validation_method='k-fold', n_splits=10, seed=42, **kwargs):

        BaseOptimizer.__init__(self, fit_data, fit_method, seed)

        if validation_method not in validation_methods.keys():
            raise ValueError('Validation method not available')
        self._validation_method = validation_method
        self._n_splits = n_splits

        self._set_kwargs(kwargs)

        # data set splitting object
        self._splitter = validation_methods[validation_method](
            n_splits=self.n_splits, random_state=seed, **self._split_kwargs)

    def train(self):
        ''' Carry out cross-validation and final fit'''

        self._validate()
        self._construct_final_model()

    def _validate(self):
        ''' Run validation '''
        valid_target, valid_predicted = [], []
        rmse_train_split = []
        rmse_validation = []
        for train_rows, test_rows in self._splitter.split(self._A):
            opt = Optimizer((self._A, self._y), self.fit_method,
                            train_rows=train_rows, test_rows=test_rows,
                            **self._fit_kwargs)
            opt.train()

            rmse_train_split.append(opt.rmse_train)
            rmse_validation.append(opt.rmse_test)
            valid_target.extend(opt.test_scatter_data.target)
            valid_predicted.extend(opt.test_scatter_data.predicted)

        self._rmse_train_split = np.array(rmse_train_split)
        self._rmse_validation = np.array(rmse_validation)
        self.validation_scatter_data = ScatterData(
            target=np.array(valid_target), predicted=np.array(valid_predicted))

    def _construct_final_model(self):
        ''' Train the final model '''
        self._fit_results = self._optimizer_function(self._A, self._y,
                                                     **self._fit_kwargs)
        self._rmse_train_final = self.compute_rmse(self._A, self._y)
        self.train_scatter_data = ScatterData(self._y, self.predict(self._A))

    def _set_kwargs(self, kwargs):
        ''' Sets up fit_kwargs and split_kwargs

        The different split methods need different keywords.
        '''
        self._fit_kwargs = {}
        self._split_kwargs = {}

        if self.validation_method == 'k-fold':
            self._fit_kwargs = kwargs
        elif self.validation_method == 'shuffle-split':
            for key, val in kwargs.items():
                if key in ['test_size', 'train_size']:
                    self._split_kwargs[key] = val
                else:
                    self._fit_kwargs[key] = val

    def __str__(self):
        s = []
        s.append('Validation method: {}'.format(self.validation_method))
        s.append(super().__str__())
        s.append('N_splits: {}'.format(self.n_splits))
        return '\n'.join(s)

    @property
    def validation_method(self):
        ''' str : Validation method name '''
        return self._validation_method

    @property
    def n_splits(self):
        ''' str : number of splits (folds) for cross validation calculation '''
        return self._n_splits

    @property
    def rmse_train_final(self):
        ''' float : root mean squared error for the final training set '''
        return self._rmse_train_final

    @property
    def rmse_train(self):
        ''' list : root mean squared training errors obtained during
                   cross validation calculation '''
        return self._rmse_train_split

    @property
    def rmse_validation(self):
        ''' list : root mean squared validation errors obtained during
                   cross validation calculation '''
        return self._rmse_validation

    @property
    def validation_score(self):
        ''' float : average root mean squared validation error '''
        return np.mean(self._rmse_validation)
