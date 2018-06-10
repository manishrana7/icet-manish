"""
Optimizer with cross validation score
"""

import numpy as np
from sklearn.model_selection import KFold, ShuffleSplit
from .base_optimizer import BaseOptimizer
from .optimizer import Optimizer
from .fit_methods import fit
from .tools import ScatterData


validation_methods = {
    'k-fold': KFold,
    'shuffle-split': ShuffleSplit,
}


class CrossValidationEstimator(BaseOptimizer):
    """
    Optimizer with cross validation.

    This optimizer first computes a cross-validation score and finally
    generates a model using the full set of input data.

    Warning
    -------
    Repeatedly setting up a CrossValidationEstimator and training
    *without* changing the seed for the random number generator will yield
    identical or correlated results, to avoid this please specify a different
    seed when setting up multiple CrossValidationEstimator instances.

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
    validation_method : string
        method to use for cross-validation; possible choices are
        "shuffle-split", "k-fold"
    number_of_splits : int
        number of times the fit data set will be split for the cross-validation
    seed : int
        seed for pseudo random number generator

    Attributes
    ----------
    train_scatter_data : ScatterData object (namedtuple)
        contains target and predicted values from each individual
        traininig set in the cross-validation split
    validation_scatter_data : ScatterData object (namedtuple)
        contains target and predicted values from each individual
        validation set in the cross-validation split
    """

    def __init__(self, fit_data, fit_method='least-squares', standardize=True,
                 validation_method='k-fold', number_of_splits=10,
                 seed=42, **kwargs):

        super().__init__(fit_data, fit_method, standardize, seed)

        if validation_method not in validation_methods.keys():
            msg = ['Validation method not available']
            msg += ['Please choose one of the following:']
            for key in validation_methods:
                msg += [' * ' + key]
            raise ValueError('\n'.join(msg))
        self._validation_method = validation_method
        self._number_of_splits = number_of_splits

        self._set_kwargs(kwargs)

        # data set splitting object
        self._splitter = validation_methods[validation_method](
            n_splits=self.number_of_splits, random_state=seed,
            **self._split_kwargs)

        self.train_scatter_data = None
        self.validation_scatter_data = None

        self._rmse_train_splits = None
        self._rmse_valid_splits = None
        self._rmse_train_final = None

    def train(self):
        """ Construct the final model using all input data available. """
        self._fit_results = fit(self._A, self._y, self.fit_method,
                                self.standardize, **self._fit_kwargs)
        self._rmse_train_final = self.compute_rmse(self._A, self._y)

    def validate(self):
        """ Run validation. """
        train_target, train_predicted = [], []
        valid_target, valid_predicted = [], []
        rmse_train_splits, rmse_valid_splits = [], []
        for train_set, test_set in self._splitter.split(self._A):
            opt = Optimizer((self._A, self._y), self.fit_method,
                            train_set=train_set,
                            test_set=test_set,
                            **self._fit_kwargs)
            opt.train()

            rmse_train_splits.append(opt.rmse_train)
            rmse_valid_splits.append(opt.rmse_test)
            train_target.extend(opt.train_scatter_data.target)
            train_predicted.extend(opt.train_scatter_data.predicted)
            valid_target.extend(opt.test_scatter_data.target)
            valid_predicted.extend(opt.test_scatter_data.predicted)

        self._rmse_train_splits = np.array(rmse_train_splits)
        self._rmse_valid_splits = np.array(rmse_valid_splits)
        self.train_scatter_data = ScatterData(
            target=np.array(train_target), predicted=np.array(train_predicted))
        self.validation_scatter_data = ScatterData(
            target=np.array(valid_target), predicted=np.array(valid_predicted))

    def _set_kwargs(self, kwargs):
        """ Sets up fit_kwargs and split_kwargs

        The different split methods need different keywords.
        """
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

    @property
    def summary(self):
        """ dict : Comprehensive information about the optimizer. """
        info = super().summary

        # Add class specific data
        info['validation_method'] = self.validation_method
        info['number_of_splits'] = self.number_of_splits
        info['rmse_train_final'] = self.rmse_train_final
        info['rmse_train'] = self.rmse_train
        info['rmse_train_splits'] = self.rmse_train_splits
        info['rmse_validation'] = self.rmse_validation
        info['rmse_validation_splits'] = self.rmse_validation_splits
        info['train_scatter_data'] = self.train_scatter_data
        info['validation_scatter_data'] = self.validation_scatter_data

        # add kwargs used for fitting and splitting
        info = {**info, **self._fit_kwargs, **self._split_kwargs}
        return info

    def __repr__(self):
        kwargs = dict()
        kwargs['fit_method'] = self.fit_method
        kwargs['validation_method'] = self.validation_method
        kwargs['number_of_splits'] = self.number_of_splits
        kwargs['seed'] = self.seed
        kwargs = {**kwargs, **self._fit_kwargs, **self._split_kwargs}
        return 'CrossValidationEstimator((A, y), {})'.format(
            ', '.join('{}={}'.format(*kwarg) for kwarg in kwargs.items()))

    @property
    def validation_method(self):
        """ string : validation method name. """
        return self._validation_method

    @property
    def number_of_splits(self):
        """ string : number of splits (folds) used for cross-validation. """
        return self._number_of_splits

    @property
    def rmse_train_final(self):
        """
        float : root mean squared error when using the full set of input data.
        """
        return self._rmse_train_final

    @property
    def rmse_train(self):
        """ float : average root mean squared training error obtained during
                    cross-validation. """
        if self._rmse_train_splits is None:
            return None
        return np.sqrt(np.mean(self._rmse_train_splits**2))

    @property
    def rmse_train_splits(self):
        """ list : root mean squared training errors obtained during
                   cross-validation. """
        return self._rmse_train_splits

    @property
    def rmse_validation(self):
        """ float : average root mean squared cross-validation error. """
        if self._rmse_valid_splits is None:
            return None
        return np.sqrt(np.mean(self._rmse_valid_splits**2))

    @property
    def rmse_validation_splits(self):
        """ list : root mean squared validation errors obtained during
                   cross-validation. """
        return self._rmse_valid_splits
