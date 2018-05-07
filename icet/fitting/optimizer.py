"""
Optimizer
"""
import numpy as np
from sklearn.model_selection import train_test_split
from .tools import ScatterData
from .base_optimizer import BaseOptimizer


class Optimizer(BaseOptimizer):
    """
    Optimizer for single `Ax = y` fit.

    One has to specify either `train_size`/`test_size` or
    `train_set`/`test_set` If either `train_set` or `test_set` (or both)
    is specified the fractions will be ignored.

    Warning
    -------
    Repeatedly setting up a Optimizer and training
    *without* changing the seed for the random number generator will yield
    identical or correlated results, to avoid this please specify a different
    seed when setting up multiple Optimizer instances.

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
    train_size : float or int
        If float represents the fraction of `fit_data` (rows) to be used for
        training. If int, represents the absolute number of rows to be used for
        training.
    test_size : float or int
        If float represents the fraction of `fit_data` (rows) to be used for
        testing. If int, represents the absolute number of rows to be used for
        testing.
    train_set : tuple/list of ints
        indices of rows of `A`/`y` to be used for training
    test_set : tuple/list of ints
        indices of rows of `A`/`y` to be used for testing
    seed : int
        seed for pseudo random number generator

    Attributes
    ----------
    train_scatter_data : ScatterData object (namedtuple)
        target and predicted value for each row in the training set
    test_scatter_data : ScatterData object (namedtuple)
        target and predicted value for each row in the test set
    """

    def __init__(self, fit_data, fit_method='least-squares', train_size=0.75,
                 test_size=None, train_set=None, test_set=None, seed=42,
                 **kwargs):

        super().__init__(fit_data, fit_method, seed)

        self._kwargs = kwargs

        # setup train and test sets
        self._setup_rows(train_size, test_size,
                         train_set, test_set)

        # will be populate once running train
        self._rmse_train = None
        self._rmse_test = None
        self.train_scatter_data = None
        self.test_scatter_data = None

    def train(self):
        """ Carry out training. """

        # select training data
        A_train = self._A[self.train_set, :]
        y_train = self._y[self.train_set]

        # perform training
        self._fit_results = self._optimizer_function(A_train, y_train,
                                                     **self._kwargs)
        self._rmse_train = self.compute_rmse(A_train, y_train)
        self.train_scatter_data = ScatterData(y_train, self.predict(A_train))

        # perform testing
        if self.test_set is not None:
            A_test = self._A[self.test_set, :]
            y_test = self._y[self.test_set]
            self._rmse_test = self.compute_rmse(A_test, y_test)
            self.test_scatter_data = ScatterData(y_test, self.predict(A_test))
        else:
            self._rmse_test = None
            self.test_scatter_data = None

    def _setup_rows(self, train_size, test_size, train_set, test_set):
        """
        Set up train and test rows depending on which arguments are
        specified.

        If `train_set` and `test_set` are `None` then `train_size` and
        `test_size` are used.
        """

        if train_set is None and test_set is None:
            train_set, test_set = self._get_rows_via_sizes(
                train_size, test_size)
        else:
            train_set, test_set = self._get_rows_from_indices(
                train_set, test_set)

        if len(train_set) == 0:
            raise ValueError('No training rows was selected from fit_data')

        if test_set is not None:  # then check overlap between train and test
            if len(np.intersect1d(train_set, test_set)):
                raise ValueError('Overlap between train set and test set')
            if len(test_set) == 0:
                test_set = None

        self._train_set = train_set
        self._test_set = test_set

    def _get_rows_via_sizes(self, train_size, test_size):
        """ Gets train and test rows via sizes. """

        # Handle special cases
        if test_size is None and train_size is None:
            raise ValueError('Both train fraction and test fraction are None')
        elif train_size is None and abs(test_size - 1.0) < 1e-10:
            raise ValueError('train rows is empty for these fractions')

        elif test_size is None:
            if train_size == self._Nrows or abs(train_size-1.0) < 1e-10:
                train_set = np.arange(self._Nrows)
                test_set = None
                return train_set, test_set

        # split
        train_set, test_set = train_test_split(
            np.arange(self._Nrows), train_size=train_size, test_size=test_size,
            random_state=self.seed)

        return train_set, test_set

    def _get_rows_from_indices(self, train_set, test_set):
        """ Gets row via indices. """
        if train_set is None and test_set is None:
            raise ValueError('Both train and test set are None')
        elif test_set is None:
            test_set = [i for i in range(self._Nrows)
                        if i not in train_set]
        elif train_set is None:
            train_set = [i for i in range(self._Nrows)
                         if i not in test_set]
        return np.array(train_set), np.array(test_set)

    @property
    def summary(self):
        """ dict : Comprehensive information about the optimizer. """
        info = super().summary

        # Add class specific data
        info['rmse_train'] = self.rmse_train
        info['rmse_test'] = self.rmse_test
        info['train_size'] = self.train_size
        info['train_set'] = self.train_set
        info['test_size'] = self.test_size
        info['test_set'] = self.test_set
        info['train_scatter_data'] = self.train_scatter_data
        info['test_scatter_data'] = self.test_scatter_data

        # add kwargs used for fitting
        info = {**info, **self._kwargs}
        return info

    def __repr__(self):
        kwargs = dict()
        kwargs['fit_method'] = self.fit_method
        kwargs['traininig_size'] = self.train_size
        kwargs['test_size'] = self.test_size
        kwargs['train_set'] = self.train_set
        kwargs['test_set'] = self.test_set
        kwargs['seed'] = self.seed
        kwargs = {**kwargs, **self._kwargs}
        return 'Optimizer((A, y), {})'.format(
            ', '.join('{}={}'.format(*kwarg) for kwarg in kwargs.items()))

    @property
    def rmse_train(self):
        """ float : root mean squared error for training set. """
        return self._rmse_train

    @property
    def rmse_test(self):
        """ float : root mean squared error for test set. """
        return self._rmse_test

    @property
    def train_set(self):
        """ list : indices of the rows included in the training set. """
        return self._train_set

    @property
    def test_set(self):
        """ list : indices of the rows included in the test set. """
        return self._test_set

    @property
    def train_size(self):
        """ int : number of rows included in training set. """
        return len(self.train_set)

    @property
    def train_fraction(self):
        """ float : fraction of rows included in training set. """
        return self.train_size / self._Nrows

    @property
    def test_size(self):
        """ int : number of rows included in test set. """
        if self.test_set is None:
            return 0
        return len(self.test_set)

    @property
    def test_fraction(self):
        """ float : fraction of rows included in test set. """
        if self.test_set is None:
            return 0.0
        return self.test_size / self._Nrows
