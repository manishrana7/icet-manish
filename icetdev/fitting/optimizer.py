'''
Optimizer
'''
import numpy as np
from sklearn.model_selection import train_test_split
from icetdev.fitting.tools import ScatterData
from icetdev.fitting.base_optimizer import BaseOptimizer


class Optimizer(BaseOptimizer):
    '''
    Optimizer for single `Ax = y` fit.


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
    train_fraction : float
        fraction of input data (=rows) to be used for training
    test_fraction : float
        fraction of input data (=rows) to be used for testing
    training_set : tuple/list of ints
        indices of rows of `A`/`y` to be used for training
    testing_set : tuple/list of ints
        indices of rows of `A`/`y` to be used for testing
    seed : int
        seed for pseudo random number generator

    Attributes
    ----------
    training_scatter_data : ScatterData object (namedtuple)
        target and predicted value for each row in the training set
    testing_scatter_data : ScatterData object (namedtuple)
        target and predicted value for each row in the testing set
    '''

    def __init__(self, fit_data, fit_method='least-squares',
                 train_fraction=0.75, test_fraction=None,
                 training_set=None, testing_set=None, seed=42, **kwargs):

        super().__init__(fit_data, fit_method, seed)

        self._seed = seed
        self._kwargs = kwargs

        # setup training and testing sets
        self._setup_rows(train_fraction, test_fraction,
                         training_set, testing_set)

        # will be populate once running train
        self._rmse_training_set = None
        self._rmse_testing_set = None
        self.training_scatter_data = None
        self.testing_scatter_data = None

    def train(self):
        ''' Carry out training. '''
        print('{s:-^{n}}'.format(s='Training', n=45))

        # select training data
        A_train = self._A[self.training_set, :]
        y_train = self._y[self.training_set]

        # perform training
        print('Fit Method {}, N_params {}'.format(self.fit_method,
                                                  self.number_of_parameters))
        print('Train size {}, Test size {} '
              .format(self.training_set_size, self.testing_set_size))

        self._fit_results = self._optimizer_function(A_train, y_train,
                                                     **self._kwargs)
        self._rmse_training_set = self.compute_rmse(A_train, y_train)
        self.training_scatter_data = ScatterData(y_train,
                                                 self.predict(A_train))
        print('Train RMSE  {:5.5f}'.format(self.rmse_training_set))

        # perform validation
        if self.testing_set is not None:
            A_test = self._A[self.testing_set, :]
            y_test = self._y[self.testing_set]
            self._rmse_testing_set = self.compute_rmse(A_test, y_test)
            self.testing_scatter_data = ScatterData(y_test,
                                                    self.predict(A_test))
            print('Test  RMSE  {:5.5f}'.format(self.rmse_testing_set))
        else:
            self._rmse_testing_set = None
            self.testing_scatter_data = None
            print('Test  RMSE  NaN')
        print('{s:-^{n}}'.format(s='Done', n=45))

    def _setup_rows(self, train_fraction, test_fraction, training_set,
                    testing_set):
        ''' Setup train and test rows depending on which arguments are
        specified and which are None.

        If `training_set` and `testing_set` are `None` then `fractions` is
        used.
        '''

        if training_set is None and testing_set is None:
            # get rows from fractions
            training_set, testing_set = \
                self._get_rows_via_fractions(train_fraction, test_fraction)
        else:  # get rows from specified rows
            training_set, testing_set = \
                self._get_rows_from_indices(training_set, testing_set)

        if len(training_set) == 0:
            raise ValueError('No training rows was selected from fit_data')

        self._training_set = training_set
        self._testing_set = testing_set

    def _get_rows_via_fractions(self, train_fraction, test_fraction):
        ''' Gets row via fractions. '''

        # Handle special cases
        if test_fraction is None and train_fraction is None:
            raise ValueError('Both train fraction and test fraction are None')
        elif train_fraction is None and abs(test_fraction - 1.0) < 1e-10:
            raise ValueError('train rows is empty for these fractions')
        elif test_fraction is None and abs(train_fraction - 1.0) < 1e-10:
            training_set = np.arange(self._Nrows)
            testing_set = None
            return training_set, testing_set

        # split
        training_set, testing_set = \
            train_test_split(np.arange(self._Nrows),
                             train_size=train_fraction,
                             test_size=test_fraction,
                             random_state=self.seed)
        if len(testing_set) == 0:
            testing_set = None
        if len(training_set) == 0:
            raise ValueError('train rows is empty, too small train_fraction')

        return training_set, testing_set

    def _get_rows_from_indices(self, training_set, testing_set):
        ''' Gets row via indices '''
        if training_set is None and testing_set is None:
            raise ValueError('Both training_set and testing_set are None')
        elif testing_set is None:
            testing_set = [i for i in range(self._Nrows)
                           if i not in training_set]
        elif training_set is None:
            training_set = [i for i in range(self._Nrows)
                            if i not in testing_set]
        return np.array(training_set), np.array(testing_set)

    def get_info(self):
        ''' Get comprehensive information concerning the optimization process.

        Returns
        -------
        dict
        '''
        info = BaseOptimizer.get_info(self)
        info['rmse training set'] = self.rmse_training_set
        info['rmse testing set'] = self.rmse_testing_set
        info['training set size'] = self.training_set_size
        info['testing set size'] = self.testing_set_size
        info['training set'] = self.training_set
        info['testing set'] = self.testing_set
        return info

    @property
    def rmse_training_set(self):
        ''' float : root mean squared error for training set '''
        return self._rmse_training_set

    @property
    def rmse_testing_set(self):
        ''' float : root mean squared error for testing set '''
        return self._rmse_testing_set

    @property
    def training_set(self):
        ''' list : indices of the rows included in the training set '''
        return self._training_set

    @property
    def testing_set(self):
        ''' list : indices of the rows included in the testing set '''
        return self._testing_set

    @property
    def training_set_size(self):
        ''' int : number of rows included in training set '''
        if self.training_set is None:
            return 0
        return len(self.training_set)

    @property
    def testing_set_size(self):
        ''' int : number of rows included in testing set '''
        if self.testing_set is None:
            return 0
        return len(self.testing_set)
