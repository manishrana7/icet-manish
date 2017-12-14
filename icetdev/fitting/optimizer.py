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

    Attributes
    ----------
    train_scatter_data : ScatterData object (namedtuple)
        contains target and predicted value for each row in the training set
    test_scatter_data : ScatterData object (namedtuple)
        contains target and predicted value for each row in the test set
    '''

    def __init__(self, fit_data, fit_method='least-squares',
                 train_fraction=0.75, test_fraction=None, train_rows=None,
                 test_rows=None, seed=42, **kwargs):
        ''' Initialize Optimizer.

        Either specify train_fraction/test_fraction or train_rows/test_rows,
        these can not be used together.

        Note that Optimizer.train_fraction will be slightly different than
        the specified value.

        Todo
        ----
        * document parameters
        '''

        super().__init__(fit_data, fit_method, seed)

        self._seed = seed
        self._kwargs = kwargs

        # setup train and test rows
        self._setup_rows(train_fraction, test_fraction, train_rows, test_rows)

        # will be populate once running train
        self._rmse_train = None
        self._rmse_test = None
        self.train_scatter_data = None
        self.test_scatter_data = None

    def train(self):
        ''' Carry out training '''
        print('{s:-^{n}}'.format(s='Training', n=45))

        # select training data
        A_train = self._A[self.train_rows, :]
        y_train = self._y[self.train_rows]

        # perform training
        print('Fit Method {}, N_params {}'.format(self.fit_method, self.Ncols))
        print('Train size {}, Test size {} '
              .format(self.train_size, self.test_size))

        self._fit_results = self._optimizer_function(A_train, y_train,
                                                     **self._kwargs)
        self._rmse_train = self.compute_rmse(A_train, y_train)
        self.train_scatter_data = ScatterData(y_train, self.predict(A_train))
        print('Train RMSE  {:5.5f}'.format(self.rmse_train))

        # perform validation
        if self.test_rows is not None:
            A_test = self._A[self.test_rows, :]
            y_test = self._y[self.test_rows]
            self._rmse_test = self.compute_rmse(A_test, y_test)
            self.test_scatter_data = ScatterData(y_test, self.predict(A_test))
            print('Test  RMSE  {:5.5f}'.format(self.rmse_test))
        else:
            self._rmse_test = None
            self.test_scatter_data = None
            print('Test  RMSE  NaN')
        print('{s:-^{n}}'.format(s='Done', n=45))

    def _setup_rows(self, train_fraction, test_fraction, train_rows,
                    test_rows):
        ''' Setup train and test rows depending on which arguments are
        specified and which are None.

        If `train_rows` and `test_rows` are `None` then `fractions` are used.
        '''

        if train_rows is None and test_rows is None:  # get rows from fraction
            train_rows, test_rows = self._get_rows_via_fractions(
                train_fraction, test_fraction)
        else:  # get rows from specified rows
            train_rows, test_rows = self._get_rows_via_rows(
                train_rows, test_rows)

        if len(train_rows) == 0:
            raise ValueError('No training rows was selected from fit_data')

        self._train_rows = train_rows
        self._test_rows = test_rows

    def _get_rows_via_fractions(self, train_fraction, test_fraction):
        ''' Gets row via fractions. '''

        # Handle special cases
        if test_fraction is None and train_fraction is None:
            raise ValueError('Both train fraction and test fraction are None')
        elif train_fraction is None and abs(test_fraction - 1.0) < 1e-10:
            raise ValueError('train rows is empty for these fractions')
        elif test_fraction is None and abs(train_fraction - 1.0) < 1e-10:
            train_rows = np.arange(self._A.shape[0])
            test_rows = None
            return train_rows, test_rows

        # split
        train_rows, test_rows = train_test_split(
            np.arange(self._A.shape[0]), train_size=train_fraction,
            test_size=test_fraction, random_state=self.seed)
        if len(test_rows) == 0:
            test_rows = None
        if len(train_rows) == 0:
            raise ValueError('train rows is empty, too small train_fraction')

        return train_rows, test_rows

    def _get_rows_via_rows(self, train_rows, test_rows):
        ''' Gets row via specified rows '''
        if train_rows is None and test_rows is None:
            raise ValueError('Both train_rows and test_rows are None')
        elif test_rows is None:
            test_rows = [i for i in range(self.Nrows) if i not in train_rows]
        elif train_rows is None:
            train_rows = [i for i in range(self.Nrows) if i not in test_rows]
        return np.array(train_rows), np.array(test_rows)

    def get_info(self):
        ''' Get a dict containing all information about optimizer '''
        info = BaseOptimizer.get_info(self)
        info['rmse_train'] = self.rmse_train
        info['rmse_test'] = self.rmse_test
        info['train_rows'] = self.train_rows
        info['test_rows'] = self.test_rows
        info['train_size'] = self.train_size
        info['test_size'] = self.test_size
        return info

    def __str__(self):
        s = []
        s.append(super().__str__())
        s.append('Train size: {}'.format(self.train_size))
        s.append('Test  size: {}'.format(self.test_size))
        return '\n'.join(s)

    @property
    def rmse_train(self):
        ''' float : root mean squared error for training set '''
        return self._rmse_train

    @property
    def rmse_test(self):
        ''' float : root mean squared error for test set '''
        return self._rmse_test

    @property
    def train_rows(self):
        ''' list : indices for the rows included in training '''
        return self._train_rows

    @property
    def test_rows(self):
        ''' list : indices for the rows included in testing '''
        return self._test_rows

    @property
    def train_fraction(self):
        ''' float : fraction of the total data included in training '''
        return len(self._train_rows) / self.Nrows

    @property
    def test_fraction(self):
        ''' float : fraction of the total data included in testing '''
        return len(self._test_rows) / self.Nrows

    @property
    def train_size(self):
        ''' int : number of rows included in training '''
        return len(self.train_rows)

    @property
    def test_size(self):
        ''' int : number of rows included in testing '''
        if self.test_rows is None:
            return 0
        return len(self.test_rows)
