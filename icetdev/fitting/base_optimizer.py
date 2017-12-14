'''
BaseOptimizer serves as base for all optimizers in HiPhive
'''

import numpy as np
from collections import OrderedDict
from icetdev.fitting.tools import compute_rmse
from icetdev.fitting.fit_methods import fit_least_squares, fit_lasso, \
    fit_bayesian_ridge, fit_ardr


fit_methods = OrderedDict([
    ('least-squares', fit_least_squares),
    ('lasso', fit_lasso),
    ('bayesian-ridge', fit_bayesian_ridge),
    ('ardr', fit_ardr),
    ])


class BaseOptimizer:
    ''' BaseOptimizer class.

    Serves as base for all Optimizers solving `Ax = y`

    Attributes
    ----------
    A : NumPy array
        fit matrix
    y : NumPy array
        target values
    optimizer_function : function
        optimizer function to be called when training
    '''

    def __init__(self, fit_data, fit_method, seed):

        if fit_method not in fit_methods.keys():
            raise ValueError('Fit method not available')

        if fit_data[0].shape[0] != fit_data[1].shape[0]:
            raise ValueError('Invalid fit data, shape did not match')

        self._A, self._y = fit_data
        self._fit_method = fit_method
        self._seed = seed
        self._optimizer_function = fit_methods[self.fit_method]
        self._fit_results = {'parameters': None}

    def compute_rmse(self, A, y):
        ''' Computes root mean square error for `A`, `y` with fitted parameters

        Parameters
        ----------
        A : NumPy array
            matrix
        y : NumPy array
            target values

        Returns
        -------
        float
            root mean squared error
        '''
        return compute_rmse(A, self.parameters, y)

    def predict(self, A):
        ''' Predict data given an input matrix `A`

        Parameters
        ----------
        A : NumPy array
            matrix

        Returns
        -------
        array
            prediction for each row in input matrix
        '''
        return np.dot(A, self.parameters)

    def contribution(self, A):
        ''' Compute the average contribution from each column (parameter)

        Parameters
        ----------
        A : NumPy array
            matrix

        Returns
        -------
        array
            average contribution from each column
        '''
        return np.mean(np.abs(np.multiply(A, self.parameters)), axis=0)

    def get_info(self):
        ''' Get a dict containing all information about optimizer '''
        info = dict()
        info['parameters'] = self.parameters
        info['fit_method'] = self.fit_method
        info['Nrows'] = self.Nrows
        info['Ncols'] = self.Ncols
        return info

    def __str__(self):
        s = []
        s.append('Fit method : {}'.format(self.fit_method))
        s.append('Nrows, Ncols : ({}, {})'.format(self.Nrows, self.Ncols))
        return '\n'.join(s)

    def __repr__(self):
        return(str(self))

    @property
    def fit_method(self):
        ''' str : fit method name '''
        return self._fit_method

    @property
    def parameters(self):
        ''' np.ndarray : copy of the parameters '''
        if self.fit_results['parameters'] is None:
            return None
        else:
            return self.fit_results['parameters'].copy()

    @property
    def Nrows(self):
        ''' int : number of rows in the A matrix '''
        return self._A.shape[0]

    @property
    def Ncols(self):
        ''' int : number of columns in the A matrix '''
        return self._A.shape[1]

    @property
    def seed(self):
        ''' int : random seed '''
        return self._seed

    @property
    def fit_results(self):
        ''' dict : dictionary containing results obtained from fitting '''
        return self._fit_results
