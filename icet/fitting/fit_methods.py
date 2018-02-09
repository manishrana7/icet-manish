'''
scikit-learn is an excellent library for training linear models and provides a
large number of useful tools.

This module provides simplified interfaces for vaiours linear model regression
methods. These methods are set up in a way that work out of the box for typical
problems in cluster expansion and force constant potential construction. This
includes slight adjustments scitkit-learn default values.

If you would like more flexibility or extended functionality or ability to
fine-tune parameters that are not included in this interface, it is of course
possible to use scikit-learn directly.
More information about the sklearn linear models can be found at
http://scikit-learn.org/stable/modules/linear_model.html

'''

import numpy as np
from ..io.logging import logger
try:
    from sklearn.linear_model import (Lasso,
                                      LassoCV,
                                      BayesianRidge,
                                      ARDRegression)
    # arrangement of logger assignments is owed to pep8 requirements
    logger = logger.getChild('fit_methods')
except Exception:
    logger = logger.getChild('fit_methods')
    logger.warning('Failed to import scitkit-learn;'
                   ' several optimizers will fail')


def fit_least_squares(X, y):
    '''
    Return the least-squares solution `a` to the linear problem `Xa=y`.

    This function is a wrapper to the `linalg.lstsq` function in NumPy.

    Parameters
    -----------
    X : matrix / array
        fit matrix
    y : array
        target array

    Returns
    ----------
    results : dict
        dict containing parameters
    '''
    results = dict()
    results['parameters'] = np.linalg.lstsq(X, y, rcond=-1)[0]
    return results


def fit_lasso(X, y, alpha=None, fit_intercept=False, max_iter=5000, tol=1e-5,
              **kwargs):
    '''
    Return the solution `a` to the linear problem `Xa=y` obtained by using
    the LASSO method as implemented in scitkit-learn.

    LASSO optimizes the following problem::

        (1 / (2 * n_samples)) * ||y - Xw||^2_2 + alpha * ||w||_1

    If `alpha` is `None` this function will conduct a grid search for an
    optimal alpha value.

    Parameters
    -----------
    X : matrix / array
        fit matrix
    y : array
        target array
    alpha : float
        alpha value
    fit_intercept : bool
        center data or not, forwarded to sklearn

    Returns
    ----------
    results : dict
        dictionary containing parameters
    '''
    if alpha is None:
        return fit_lassoCV(X, y, fit_intercept=fit_intercept,
                           max_iter=max_iter, tol=tol, **kwargs)
    else:
        lasso = Lasso(alpha=alpha, fit_intercept=fit_intercept,
                      max_iter=max_iter, tol=tol**kwargs)
        lasso.fit(X, y)
        results = dict()
        results['parameters'] = lasso.coef_
        return results


def fit_lassoCV(X, y, alphas=None, fit_intercept=False, cv=10, max_iter=5000,
                tol=1e-5, **kwargs):
    '''
    Return the solution `a` to the linear problem `Xa=y` obtained by using
    the LassoCV method as implemented in scitkit-learn.

    Parameters
    -----------
    X : matrix / array
        fit matrix
    y : array
        target array
    alphas : list / array
        list of alpha values to be evaluated during regularization path
    fit_intercept : bool
        center data or not, forwarded to sklearn
    cv : int
        how many folds to carry out in cross-validation

    Returns
    -------
    results : dict
        dictionary containing parameters,
        alpha_path (all tested alpha values),
        mse_path (mse for validation set for each alpha),
        alpha_optimal (alpha value that yields the lowest validation rmse)

    '''

    if alphas is None:
        alphas = np.logspace(-9, 0.5, 100)

    lassoCV = LassoCV(alphas=alphas, fit_intercept=False, cv=cv,
                      max_iter=max_iter, tol=tol)
    lassoCV.fit(X, y)
    results = dict()
    results['parameters'] = lassoCV.coef_
    results['alpha_optimal'] = lassoCV.alpha_
    results['alpha_path'] = lassoCV.alphas_
    results['mse_path'] = lassoCV.mse_path_.mean(axis=1)
    return results


def fit_bayesian_ridge(X, y, fit_intercept=False, **kwargs):
    '''
    Return the solution `a` to the linear problem `Xa=y` obtained by using
    Bayesian ridge regression as implemented in scitkit-learn.

    Parameters
    -----------
    X : matrix / array
        fit matrix
    y : array
        target array
    fit_intercept : bool
        center data or not, forwarded to sklearn

    Returns
    ----------
    results : dict
        dict containing parameters, covariance matrix
    '''
    brr = BayesianRidge(fit_intercept=fit_intercept, **kwargs)
    brr.fit(X, y)
    results = dict()
    results['parameters'] = brr.coef_
    results['covariance-matrix'] = brr.sigma_
    return results


def fit_ardr(X, y, threshold_lambda=1e8, fit_intercept=False, **kwargs):
    '''
    Return the solution `a` to the linear problem `Xa=y` obtained by using
    the automatic relevance determination regression (ARDR) method as
    implemented in scitkit-learn.

    Parameters
    -----------
    X : matrix / array
        fit matrix
    y : array
        target array
    threshold_lambda : float
        threshold lambda parameter forwarded to sklearn
    fit_intercept : bool
        center data or not, forwarded to sklearn

    Returns
    ----------
    results : dict
        dictionary containing parameters, covariance matrix
    '''
    ardr = ARDRegression(threshold_lambda=threshold_lambda,
                         fit_intercept=fit_intercept, **kwargs)
    ardr.fit(X, y)
    results = dict()
    results['parameters'] = ardr.coef_
    results['covariance-matrix'] = ardr.sigma_
    return results
