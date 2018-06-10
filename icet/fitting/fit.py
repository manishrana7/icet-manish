"""
Wrapper module to different fit_functions provided in fit_methods.py
"""

from collections import OrderedDict
from sklearn.preprocessing import StandardScaler
from .fit_methods import (fit_least_squares,
                          fit_lasso,
                          fit_elasticnet,
                          fit_bayesian_ridge,
                          fit_ardr)


fit_methods = OrderedDict([
    ('least-squares', fit_least_squares),
    ('lasso', fit_lasso),
    ('elasticnet', fit_elasticnet),
    ('bayesian-ridge', fit_bayesian_ridge),
    ('ardr', fit_ardr),
    ])
available_fit_methods = list(fit_methods.keys())


def fit(X, y, fit_method, standardize=True, **kwargs):
    """ Wrapper function for all available fit methods.

    Parameters
    -----------
    X : matrix / array
        fit matrix
    y : array
        target array
    fit_method : string
        method to be used for training; possible choice are
        "least-squares", "lasso", "elasticnet", "bayesian-ridge", "ardr"
    standardize : bool
        whether or not to standardize the fit matrix before fitting

    Returns
    ----------
    results : dict
        dict containing parameters and possibly pther information obtained by
        the fit_method
    """

    if fit_method not in available_fit_methods:
        msg = ['Fit method not available']
        msg += ['Please choose one of the following:']
        for key in available_fit_methods:
            msg += [' * ' + key]
        raise ValueError('\n'.join(msg))

    if standardize:
        ss = StandardScaler(copy=True, with_mean=False, with_std=True)
        ss.fit_transform(X)  # change in place
        results = fit_methods[fit_method](X, y, **kwargs)
        ss.inverse_transform(X)  # change in place
        ss.transform(results['parameters'].reshape(1, -1)).reshape(-1,)
    else:
        results = fit_methods[fit_method](X, y, **kwargs)
    return results
