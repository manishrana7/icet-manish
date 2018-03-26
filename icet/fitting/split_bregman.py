"""
Split Bregman Algorithm as defined on p. 5 in
Nelson, Hart (Compressive sensing as a new paradigm for model building)
"""
import numpy as np
from scipy.optimize import minimize


def split_bregman(A, y, mu=1e-3, lmbda=100, n_iters=1000, tol=1e-6, verbose=0):
    """
    Split Bregman Algorithm as defined
    on p. 5 in Nelson, Hart (Compressive
    sensing as a new paradigm for model building)

    Parameters:
    A : matrix
        Sensing matrix.
    y : numpy array
        solution vector
    mu : float
        Sparseness parameter
    lmbda : float
        Split bregmar parameter
    n_iters : int
        maximal number of split bregman iterations.
    tol : float
        tolerance for when stopping split bregman
        iterations.
    """
    n_cols = A.shape[1]
    d = np.squeeze(np.zeros((n_cols, 1)))
    b = np.squeeze(np.zeros((n_cols, 1)))
    x = np.squeeze(np.zeros((n_cols, 1)))

    old_norm = 0.0

    # Precompute for speed.
    AtA = np.dot(A.conj().transpose(), A)
    ftA = np.dot(y.conj().transpose(), A)
    ii = 0
    for i in range(n_iters):
        if verbose:
            print('Iteration ', i)
        args = (A, y, mu, lmbda, d, b, AtA, ftA)
        res = minimize(objective_function, x, args, method="BFGS", options={
                       'disp': False}, jac=objective_function_derivative)
        x = res.x

        d = shrink(mu*x + b, 1.0/lmbda)
        b = b + mu*x - d

        new_norm = np.linalg.norm(x)
        ii = ii + 1

        if verbose:
            print('|new_norm-old_norm| = ', abs(new_norm-old_norm))
        if abs(new_norm-old_norm) < tol:
            break

        old_norm = new_norm
    else:
        print("Warning: split bregman ran for max iters")

    tmp = np.dot(A, x) - y
    res = np.dot(tmp.conj().transpose(), tmp)
    fit_results = {'parameters': x}
    return fit_results


def objective_function(x, A, y, mu, lmbda, d, b):
    """
    Objective function to minimize by BFGS.

    Parameters:
    x : numpy array
        solution vector
    y : numpy array
        solution vector
    mu : float
        the parameter that adjusts sparseness.
    lmbda : float
        Split Bregman parameter
    d : numpy array
        same notation as Nelson, Hart paper
    b : numpy array
        same notation as Nelson, Hart paper
    """

    error_vector = np.dot(A, x) - y

    obj_function = 0.5*np.vdot(error_vector, error_vector)

    if obj_function.imag > 0.0:
        raise RuntimeError(
            "Objective function contains non-zero imaginary part.)")

    sparseness_correction = d - b - mu*x
    obj_function += 0.5*lmbda * \
        np.vdot(sparseness_correction, sparseness_correction)

    if obj_function.imag > 0.0:
        raise RuntimeError(
            "Objective function contains non-zero imaginary part.)")

    return obj_function


def objective_function_derivative(x, A, y, mu, lmbda, d, b, AtA, ftA):
    """
    The derivative of the objective function.

    More things can be precomputed in this
    function, but you wont gain that much
    since the heavy matrix operations are alread precomputed.

    Parameters:
    x : numpy array
        solution vector
    y : numpy array
        solution vector
    mu : float
        the parameter that adjusts sparseness.
    lmbda : float
        Split Bregman parameter
    d : numpy array
        same notation as Nelson, Hart paper
    b : numpy array
        same notation as Nelson, Hart paper
    AtA : matrix
        sensing matrix transpose times sensing matrix.
    ftA : matrix
        np.dot(y.conj().transpose(), A)

    """
    ret = np.squeeze(np.dot(x[np.newaxis, :], AtA) -
                     ftA - lmbda*mu*(d - mu * x - b))
    return ret


def shrink(y, alpha):
    """
    Shrink operator as defined in Eq. (11) (p. 5)
    in Nelson, Hart (Compressive sensing as a new
    paradigm for model building).

    y : numpy array
    alpha : float
    """
    return np.sign(y)*np.maximum(np.abs(y)-alpha, 0.0)
