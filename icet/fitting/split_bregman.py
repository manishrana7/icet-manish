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
    sensing as a new paradigm for model building) """
    nCols = A.shape[1]
    d = np.squeeze(np.zeros((nCols, 1)))
    b = np.squeeze(np.zeros((nCols, 1)))
    x = np.squeeze(np.zeros((nCols, 1)))

    oldNorm = 0.0

    # Precompute for speed.
    AtA = np.dot(A.conj().transpose(), A)
    ftA = np.dot(y.conj().transpose(), A)
    ii = 0
    for i in range(n_iters):
        if verbose:
            print('Iteration ', i)
        args = (A, y, mu, lmbda, d, b, AtA, ftA)
        res = minimize(objectivefunction, x, args, method="BFGS", options={
                       'disp': False}, jac=objectivefunctionDer)
        x = res.x

        d = shrink(mu*x + b, 1.0/lmbda)
        b = b + mu*x - d

        newNorm = np.linalg.norm(x)
        ii = ii + 1

        if verbose:
            print('|newNorm-oldNorm| = ', abs(newNorm-oldNorm))
        if abs(newNorm-oldNorm) < tol:
            break

        oldNorm = newNorm
    else:
        print("Warning: split bregman ran for max iters")

    tmp = np.dot(A, x) - y
    res = np.dot(tmp.conj().transpose(), tmp)
    fit_results = {'parameters': x,
                   'iters_run': ii
                   }
    return fit_results

# objective function to minimize by BFGS


def objectivefunction(x, A, y, mu, lmbda, d, b, AtA, ftA):
    tmp1 = np.dot(A, x) - y
    term1 = 0.5*np.vdot(tmp1, tmp1)

    tmp2 = d - b - mu*x
    term2 = 0.5*lmbda*np.vdot(tmp2, tmp2)
    if term1.imag > 0.0:
        print(term1.imag)
        print('Objective function contains non-zero imaginary part. ')
    if term2.imag > 0.0:
        print(term2.imag)
        print('Objective function contains non-zero imaginary part. ')

    return term1 + term2

# The derivative of the analytical function above


def objectivefunctionDer(x, A, y, mu, lmbda, d, b, AtA, ftA):
    """
    More things can be precomputed in this
    function, but you wont gain that much
    since the heavy matrix operations are alread precomputed.
    """
    ret = np.squeeze(np.dot(x[np.newaxis, :], AtA) -
                     ftA - lmbda*mu*(d - mu * x - b))
    return ret


"""
Shrink operator as defined in Eq. (11) (p. 5)
in Nelson, Hart (Compressive sensing as a new
paradigm for model building)
"""


def shrink(y, alpha):
    return np.sign(y)*np.maximum(np.abs(y)-alpha, 0.0)
