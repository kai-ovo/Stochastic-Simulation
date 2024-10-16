import numpy as np

def hermitePoly(x, degree):
    """
    psi = hermitePoly(x, degree)
    Evaluates hermite polynomials up to degree `degree` on all points x

    Arguments:
        x: points to evaluate at
        degree: highest degree to evaluate to

    Returns
        psi: (degree+1,length(x)) Evaluate degree+1 hermite polynomials on x
    """
    x = x[:]
    c = np.identity(degree+1)
    psi = np.polynomial.hermite_e.hermeval(x, c)
    return psi