# Copyright 2023 Daniel Sharp
# Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation
# files (the “Software”), to deal in the Software without
# restriction, including without limitation the rights to use,
# copy, modify, merge, publish, distribute, sublicense, and/or
# sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following
# conditions:
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.

import numpy as np

def covfun(x1, x2, L=0.3, var_y=0.3):
    return var_y*np.exp(-np.abs(x1-x2)/L)

def form_KL_uniform(covfun, lo, hi, num_quad):
    pts, wts = np.polynomial.legendre.leggauss(num_quad)
    pts = (pts * (hi-lo)/2) + (hi + lo)/2
    wts *= (hi-lo)/2
    mesh1, mesh2 = np.meshgrid(pts, pts)
    C_mat = covfun(mesh1, mesh2)
    W_half = np.diag(np.sqrt(wts))
    A = W_half @ C_mat @ W_half
    lam, phi = np.linalg.eigh(A)
    psi = np.diag(1 / np.sqrt(wts)) @ phi
    lam = lam[::-1]
    psi = psi[:,::-1]
    psi = psi * (2*(psi[0,:] > 0) - 1)
    return lam, psi, pts, wts

class KLE:
    def __init__(self, lam, psi, pts, wts, covfun):
        self.lam = lam
        self.psi = psi
        self.pts = pts
        self.wts = wts
        self.covfun = covfun

def default_KLE(N_trunc, num_quad = 100):
    lam, psi, pts, wts = form_KL_uniform(covfun, 0., 1., num_quad)
    return KLE(lam[:N_trunc], psi[:,:N_trunc], pts, wts, covfun)

def eval_KLE(kl, xgrid, z = None):
    if z is None:
        rng = np.random.default_rng()
        z = rng.standard_normal(len(kl.lam))
    X,Y = np.meshgrid(kl.pts, xgrid)
    covmat = kl.covfun(X,Y) @ np.diag(kl.wts)
    psi_evals = covmat @ kl.psi @ np.diag(1/np.sqrt(kl.lam))
    return psi_evals @ z

def ExampleKLE():
    xgrid = np.linspace(0,1,101)
    n_trunc = 100
    kl = default_KLE(n_trunc, max(100, round(1.5*n_trunc)))
    rng = np.random.default_rng()
    z = rng.standard_normal(n_trunc)
    kl_eval = eval_KLE(kl, xgrid, z)
    return xgrid, kl_eval, z

xgrid, kl_eval, z = ExampleKLE()