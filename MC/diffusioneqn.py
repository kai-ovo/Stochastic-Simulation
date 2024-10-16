import numpy as np
"""
Solve 1-D diffusion equation with given diffusivity field k
and left-hand flux F.

ARGUMENTS: 
    xgrid = vector with equidistant grid points
        F = flux at left-hand boundary, k*du/dx = -F 
   source = source term, either a vector of values at points in xgrid
            or a constant
  rightbc = Dirichlet BC on right-hand boundary
Domain is given by xgrid (should be [0,1])
"""
def diffusioneqn(xgrid, F, k, source, rightbc):
    N = len(xgrid)
    h = xgrid[N-1] - xgrid[N-2]
    
    A = np.zeros((N-1, N-1))
    b = np.zeros(N-1)
    
    if isinstance(source, (int, float)):
        f = -source * np.ones(N-1)
    else:
        f = -source[:N-1]
    
    A -= np.diag(2*k[:-1] + k[1:] + np.insert(k[:-2],0,k[0]))
    A += np.diag(  k[:-2] + k[1:-1], 1)
    A += np.diag(  k[:-2] + k[1:-1],-1)
    A /= 2 * h**2
    
    A[0, 1] += k[0] / h**2
    b[0] = 2 * F / h
    
    b[N-2] = rightbc * (k[N-1] + k[N-2]) / (2 * h**2)
    
    uinternal = np.linalg.solve(A, f - b)
    usolution = np.append(uinternal, rightbc)
    
    return usolution