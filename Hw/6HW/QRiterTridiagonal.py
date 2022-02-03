'''
Performing QR at each step is computationally expensive. We observe at each iteration if the matrix Ai is a tridiagonal matrix then Ai+1 would also be a tridiagonal matrix and the QR factorization of a tridiagonal matrix can be computed using Givens Rotations.

In this problem you will compute the eigenvalues and eigenvectors of a symmetric matrix A, for which you would first use Householder transformations to reduce the matrix A to a tridiagonal matrix T. Then you will run QR iteration with a shift on T. The shift at every iteration is equal to the (n,n)-th entry of the matrix Ai. Stop your iterative process when the maximum change in the eigenvalues between iterations falls below 10−10.

For your assistance the setup code also includes a function for computing the householder vector called Householder(A, col, row_start, row_end). You may readily invoke this function in your solution.
'''
import numpy as np
import scipy.linalg as spla
import numpy.linalg as la

n = 4
A = np.random.randn(n, n)
A = A + A.T  # making the matrix symmetric

def Householder(H, col, row_start, row_end):
    """
    Returns an instance of ``numpy.ndarray`` of size equal to
    ``(row_end-row_start+1)``, corresponding to the Householder reflection vector
    for the given column vector of the matrix ``H``.

    :arg H: An instance of ``numpy.ndarray`` correspsonding to the matrix
    whose householder reflector is to be computed.

    :arg col: An instance of ``int`` corresponding to the column index of
    the vector whose householder reflector is to be found.

    :arg row_start, row_end: Denotes the start and end point indices(both
    inclusive) for the column vector.

    :returns v: the householder vector correspsonding to the given parameters.
    """

    v = np.zeros(row_end-row_start+1)
    v[0] = -np.sign(H[row_start, col])*np.linalg.norm(H[row_start:row_end+1, col])
    v += H[row_start:row_end+1, col]
    return v

import numpy as np
import math
def make_tridiagonal_through_hh(A):
    n = A.shape[0]
    Q = np.eye(n)
    M = A.copy()
    for i in range(n-2):
        v = Householder(M,i,i+1,n-1)
        H = np.eye(n-i-1) - 2*np.outer(v,v)/np.inner(v,v) 
        Hn = np.eye(n)
        Hn[i+1:,i+1:] = H
        Q = Hn@Q
        M = Hn@M@Hn.T
    Qret = Q.T
    Tret = Q@A@Q.T
    return Tret, Qret

def qr_of_tridiagonal_through_Givens(H):
    n = H.shape[0]
    M = H.copy()
    Qret = np.eye(n)
    for i in range(n-1):
        Q = np.eye(n) 
        r = math.sqrt(M[i][i]**2 + M[i+1][i]**2)
        Q[i][i] =Q[i+1][i+1]= M[i][i]/r
        Q[i][i+1] = M[i+1][i]/r
        Q[i+1][i] = - M[i+1][i]/r
        M = Q@M
        Qret = Q@Qret
    return Qret,M
# checking functions work properly
T_hh, Q_hh = make_tridiagonal_through_hh(A)
Q_Givens, R_Givens = qr_of_tridiagonal_through_Givens(T_hh)

# write code for QR iteration below for computing eigenpairs
n = A.shape[0]
def getDia(M):
    n =  M.shape[0]
    ret = np.zeros(n)
    for i in range(n):
        ret[i] = M[i][i]
    return ret

Mold = T_hh.copy()
I = np.eye(n)
diff = 1
while diff > 1e-10:
    eigold = getDia(Mold)
    shft = Mold[n-1][n-1]
    Mold = Mold - shft*I
    Q,R = qr_of_tridiagonal_through_Givens(Mold)
    Mnew = R@Q + shft*I
    eignew = getDia(Mnew)
    Mold = Mnew.copy()
    diff = abs(la.norm(eignew) - la.norm(eigold))

