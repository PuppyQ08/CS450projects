'''
Given a matrix A and vector x, use LU to calculate y=A−1x. You may use scipy.linalg.lu scipy.linalg.solve_triangular.

INPUT:

A: a square matrix.
x: a vector
OUTPUT:

y: output vector equivalent to A−1x.
'''
import scipy.linalg as la
import numpy as np
P,L,U = la.lu(A)
Pp = P.transpose()#remember permute matrix : PP^T = I
Xx = np.dot(Pp,x)
M = la.solve_triangular(L,Xx,lower = True)
y = la.solve_triangular(U,M)
