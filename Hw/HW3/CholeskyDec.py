import numpy as np
# Python3 program to decompose
# a matrix using Cholesky
# Decomposition
def Cholesky_Decomposition(matrix, n):
    lower = np.zeros((n,n))
    # Decomposing a matrix
    # into Lower Triangular
    for i in range(n):
        for j in range(i + 1):
            sum1 = 0
            if (j == i):
                for k in range(j):
                    sum1 += pow(lower[j][k], 2)
                lower[j][j] = math.sqrt(matrix[j][j] - sum1)
            else:
                for k in range(j):
                    sum1 += (lower[i][k] *lower[j][k])
                lower[i][j] = (matrix[i][j] - sum1) /lower[j][j]
    print(lower)

n = np.random.randint(100,150)
A = np.random.randn(n,n)
A = A @ A.T + n * np.eye(n,n)

def not_allowed(*args, **kwargs):
    raise RuntimeError("You called an illegal function.")

import scipy.linalg as sla
import numpy.linalg as la
import math
num  = A.shape[0]
L = np.zeros((num,num))
Cholesky_Decomposition(A,num)
print(A)
for i in range(num):
    if(i != 0):
        for j in range(i):
            temp = A[i][j]
            for k in range(j):
                temp -= L[i][k]*L[j][k]
            L[i][j] = temp/L[j][j]
    insqrt = A[i][i]
    for k in range(i):
        insqrt -= L[i][k]*L[i][k] 
    L[i][i] = math.sqrt(insqrt)
print(L)

