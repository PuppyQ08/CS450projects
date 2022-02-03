'''
In this problem, you will implement QR factorization with column pivoting
AP=QR,
where P is a permutation matrix, Q is orthogonal, and R is upper triangular.
If A is rank-deficient (i.e. does not have full rank) then instead of computing the entire QR factorization in the usual manner, we want to compute a QR factorization approximately spans the same column space as A while using fewer columns in Q than A.
This is accomplished by column pivoting. We will want to pivot columns so that the first column of the submatrix that we are reducing has the largest Euclidean norm.
At the i-th step of QR with column pivoting:
Compute the 2-norms of each column of the i-th submatrix of A.
Swap the entire column with the largest norm into the first column of the submatrix.
Apply a Householder reflector to reduce the i-th column into upper triangular form.
'''
import numpy as np
#n = np.random.randint(100,150)
#r = np.random.randint(5,15)
n = 3
r = 2
#construct rank-deficient A
A = np.random.rand(n,n)
u,s,vt = np.linalg.svd(A)
s1 = np.array([10**(-i) for i in range(r)])
s[-r:] = s1
A = u @ np.diag(s) @ vt

print(A)
def apply_householder(Q, R, col):
    #compute Householder vector
    v = np.zeros(R[col:, col].shape)
    v[0] += np.sign(R[col, col])*np.linalg.norm(R[col:, col])#that np.sing is to make sure no cancelation happened
    v += R[col:,col]

    #apply Householder vector
    k = 2./np.inner(v,v)
    R[col:,col:] -= k * np.outer(v, v @ R[col:,col:])
    Q[:,col:] -= k * np.outer(Q[:,col:] @ v, v)


import scipy.linalg as sla
import math
QQ,RR,PP = sla.qr(A,pivoting=True)
print(PP)
print(QQ)
print(RR)

def swapM(A,col1,col2):
    if col1 != col2:
        temp = A[:,col1].copy()
        A[:,col1] = A[:,col2]
        A[:,col2] =temp

def swapV(A,col1,col2):
    if col1 != col2:
        temp = A[col1].copy()
        A[col1] = A[col2]
        A[col2] =temp

def norm2(A,col,start,rown):
    rt = 0
    for i in range(start,rown):
        rt += pow(A[i,col],2)
    return rt

numr = A.shape[0]
numc = A.shape[1]
P = np.arange(numr)
Q = np.eye(numr) 
R = A.copy() 
for i in range(numc-1):
    comp = norm2(R,i,i,numr)
    tgtcol = i
    for j in range(i+1,numc):
        comp2 = norm2(R,j,i,numr)
        if(i==0):
            print(comp2," ",j)
        if comp2 > comp:
            comp = comp2
            tgtcol = j
    swapV(P,i,tgtcol)
    swapM(R,i,tgtcol)
    apply_householder(Q,R,i)

print(P)            
print(Q)
print(R)
print(A[:,P])
print(Q@R)
