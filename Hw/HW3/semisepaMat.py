'''
The interaction matrix comprised that computes the potential at q1,…,qn given a vector of charges Q1,...,Qn has a special structure known as a semiseparable (or SS for short).
SS matrices show up in a variety of different applications from computing potentials, to PDE solvers, to multivariate statistics in the form of the covariance matrix.
In this problem, you will solve a linear system with an SS matrix of the following form using the Sherman-Morrison formula
[A11A21A12A22][x1x2]=[b1b2].
In general, semiseparable matrices have diagonal blocks which have full rank and off-diagonal blocks which have (approximately) low rank.
In this particular problem, we assume that the SS matrix has off-diagonal blocks with rank 1 (i.e. an outer product of two vectors). In this particular problem, A11 is a n×n matrix, and A22 is a m×m matrix. x1 and b1 are lenght n vectors.
x2 and b2 are length m vectors. A12=uvT where u is a length n vector and v is a length m vector. Further, A21=wpT, where w is a length m vector and p is a length n vector.
What do I need to do?
Using Gaussian-Elimination on the block-structure of the matrix and applying the Sherman-Morrison formula, solve the linear system for x_1 and x_2.
Solving a semiseparable system in this fashion has some perfomance benefits. 
In a comment block, please describe the expected cost, in terms of n and m, for solving the entire system without taking advantage of the semiseparable structure of the matrix, such as by calling la.solve.
In another comment block, specify the expected cost, in terms of n and m, for solving the system by taking advantage of the semiseparable structure of the matrix.
NOTE: You may assume that the cost of solving a generic p×p system is cost≈Cp3, where C is some constant that does not depend on the problem size. 
(You may also find a detailed operation count for LU in our textbook.) You need not specify C, you may simply use it as a known constant. 
For your cost analysis, please only consider the leading order terms, e.g., you may simplify 2nm2+3n2m+n3+m3+457nm+1589n as 2nm2+3n2m+n3+m3.
'''

import numpy as np
import scipy.linalg as sla

def not_allowed(*args, **kwargs):
    raise RuntimeError("You called an illegal function.")

n = np.random.randint(75,100)
m = np.random.randint(75,100)
#Setup A
A11 = np.random.randn(n,n)
A11=A11.T @ A11 + n * np.eye(n)
A22 = np.random.randn(m,m)
A22=A22.T @ A22 + m * np.eye(m)

u = np.random.randn(n)
v = np.random.randn(m)
w = np.random.randn(m)
p = np.random.randn(n)

#right-hand-sides
b = np.random.randn(n+m)

A22_inv = sla.inv(A22)
A11_inv = sla.inv(A11)

for attr in dir(sla):
    setattr(sla, attr, not_allowed)
import numpy.linalg as la
la.solve = not_allowed
la.inv = not_allowed

n = A11.shape[0]
m = A22.shape[0]
uu = np.dot(w,np.dot(p.T,np.dot(A11_inv,u)))
#The x2 function is :
#uu is w*pT*A11-1*u
#A22 - uu*vT) x2 = b2 - A21A11-1*b1(B)
#Solve  A22 z = uu
#Then Solve A22 y = B
#So x2 = y + ((vT y)/(1 - vT z))z

#1: z = A22-1 * uu
z = np.dot(A22_inv,uu)
b1 = b[0:n]
b2 = b[n:n+m]
#2: B = b2 - A21A11-1*b1
B =  b2 - np.dot(w,np.dot(p.T,np.dot(A11_inv,b1)))

#3: y = A22-1 * B  
y = np.dot(A22_inv,B)

#4: x2 = y + ((vT y)/(1 - vT z))z
denom = 1 - np.dot(v.T, z)
vty = np.dot(v.T,y)
x2 = y +  np.dot(vty/denom,z)

#for x1
ww = np.dot(u,np.dot(v.T,np.dot(A22_inv,w)))
#The x1 function is :
#A11 - ww*pT) x1 = b1 - A12A22-1*b2(B1)
#Solve  A11 zz = ww
#Then Solve A11 yy = B1
#So x1 = yy + ((pT yy)/(1 - pT zz))zz

#1: z = A11-1 * ww
zz = np.dot(A11_inv,ww)
#2: B1 = b1 - A12A22-1*b2
B1 =  b1 - np.dot(u,np.dot(v.T,np.dot(A22_inv,b2)))

#3: yy = A11-1 * B  
yy = np.dot(A11_inv,B1)

#4: x1 = y + ((pT y)/(1 - pT z))z
denom = 1 - np.dot(p.T, zz)
pty = np.dot(p.T,yy)
x1 = yy +  np.dot(pty/denom,zz)
