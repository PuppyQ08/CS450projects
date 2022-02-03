'''
In general, the LU factors for a matrix A of bandwidth k will also have bandwidth k (assuming no pivoting). While you could solve Ax=b with Gaussian elimination, it would be inefficient and require O(n3) operations, where n is the number of rows or columns.



In this problem, you will write a banded solver that takes advantage of the structure of A. You will be given a list of matrices, along with the associated right-hand side vectors b, as well as the bandwith, k. These are stored in A_list, b_list, and k, respectively.



In taking advantage of the banded structure, make sure that your code doesn't perform any operations on the parts of the matrix that are zero (if you are stuck, see section 2.5.3 of Heath). No pivoting is required. You should not use scipy.linalg.lu or the like.



Note, if you want to re-use the storage space of a particular A matrix to hold L in the lower triangular portion (without the diagonal) and U in the upper triangular portion, that is fine. The autograder will only check the important diagonals of each L and U inside L_list and U_list.

To demonstrate the performance benefits from using the banded structure of a given matrix, you will time the execution of your banded solver for each matrix in A_list. You will also compute the time required to solve Ax=b for each matrix in A_list using la.solve. Make a plot comparing both sets of your timing data versus the number of rows in each matrix from A_list.


1.For each of the matrices in A_list and right-hand sides in b_list:

Compute a banded L and U and store it in the lists L_list and U_list.
Solve LUx=b for x and store the solution in x_lists, taking advantage of the banded structure of L and U.
Compute the time it takes to both compute LU and solve for x.
Compute the time it takes to execute la.solve(A,b)
2.Generate a plot containing the timing data for your banded solver and la.solve.



3.In a print statement, please answer the following questions.

Can you explain the performance trends that you observe?
How does the time complexity of your banded solver scale as a function of the number of rows, N?
INPUT:

A_list: A list of numpy arrays
b_list: List right hand sides
k: The bandwidth, an integer constant
OUTPUT:

L_list, U_list: List of LU factors from each of the matrices in A_list.
x_list: List of solutions to each system from A_list and b_list.
'''
import numpy as np
import scipy.linalg
import scipy.sparse as sp
import numpy.linalg as la
def rowrange(row, band, size):
    dowlmt = 0 if row - band < 0 else row - band
    uprlmt = size - 1 if row + band > size - 1 else row + band
    return dowlmt, uprlmt+1 
def bandLU(A, b, k):
    print(A)
    m =  A.shape[0]
    Ux = A.copy() 
    Lx = A.copy()
    #upper triangle M
    '''
    for i in range(m-1):#each leading row
        for l in range(i+1,m): #each row 
            rang = rowrange(i,k,m)
            lead = Ux[l][i]
            for j in range(rang[0],rang[1]):
                Ux[l][j] = Ux[l][j]- lead * Ux[i][j] / float(Ux[i][i])
    #lower triangle
    for i in range(m-1):#each leading col 
        leadr = rowrange(i,k,m)
        diagonl = Lx[i][i]
        for idx in range(leadr[0],leadr[1]):
            Lx[idx][i] *= 1/ diagonl
        for l in range(i+1,m): #each col 
            rang = rowrange(i,k,m)
            lead = Lx[i][l]
            for j in range(rang[0],rang[1]):
                Lx[j][l] = Lx[j][l]- lead * Lx[j][i] / float(Lx[i][i])
    '''
    #lower triangle
    for i in range(m-1):#each leading col 
        rang = rowrange(i,k,m)
        diagonl = Lx[i][i]
        for idx in range(rang[0],rang[1]):
            Lx[idx][i] *= 1/ diagonl
        for l in range(i+1,rang[1]): #each col 
            leadL = Lx[i][l]
            leadU = Ux[l][i]
            for j in range(rang[0],rang[1]):
                Lx[j][l] = Lx[j][l]- leadL * Lx[j][i] / float(Lx[i][i])
                Ux[l][j] = Ux[l][j]- leadU * Ux[i][j] / float(Ux[i][i])
    Lx[m-1][m-1] = 1

    #solve LUx = b:
    #solve Ly = b
    Yx = b.copy()
    Xx = b.copy()
    for i in range(m):
        Yx[i] = b[i]/Lx[i][i]
        size = rowrange(i,k,m)
        for j in range(i+1,size[1]):
            b[j] = b[j] - Lx[j][i]*Yx[i]
    #solve Ux = y
    for i in reversed(range(m)):
        Xx[i] = Yx[i]/Ux[i][i]
        size = rowrange(i,k,m)
        for j in range(i+1,size[1]):
            Yx[j] = Yx[j] - Ux[j][i]*Xx[i]
        
        
    return Ux,Lx,Xx

el = 11
n = 2**el

A = sp.diags([-1,2,-1],[-1,0,1],shape = (n,n)).toarray()

A_list = []
b_list= []
k = 1

for ells in range(6,el+1):
    A1 = A[:2**ells,:2**ells]
    A_list.append(A1)

    b1 = np.ones(2**ells)
    b_list.append(b1)
      
L_list = []
U_list = []
x_list = []
for A,b in zip(A_list,b_list):
    result = bandLU(A,b,k)
    L_list.append(result[0])
    U_list.append(result[1])
    x_list.append(result[2])
L_list = np.array(L_list)
U_list = np.array(U_list)
x_list = np.array(x_list)
