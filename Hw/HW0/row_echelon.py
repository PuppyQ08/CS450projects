import numpy as np
def fnn(A):
    num = A.shape[0]
    B = A
    for i in range(num-1):
        for j in range(i+1, num):
            B[j,:] = (A[j,:]) -A[i,:]*A[j,i]/A[i,i]
            B[i,:] = B[i,:]/A[i,i]
        
A = np.array([[3, -9, 12],[3,-7,8],[0,3,-6]])
#A = np.array([[2.,1.,-1.],[-3.,-1.,2.],[-2.,1.,2.]])
fnn(A)
