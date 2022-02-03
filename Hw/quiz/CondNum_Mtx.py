import numpy as np
B = np.linalg.inv(A)
num = np.shape(A)
Alist = []
Blist = []
for i in range(num[0]):
    Alist.append(abs(A[i][i]))
    Blist.append(abs(B[i][i]))
cond = max(Alist)*max(Blist)
