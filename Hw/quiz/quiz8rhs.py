import numpy as np
import numpy.linalg as la
a = np.array([[10,-9],[3,0],[8,6]])
lst = []
lst.append(np.array([[-16],[-3],[-1]]))
lst.append(np.array([[2],[-15],[3]]))
lst.append(np.array([[1],[3],[16]]))
lst.append(np.array([[3],[-11],[19]]))
lst.append(np.array([[-1.9],[0],[1]]))
inve = la.inv(np.dot(a.T,a))
P = np.dot(a,np.dot(inve,a.T))
for b in lst:
    up = la.norm(np.dot(P,b),2)
    down = la.norm(b,2)
    print(up/down)

