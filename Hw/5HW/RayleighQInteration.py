#get the eigenvalue of A with b
'''
Attention:
    This code is not available to run since householder function is not included
'''
eig_vals = []
def eigv(A,b):
    return np.inner(b,A@b)/np.inner(b,b)

def geteigv(A):#get the first eigenvalue of A matrix
    n = A.shape[0]
    xnow = np.random.rand(n)
    vnow = eigv(A,xnow)
    vprev= 0
    I = np.eye(n)
    while la.norm(A@xnow - vnow*xnow) > 1e-8:
        xprev = xnow.copy()
        vprev = vnow
        xnow = la.solve(A - vprev*I,xprev)
        xnow = xnow/la.norm(xnow)
        vnow = eigv(A,xnow)
    return vnow,xnow

def getH(x):#get the H of corresponding x
    H = []
    for i in range(x.shape[0]):
        Hn = np.zeros(x.shape[0])
        Hn[i] = 1
        H.append(householder_vec(Hn,x))
    H = np.array(H).T
    return H

def getvecval(M,x1):#pass in a submatrix, eigenvector and eigenvalue of M
    HMH = apply_householder(M,x1)
    #get value of B, b and L1
    submatx = HMH[1:,1:].copy()
    if submatx.shape[0] == 1:
        lamd = HMH[1:,1:]
        eigvec = np.zeros(2)
        eig_vals.append(lamd)
        eigvec[1] = 1
        eigvec[0] = -M[0,1]/(M[0,0] - lamd)
        eigvec = eigvec/la.norm(eigvec)
        return np.append(np.array([x1]),np.array([eigvec]),axis = 0)
    b = HMH[0,1:].copy() 
    L1 = HMH[0,0]
    #diagonalize submatx By = L2y
    L2,y = geteigv(submatx)
    eig_vals.append(L2)
    #get alpha =  b.Ty/L2-L1
    alpha = np.inner(b,y)/(L2-L1)
    #get H of M:
    H = getH(x1) 
    #here enter the iteration, we got an array of submatrix's eigenvector
    retrn = getvecval(submatx[:,:],y)
    #merge y into return
    #retrn = np.append(retrn,np.array([y]),axis=0)
    #form right hand side from alpha and vector
    nextretrn = np.array([x1])
    for each in retrn:
        rhs = np.concatenate((np.array([alpha]),each))
        #get x2,which is another eigen vector of M
        x2 = la.solve(H,rhs) 
        x2 = x2/la.norm(x2)
        nextretrn = np.append(nextretrn,np.array([x2[:]]),axis = 0)
    return nextretrn[:,:] 
    
n = A.shape[0]
L,x = geteigv(A[:,:])
eig_vals.append(L)
eig_vecs = getvecval(A,x).T
