import numpy as np
import scipy.sparse as ss
import matplotlib.pyplot as plt
import numpy.linalg as la

#def open(filename, mode="r"):
#    """
#    Only used in the autograder. No need to use this locally.
#    """
#    try:
#        data = data_files["data/"+filename]
#    except KeyError:
#        raise IOError("file not found")
#
#    # 'data' is a 'bytes' object at this point.
#
#    from io import StringIO
#    return StringIO(data.decode("utf-8"))

def readmesh(fname):
    """
    Read a mesh file and return vertics as a (npts, 2)
    numpy array and triangles as (ntriangles, 3) numpy
    array. `npts` is the number of vertices of the mesh
    and `ntriangles` is the number of triangles of the
    mesh.
    """
    with open(fname, "r") as f:
        npoints = int(next(f))
        points = np.zeros((npoints, 2))

        for i in range(npoints):
            points[i, :] = [float(x) for x in next(f).split()]

        ntriangles = int(next(f))
        triangles = np.zeros((ntriangles, 3), dtype=int)

        for i in range(ntriangles):
            triangles[i, :] = [int(x)-1 for x in next(f).split()]

    return points, triangles

def plotmesh(points, triangles, tricolors = None):
    """
    Given a list of points (shape: (npts, 2)) and triangles
    (shape: (ntriangles, 3)), plot the mesh.
    """
    plt.figure()
    plt.gca().set_aspect('equal')

    if (tricolors is None):
        plt.triplot(points[:, 0], points[:, 1], triangles, "bo-", lw = 1.0)
    else:
        plt.tripcolor(points[:, 0], points[:, 1], triangles, facecolors = tricolors, edgecolors = "k")

    return

def mesh2dualgraph(triangles):
    """
    Calculate the graph laplacian of the dual graph associated
    with the mesh given by numpy array traingles.
    """
    n, m = triangles.shape

    assert(m == 3), "Triangle should have exactly three points !!"

    G = np.zeros((n, n))

    for i, ti in enumerate(triangles):
        for j, tj in enumerate(triangles):
            ## If there is a common edge
            if (len( set(ti) - set(tj) ) == 1):
                G[i, j] = G[j, i] = -1

    for i in range(n):
        G[i, i] = -np.sum(G[i, :])

    return ss.csr_matrix(G)

def lanczos1(A, r0, m):
    n = len(r0)
    #Acknowledgement of the source: https://en.wikipedia.org/wiki/Lanczos_algorithm
    T = np.zeros((m,m))
    Q = np.zeros((n,m))
    q0 = np.zeros(n)
    bv = np.ones(n)
    r0 = r0 - (r0@bv)/(bv@bv)*bv
    beta0 = la.norm(r0)
    for k in range(m):
        q = r0/beta0
        u = A @ q
        r = u - beta0*q0
        alpha = q @ r
        r = r - alpha * q
        beta = la.norm(r)
        beta0 = beta.copy()
        r0 = r.copy()
        q0 = q.copy()
        T[k,k] = alpha
        if k != m-1:
            T[k+1,k] = beta
            T[k,k+1] = beta
        Q[:,k] = q

    return Q, T

def lanczos2(L, x0, niter):
    beta = la.norm(x0) 
    A = L.copy()
    n = x0.shape[0]
    T = np.zeros((niter,niter))
    Q = np.zeros((n,niter))
    ones = np.ones(n)
    #orthogolize
    #x0 = x0 - np.inner(x0,ones)*ones/np.inner(ones,ones)
    x0 = x0 - (x0@ones)/(ones@ones)*ones
    qlast = np.zeros(n)
    qnow = x0.copy()/la.norm(x0)
    for i in range(niter):
        unow = L@qnow
        alfa = qnow.T@unow
        unow = unow.copy() - beta*qlast - alfa*qnow
        beta = la.norm(unow)
        print(beta)
        Q[:,i] = qnow.T
        T[i][i] = alfa
        qlast = qnow.copy()
        qnow = unow.copy()/beta
        if i != niter - 1:
            T[i+1][i] = T[i][i+1] = beta
    print(T)
    w,v = la.eig(T)
    #Implement the body of this function
    return Q, T

def fiedler_ritz(Q,T):
    #Implement the body of this function
    w,v = la.eigh(T)
    #idx = -1
    #smallv = 1000000
    #for i in range(len(w)):
    #    if w[i]<= smallv:
    #        smallv = w[i]
    #        idx = i
    idx = np.argsort(w)
    fiedlerVec = Q@v[:,idx[0]]
    return fiedlerVec

def fiedler(G, k):
    """
    Calculate the fiedler vector of the graph Laplacian matrix
    'G' using 'k' niter of Lanczos algorithm.
    """
    n, m = G.shape

    assert (n == m), "Matrix should be square !!"

    x0 = np.linspace(1, n, num = n)

    ## You should complete this Lanczos function
    Q, T = lanczos2(G, x0, k)

    ## You should complete this Fiedler vector computation function
    fiedlerVec = fiedler_ritz(Q,T)

    partitionVec = np.zeros_like(fiedlerVec)
    mfeidler = np.ma.median(fiedlerVec)

    for i in range(n):
        if (fiedlerVec[i] >= mfeidler):
            partitionVec[i] = 1
        else:
            partitionVec[i] = -1

    return partitionVec

points, triangles = readmesh("mesh.1")
plotmesh(points, triangles)
G = mesh2dualgraph(triangles)
partitionVec = fiedler(G, 150)
plotmesh(points, triangles, partitionVec)
