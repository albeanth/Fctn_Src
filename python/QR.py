import sys
import math
import numpy as np

def norm(X,p):
    tmp = 0.0
    for elem in X:
        tmp += (abs(elem))**p
    tmp = tmp**(1/p)
    return(tmp)

def MatrixInfNorm(A):
    val = 0.0
    for row in A:
        for elem in row:
            if elem > val:
                val = elem
    return(val)

def QR(A):
    tmp = np.shape(A)
    m = tmp[0]; n = tmp[1]
    R = np.array(A); Q = np.eye(m,m)

    for k in np.arange(0,m-1):
        x = R[k:m,k]
        if np.sign(x[0]) == 0:
            tmpsign = 1
        else:
            tmpsign = np.sign(x[0])
        v = tmpsign*norm(x,2)*np.eye(1,len(x))[0] + x
        v = v/norm(v,2)
        R[k:m,k:n] = R[k:m,k:n] - 2*np.outer(v,np.dot(v,R[k:m,k:n]))
        tmpI = np.eye(m,m)
        tmpI[k:m,k:m] = tmpI[k:m,k:m] - 2*np.outer(v,v) # parenthesis around v*v' ensures it remains Hermitian
        Q = np.dot(Q,tmpI)

    err = MatrixInfNorm(np.dot(Q,R)-A)
    return(Q,R,err)
