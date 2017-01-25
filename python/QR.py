import sys
import math
import numpy as np

def norm(X,p):
    # calculates the p-norm of a vector
    tmp = 0.0
    for elem in X:
        tmp += (abs(elem))**p
    tmp = tmp**(1/p)
    return(tmp)

def MatrixInfNorm(A):
    # Calculates the maximum value in a matrix.
    val = 0.0
    for row in A:
        for elem in row:
            if elem > val:
                val = elem
    return(val)

def QR(A):
    # This function completes the QR factorization for any mxn matrix.
    # Q, R, and the inf-norm of Q*R-A are returned.
    tmp = np.shape(A)
    m = tmp[0]; n = tmp[1]
    R = np.array(A); Q = np.eye(m,m)

    # factorization for R uses Alg 10.1 in Trefethen & Bau
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

def QRsolve(A,b):
    # This function uses QR factoriation of A to solve a system of equations Ax=b.
    # the solution to Ax=b, x, is returned.
    tmp = np.shape(A)
    m = tmp[0]; n = tmp[1]
    R = np.array(A); Q = np.eye(m,m)
    y = np.array(b); soln = np.ones(len(b))

    # factorization for R. Alg 10.1 in Trefethen & Bau
    for k in np.arange(0,m-1):
        x = R[k:m,k]
        if np.sign(x[0]) == 0:
            tmpsign = 1
        else:
            tmpsign = np.sign(x[0])
        v = tmpsign*norm(x,2)*np.eye(1,len(x))[0] + x
        v = v/norm(v,2)
        R[k:m,k:n] = R[k:m,k:n] - 2*np.outer(v,np.dot(v,R[k:m,k:n]))
        # Implicit calculation for Q^* b . Alg 10.2 in Trefethen & Bau
        y[k:m] = y[k:m] - 2*np.dot(v,np.dot(v,y[k:m]))

    # Solve Rx=y for solution via back substitution
    for j in np.arange(m-1,0-1,-1):
        tmp = 0.0
        for k in np.arange(j,m-1):
            tmp += soln[k+1]*R[j,k+1]
        soln[j] = (y[j]-tmp)/R[j,j]

    return(soln)
