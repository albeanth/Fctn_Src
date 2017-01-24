import math
import numpy as np

def norm(Vec,size):
    tmp = 0.0
    for idx in np.arange(0,size):
        tmp += Vec[idx]**2
    tmp = math.sqrt(tmp)
    return(tmp)

def hess(A):
    '''
    This function reduces any symmetric, mxm matrix to Hessenberg form
    through orthogonality similarity transforms. For symmetric A, the
    Hessenberg form will be tridiagonal.

    It will take the original matrix, A, and return the tridiagonal form, T,
    and associated unitary matrix Q.
    '''

    dumA = np.array(A)
    n = np.shape(dumA)[1] # obtains n in (m,n) matrix
    v = []

    for k in np.arange(0,n-2):
        print(k)
        x = np.zeros((n-(k+1),1))
        for idx in np.arange(0,len(x)):
            x[idx]=dumA[k+1+idx,k]
        tmp = np.shape(x)[0]
        if np.sign(x[0][0]) == 0:
            tmpsign = 1;
        else:
            tmpsign = np.sign(x[0][0]);
        dum = tmpsign*norm(x,tmp)*np.eye(tmp,1) + x;
        v.append(dum/norm(dum,tmp));
        dumA[k+1:,k:] = dumA[k+1:,k:] - 2*np.dot(v[k],np.dot(np.transpose(v[k]),dumA[k+1:,k:]));
        dumA[:,k+1:] = dumA[:,k+1:] - 2*np.dot(np.dot(dumA[:,k+1:],v[k]),np.transpose(v[k]));

    return(dumA)
