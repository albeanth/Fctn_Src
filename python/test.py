import Hess
import QR
import numpy as np
import scipy
from scipy import linalg

A = [[2.,-1.,0.,0.],[-1.,2.,-1.,0.],[0.,-1.,2.,-1.],[0.,0.,-1.,2.]]
b = [0.,5.,5.,5.]
# test = Hess.hess(A)
# np.set_printoptions(precision=15)
# for row in test:
#     for elem in row:
#         print('{:.12e}'.format(elem))
    # print('\n')

# Q,R,err = QR.QR(A)
# print('Error in QR factorization = '+str(err))

x = QR.QRsolve(A,b)
print(x)
