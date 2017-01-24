import Hess
import QR
import numpy as np
import scipy
from scipy import linalg

# A = [[1.,2.,3.],[4.,5.,6.],[7.,8.,9.]]
A = [[1.,2.,3.,4.],[5.,6.,7.,8.],[9.,10.,11.,12.],[13.,14.,15.,16.]]
# test = Hess.hess(A)
# np.set_printoptions(precision=15)
# for row in test:
#     for elem in row:
#         print('{:.12e}'.format(elem))
    # print('\n')

Q,R,err = QR.QR(A)
print('Error in QR factorization = '+str(err))
