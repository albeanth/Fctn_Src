import fctns
import numpy as np
import scipy
from scipy import linalg

# A = [[1.,2.,3.],[4.,5.,6.],[7.,8.,9.]]
A = [[1.,2.,3.,4.],[5.,6.,7.,8.],[9.,10.,11.,12.],[13.,14.,15.,16.]]
test = fctns.hess(A)
# test_stckex = fctns.hess_stackex(A)
np.set_printoptions(precision=15)
for row in test:
    for elem in row:
        print('{:.12e}'.format(elem))
    print('\n')
# print(1*test[0][-1])
# sc = scipy.linalg.hessenberg(A)
# print(sc)
# print(1*sc[0][-1])
# print(test-sc)
