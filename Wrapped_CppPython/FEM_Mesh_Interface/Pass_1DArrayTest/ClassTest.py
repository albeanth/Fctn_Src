import sys,os
import math
import numpy as np
import GaussInt

dum1 = GaussInt.GaussianIntegration()
a = GaussInt.Line(); b = GaussInt.Line()
a = np.array([1.,2.,3.,4.])
b = np.array([5.,6.,7.,8.])
c = 500
dum1.ImportMeshGridInfo(a,b,c)

# class Test():
#     def __init__(self, a,b):
#         self.vec = np.linspace(a,b)
#
# a = Test(0,1)
# print(a.vec)
