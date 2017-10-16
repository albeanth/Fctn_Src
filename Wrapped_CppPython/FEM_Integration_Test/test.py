import FEMInt
from math import sin,cos,pi
import sys
import numpy as np
sys.path.append('/Users/tonyalberti/Documents/projects/Fctn_Src/python/FEM_test')
import SetUpGrid
import BVP_Test as BVP

def Cpp_GridPrep(grid):
    b = FEMInt.LineInt(); c = FEMInt.ArrayInt(); d = FEMInt.LineDouble()
    a = grid.nels
    b = grid.order.tolist()
    c = grid.nod.tolist()
    tmp2 = grid.xnod
    tmp3 = np.reshape(tmp2, len(tmp2))
    d = tmp3.tolist()
    e = grid.maxord.tolist()
    return(a,b,c,d,e)

#=====================#
#     CFEM BVP TEST   #
#=====================#
u = lambda x: sin(x) #cos(x**2) #x**2 #
up = lambda x: cos(x) #-sin(x**2)*2*x # 2*x #
upp = lambda x: -sin(x) #(-cos(x**2)*(2*x)**2 - sin(x**2)*2)
k = lambda x: x+1
kp = 1
source = lambda x: -( kp*up(x) + k(x)*upp(x) )
NumOfElem = [200]

testGrid = SetUpGrid.CFEMGrid1D([0.0,pi],NumOfElem,1)
FESoln = BVP.CFEM_BVP_1D(testGrid,u,source,k)

Xa,Xb,Xc,Xd,Xe = Cpp_GridPrep(testGrid)
err = FEMInt.GaussianIntegration()
Error = err.Error_Integrate1D(FESoln,Xa,Xb,Xc,Xd,Xe)
print(Error)
print()

TestIntegral = err.FEM_Func_Integrate_1D(FESoln,Xa,Xb,Xc,Xd,Xe)
print(TestIntegral)
print()

TestIntegral = err.FEM_Func_Integrate_2D(FESoln,Xa,Xb,Xc,Xd,Xe,Xa,Xb,Xc,Xd,Xe)
print(TestIntegral)
print()

# BVP.plot(testGrid,FESoln)
# sys.exit()
# l2Err, h1Err = BVP.Error(testGrid,FESoln,u,up)
# print(l2Err)
# integL2 = FunT.L2_Norm(testGrid,FESoln)
