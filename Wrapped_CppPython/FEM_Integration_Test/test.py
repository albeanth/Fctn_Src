import FEMInt
from math import sin,cos,pi
import sys
import time
import numpy as np
sys.path.append('/nfs/depot/nerc_u1/albertia/projects/PhD/Fctn_Src/python/FEM_test')
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
u = lambda x: cos(x**2) #sin(x) #x**2 #
up = lambda x: -sin(x**2)*2*x #cos(x) # 2*x #
upp = lambda x: (-cos(x**2)*(2*x)**2 - sin(x**2)*2) #-sin(x) 
k = lambda x: x+1
kp = 1
source = lambda x: -( kp*up(x) + k(x)*upp(x) )
NumOfElem = [500]

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

start = time.perf_counter()
TestIntegral_Serial = err.FEM_Func_Integrate_2D_Serial(FESoln,Xa,Xb,Xc,Xd,Xe,Xa,Xb,Xc,Xd,Xe)
end = time.perf_counter()
SerialTime = end-start
print('Serial Calulation took {0:.4f}'.format(SerialTime))

Threads = [2,4,6,8,12,16,20,24,32]
ParallelTime = np.zeros(len(Threads))
TestIntegral_Parallel = np.zeros(len(Threads))

for idx,NUMTHREADS in enumerate(Threads):
  start = time.perf_counter()
  TestIntegral_Parallel[idx] = err.FEM_Func_Integrate_2D_Parallel(NUMTHREADS,FESoln,Xa,Xb,Xc,Xd,Xe,Xa,Xb,Xc,Xd,Xe)
  end = time.perf_counter()
  ParallelTime[idx] = end-start

print('\n\n----- Summary Table -----')
print('# Threads    Time       Speedup   F_parallel   Solution ')
print('     {0: d}      {1:.4f}       --         --      {2:.12f}'.format(1, SerialTime, TestIntegral_Serial ))
for idx in range(0,len(Threads)):
  print('     {0: d}      {1:.4f}     {2:.4f}     {3:.4f}    {4:.12f}'.format(Threads[idx], ParallelTime[idx], SerialTime/ParallelTime[idx], ((ParallelTime[idx]/SerialTime) - 1)*(Threads[idx]/(1-Threads[idx])),TestIntegral_Parallel [idx]))



# BVP.plot(testGrid,FESoln)
# sys.exit()
# l2Err, h1Err = BVP.Error(testGrid,FESoln,u,up)
# print(l2Err)
# integL2 = FunT.L2_Norm(testGrid,FESoln)
