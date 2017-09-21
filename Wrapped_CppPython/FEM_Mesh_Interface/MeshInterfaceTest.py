import sys,os
import math
import numpy as np
import SetUpGrid as mesh
import GaussInt

# dum1 = GaussInt.GaussianIntegration()
# a = GaussInt.ArrayInt();
# tmp1 = np.array([[1,2,3],[4,5,6]])
# a = tmp1.tolist()
# dum1.Read2D(a)

def Cpp_GridPrep(grid):
    b = GaussInt.LineInt(); c = GaussInt.ArrayInt(); d = GaussInt.LineDouble()
    a = grid.nels
    b = grid.order.tolist()
    c = grid.nod.tolist()
    tmp2 = grid.xnod
    tmp3 = np.reshape(tmp2, len(tmp2))
    d = tmp3.tolist()
    e = grid.maxord.tolist()
    return(a,b,c,d,e)

def Cpp_PGDPrep(tmpX,tmpY,tmpT):
    X = GaussInt.ArrayInt();Y = GaussInt.ArrayInt();T = GaussInt.ArrayInt()
    X = tmpX.tolist()
    Y = tmpY.tolist()
    T = tmpT.tolist()
    return(X,Y,T)

Xbnds = [0.0, math.pi]; xnel = [100]; Xorder = 1
spaceX = mesh.CFEMGrid1D(Xbnds,xnel,Xorder)

Ybnds = [0.0, math.pi]; ynel = [100]; Yorder = 1
spaceY = mesh.CFEMGrid1D(Ybnds,ynel,Yorder)

Tbnds = [0.0, 1.0]; tnel = [100]; Torder = 1; delta = 1E-8
if Torder != 1: print(Red+'High order elements are not supported for time domain.'); sys.exit()
time = mesh.DFEMGrid1D(Tbnds,tnel,Torder,delta)

# Xa,Xb,Xc,Xd,Xe = Cpp_Prep(spaceX)
# tmp1 = GaussInt.GaussianIntegration()
# tmp1.Get2DInfo(Xa,Xb,Xc,Xd,Xe)

# SourceX = GaussInt.GaussianIntegration() # create instance of c++ class
Xa,Xb,Xc,Xd,Xe = Cpp_GridPrep(spaceX)
Ya,Yb,Yc,Yd,Ye = Cpp_GridPrep(spaceY)
Ta,Tb,Tc,Td,Te = Cpp_GridPrep(time)
# Sx = SourceX.Source_Integrate('t', Xa,Xb,Xc,Xd,Xe, Ya,Yb,Yc,Yd,Ye, Ta,Tb,Tc,Td,Te)
# print(np.array(list(Sx)))


tmpX = np.random.rand(2,spaceX.nels)
tmpY = np.random.rand(2,spaceY.nels)
tmpT = np.random.rand(2,time.nels)
X,Y,T = Cpp_PGDPrep(tmpX,tmpY,tmpT)

err = GaussInt.GaussianIntegration()
Error = err.Error_Integrate3D(X,Y,T, Ta,Tb,Tc,Td,Te, Xa,Xb,Xc,Xd,Xe, Ya,Yb,Yc,Yd,Ye)
print(np.array(list(Error)))
