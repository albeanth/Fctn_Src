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

def Cpp_Prep(grid):
    b = GaussInt.LineInt(); c = GaussInt.ArrayInt(); d = GaussInt.LineDouble()
    a = grid.nels
    b = grid.order.tolist()
    c = grid.nod.tolist()
    tmp2 = grid.xnod
    tmp3 = np.reshape(tmp2, len(tmp2))
    d = tmp3.tolist()
    e = grid.maxord.tolist()
    return(a,b,c,d,e)


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

SourceX = GaussInt.GaussianIntegration() # create instance of c++ class
Xa,Xb,Xc,Xd,Xe = Cpp_Prep(spaceX)
Ya,Yb,Yc,Yd,Ye = Cpp_Prep(spaceY)
Ta,Tb,Tc,Td,Te = Cpp_Prep(time)
Sx = SourceX.Source_Integrate('t', Xa,Xb,Xc,Xd,Xe, Ya,Yb,Yc,Yd,Ye, Ta,Tb,Tc,Td,Te)
print(np.array(list(Sx)))
