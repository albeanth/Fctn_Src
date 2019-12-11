import sys
import os
import math
import numpy as np
sys.path.append('./FEM/src/')
import QuadParams
import ShapeFuncs

fun = [lambda x,t: (x*t), lambda x,t: (x**2 * t/2), lambda x,t: ((x+9) * t**3)]
xRange = np.linspace(0,5)
order = 2
nw,xw,w = QuadParams.QP(order)
''' General numerical integration with known explicit funtion below. '''
u_h = 0.0
t = 0.5
for idx in range(0,len(xRange)-1):
    xL = xRange[idx]
    xR = xRange[idx+1]
    dx = (xR-xL)/2.
    for l1 in range(0,nw):
        x = xL + (1 + xw[l1])*dx
        sumfun = 0.0
        for elem in fun:
            sumfun += elem(x,t)
        #complete integrations
        u_h += math.sin(x)*(sumfun**2 - sumfun**3) *w[l1]*dx

print('{0: .8f}'.format(u_h))
