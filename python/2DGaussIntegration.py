import sys
import os
import math
import numpy as np
sys.path.append('./FEM/src/')
import QuadParams
import ShapeFuncs


class function(object):
    def __init__(self, fun, exact):
        self.fun = fun
        self.exact = exact

''' below are random test functions '''
# test = function(lambda x,y: x+y, 1.0)
# test = function(lambda x,y: x**2+y**2, 2/3)
# test = function(lambda x,y: x**4+y**(1/3), 19/20)
# test = function(lambda x,y: math.sin(x**2)+math.cos(y/5), 1.30361)

''' below are the functions used in PGD work '''
# test = function(lambda x,y: y*math.sin(x), 1.0)
# xRange = np.linspace(0,math.pi)

# test = function(lambda x,y: y*math.sin(x**2), 0.215203862334)
# xRange = np.linspace(0,math.sqrt(2.)*math.sqrt(math.pi))

test = function(lambda x,y: y*(1 + math.sin( (x-math.sqrt(math.pi))**2 ) ), 2.667285320389)
xRange = np.linspace(0,2.*math.sqrt(math.pi))

yRange = np.linspace(0,1)

order = 2
nw,xw,w = QuadParams.QP(order)

''' General numerical integration with known explicit funtion below. '''
u_h = 0.0
du_h = 0.0
for idx in range(0,len(xRange)-1):
    xL = xRange[idx]
    xR = xRange[idx+1]
    dx = (xR-xL)/2.

    for l1 in range(0,nw):
        x = xL + (1 + xw[l1])*dx

        for idy in range(0,len(yRange)-1):
            yL = yRange[idy]
            yR = yRange[idy+1]
            dy = (yR-yL)/2.

            for l2 in range(0,nw):
                y = yL + (1 + xw[l2])*dy

                #complete integrations
                u_h = u_h + test.fun(x,y) *w[l1]*dx *w[l2]*dy

print('{0: .8f}, {1: .8f}'.format(u_h,test.exact))



''' Possible FEM implementation below '''
# u = 0.0
# u_h = 0.0
# du_h = 0.0
# for idx in range(0,len(xRange)-1):
#     xL = xRange[idx]
#     xR = xRange[idx+1]
#     dx = (xR-xL)/2.
#     # print('{0: .4f}, {1: .4f}, {2: .4f}'.format(xL,xR,dx))
#
#     for l1 in range(0,nw):
#         x = xL + (1 + xw[l1])*dx
#         xpsi,xdpsi = shape.shape(xw[l1],order)
#
#         for k1 in range(0,order):
#             for idy in range(0,len(yRange)-1):
#                 yL = yRange[idy]
#                 yR = yRange[idy+1]
#                 dy = (yR-yL)/2.
#
#                 for l2 in range(0,nw):
#                     y = yL + (1 + xw[l2])*dy
#                     ypsi,ydpsi = shape.shape(xw[l2],order)
#
#                     tmp_hval = 0.0
#                     for k2 in range(0,order):
#                         tmp_hval = tmp_hval + test.fun(x,y) *xpsi[k1]*ypsi[k2]
#
#                     #complete integrations
#                     u_h = u_h + test.fun(x,y)      *w[l1]*dx *w[l2]*dy
#
# print('{0: .8f}, {1: .8f}'.format(u_h,test.exact))
