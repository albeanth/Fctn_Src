# import python packages
import sys
import os
from numpy import ones, zeros, linspace
from math import pi
import matplotlib.pyplot as plt
try:
    from colorama import Fore, Style, init
    init(autoreset=True)
    Yellow = Fore.YELLOW; Red = Fore.RED; Green = Fore.GREEN; Cyan = Fore.CYAN; Magenta = Fore.MAGENTA
    StyDim = Style.DIM
except ImportError:
    print('\nYou should get colorama. It\'s pretty sweet.\n')
    Yellow = ''; Red = ''; Green = ''; Cyan = ''; Magenta = ''
    StyDim = ''

# import user defined class
from src.BVP_Solvers import BVP_Solvers as BVP
NumOfElem = [4,8,16,32,64,128,256,512,1024]
# NumOfElem = [16]
h = zeros(len(NumOfElem))
Problem = 2
Hetero = 2
#====================#
#    CFEM BVP TEST   #
#=====================
L2Error = ones(len(NumOfElem))
H1Error = ones(len(NumOfElem))
InfError = ones(len(NumOfElem))
InfErr_Loc = ones(len(NumOfElem))
CFEM = BVP(Problem, Hetero)  # create instance of BVP Solver for CFEM work
for i,ne in enumerate(NumOfElem):
    if Problem == 1:
        CFEM.CFEMGrid1D(0, pi, ne, 1) # set up CFEM grid
    elif Problem == 2:
    CFEM.CFEMGrid1D(0.0, 3.0*pi/2.0, ne, 1)  # set up CFEM grid
    elif Problem == 3:
        CFEM.CFEMGrid1D(0, 6.0, ne, 1) # set up CFEM grid

    CFEM.General_1D()                  # call BVP CFEM solver
    CFEM.L2Error()                    # compute L2 and H1 error
    CFEM.LinfError()
    L2Error[i] = CFEM.l2Err
    H1Error[i] = CFEM.h1Err
    InfError[i] = CFEM.linf_ErrVal
    InfErr_Loc[i] = CFEM.linf_Loc
    h[i] = CFEM.hmax

    # CFEM.Plot(2, "CFEM")
CFEM.Spatial_Convergence(L2Error, H1Error, h, False)

#====================#
#    DFEM BVP TEST   #
#====================#
# DG_L2Error = ones(len(NumOfElem))
# DG_H1Error = ones(len(NumOfElem))
# DG_InfError = ones(len(NumOfElem))
# DG_InfErr_Loc = ones(len(NumOfElem))
# DFEM = BVP()
# for i,ne in enumerate(NumOfElem):
# #     # DFEM.DFEMGrid1D(0.0, pi, ne, 1)
#     DFEM.DFEMGrid1D(0.0, 3.0*pi/2.0, ne, 1)
#     # DFEM.DFEMGrid1D(0, 6.0, ne, 1)  # set up CFEM grid
#     DFEM.General_1D()
#     DFEM.L2Error()
#     DFEM.LinfError()
#     DG_L2Error[i] = DFEM.l2Err
#     DG_H1Error[i] = DFEM.h1Err
#     DG_InfError[i] = DFEM.linf_ErrVal
#     DG_InfErr_Loc[i] = DFEM.linf_Loc
#     h[i] = DFEM.hmax

# # print(DG_InfErr_Loc)
# # DFEM.Plot(2, "DFEM")
# DFEM.Spatial_Convergence(DG_L2Error, DG_H1Error, h, False)
# # plt.figure(2)
# # x = linspace(0, pi, 100)#3.0*pi/2.0)
# # plt.plot(x, CFEM.u(x), '-', linewidth=1, color='black')
# # plt.plot(CFEM.xnod, CFEM.soln, '-', linewidth=1, color='blue')
# # for idx in range(0, DFEM.nnodes, 2):
# #     plt.plot(DFEM.xnod[idx:idx+2], DFEM.soln[idx:idx+2], '-', linewidth=1, color='red')
# # plt.xlabel('Space')
# # plt.grid(True)
# # plt.show()

# plt.figure(3)
# plt.loglog(h, InfError, 'x-', linewidth=2, color='red',label='CFEM')
# plt.loglog(h, h**2/100, ':', linewidth=2, color='red', alpha=0.5)
# plt.loglog(h, DG_InfError, 'x-', linewidth=2, color='blue',label='DFEM')
# plt.loglog(h, h/10, ':', linewidth=2, color='blue', alpha=0.5)
# plt.xlabel('dx')
# plt.grid(True)
# plt.legend()
# plt.show()

