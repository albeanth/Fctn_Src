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
# NumOfElem = [4,8,16,32,64,128,256,512,1024]
NumOfElem = [128]
h = zeros(len(NumOfElem))
Problem = 5
Hetero = True
#====================#
#    CFEM BVP TEST   #
#=====================
L2Error = ones(len(NumOfElem))
H1Error = ones(len(NumOfElem))
InfError = ones(len(NumOfElem))
InfErr_Loc = ones(len(NumOfElem))
LocalError = ones(len(NumOfElem))
CFEM = BVP(Problem, Hetero)  # create instance of BVP Solver for CFEM work
for i,ne in enumerate(NumOfElem):
    if ((Problem == 0) or (Problem == 1)):
        CFEM.CFEMGrid1D(0.0, 1, ne, 1)
    elif ((Problem == 2) or (Problem == 5)):
        CFEM.CFEMGrid1D(0.0, 2, ne, 1)
    elif Problem == 3:
        CFEM.CFEMGrid1D(0.0, 3, ne, 1)
    elif Problem == 4:
        CFEM.CFEMGrid1D(0, pi, ne, 1) # set up CFEM grid
    elif Problem == 6:
        CFEM.CFEMGrid1D(0.0, 3.0*pi/2.0, ne, 1)  # set up CFEM grid
    elif Problem == 7:
        CFEM.CFEMGrid1D(0, 6.0, ne, 1) # set up CFEM grid

    CFEM.CFEM_1D()                  # call BVP CFEM solver
    # for val in zip(CFEM.xnod, CFEM.soln):
    #     print("{0:.12e} {1:.12e}".format(val[0], val[1]))
    if Problem != 6:
        CFEM.L2Error()                    # compute L2 and H1 error
        CFEM.LinfError()
        CFEM.LocalError(0.85)
        L2Error[i] = CFEM.l2Err
        H1Error[i] = CFEM.h1Err
        InfError[i] = CFEM.linf_ErrVal
        InfErr_Loc[i] = CFEM.linf_Loc
        LocalError[i] = CFEM.loc_err
        h[i] = CFEM.hmax

print("    h       L2Error   LocalError")
for val in zip(h, L2Error,LocalError):
    print("{0:.4e} {1:.4e} {2:.4e}".format(val[0], val[1], val[2]))

if len(NumOfElem) == 1:
    CFEM.Plot_Both(2, "CFEM")

if ((Problem != 5) or (Problem != 7)):
    CFEM.Spatial_Convergence(L2Error, H1Error, LocalError, h, False)

#====================#
#    DFEM BVP TEST   #
#====================#
DG_L2Error = ones(len(NumOfElem))
DG_H1Error = ones(len(NumOfElem))
DG_InfError = ones(len(NumOfElem))
DG_InfErr_Loc = ones(len(NumOfElem))
DG_LocalError = ones(len(NumOfElem))
DFEM = BVP(Problem, Hetero)
for i,ne in enumerate(NumOfElem):
    if ((Problem == 0) or (Problem == 1)):
        DFEM.DFEMGrid1D(0.0, 1, ne, 1)
    elif ((Problem == 2) or (Problem == 5)):
        DFEM.DFEMGrid1D(0.0, 2, ne, 1)
    elif Problem == 3:
        DFEM.DFEMGrid1D(0.0, 3, ne, 1)
    elif Problem == 4:
        DFEM.DFEMGrid1D(0.0, pi, ne, 1)
    elif Problem == 6:
        DFEM.DFEMGrid1D(0.0, 3.0*pi/2.0, ne, 1)
    elif Problem == 7:
        DFEM.DFEMGrid1D(0, 6.0, ne, 1) 
    DFEM.DFEM_1D()
    if Problem != 7:
        DFEM.L2Error()
        DFEM.LinfError()
        DFEM.LocalError(0.85)
        DG_L2Error[i] = DFEM.l2Err
        DG_H1Error[i] = DFEM.h1Err
        DG_InfError[i] = DFEM.linf_ErrVal
        DG_InfErr_Loc[i] = DFEM.linf_Loc
        DG_LocalError[i] = DFEM.loc_err
        h[i] = DFEM.hmax

print("    h       L2Error   LocalError")
for val in zip(h, DG_L2Error, DG_LocalError):
    print("{0:.4e} {1:.4e} {2:.4e}".format(val[0], val[1], val[2]))

if len(NumOfElem) == 1:
    DFEM.Plot_Both(2, "DFEM")

if ((Problem != 5) or (Problem != 7)):
    DFEM.Spatial_Convergence(DG_L2Error, DG_H1Error, DG_LocalError, h, True)
