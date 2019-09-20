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
NumOfElem = [32]
h = zeros(len(NumOfElem))
Problem = 4
Hetero = 4
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
    elif Problem == 4:
        CFEM.CFEMGrid1D(0, 6.0, ne, 1) # set up CFEM grid
    CFEM.General_1D()                  # call BVP CFEM solver

#====================#
#    DFEM BVP TEST   #
#====================#
DG_L2Error = ones(len(NumOfElem))
DG_H1Error = ones(len(NumOfElem))
DG_InfError = ones(len(NumOfElem))
DG_InfErr_Loc = ones(len(NumOfElem))
DFEM = BVP(Problem, Hetero)
for i,ne in enumerate(NumOfElem):
    if Problem == 1:
        DFEM.DFEMGrid1D(0.0, pi, ne, 1)
    elif Problem == 2:
        DFEM.DFEMGrid1D(0.0, 3.0*pi/2.0, ne, 1)
    elif Problem == 4:
        DFEM.DFEMGrid1D(0, 6.0, ne, 1) 
    DFEM.General_1D()

plt.figure(2)
plt.plot(CFEM.xnod, CFEM.soln, '.-', linewidth=1, color='red')
for idx in range(0, DFEM.nnodes, 2):
    plt.plot(DFEM.xnod[idx:idx+2], DFEM.soln[idx:idx+2], '.-', linewidth=1, color='blue')
plt.xlabel('Space')
plt.grid(True)
plt.show()

