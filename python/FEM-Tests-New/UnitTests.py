# import python packages
import sys
from numpy import ones, zeros
from math import pi
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
from BVP_Solvers import BVP_Solvers as BVP

pltCnt = 1 # counter for plots
# NumOfElem = [4,8,16,32,64,128,256]
NumOfElem = [8]
h = zeros(len(NumOfElem))
#====================#
#    CFEM BVP TEST   #
#=====================
L2Error = ones(len(NumOfElem))
H1Error = ones(len(NumOfElem))
CFEM = BVP()  # create instance of BVP Solver for CFEM work
for i,ne in enumerate(NumOfElem):
    CFEM.CFEMGrid1D(0.0, pi, ne, 1) # set up CFEM grid
    CFEM.CFEM_1D()                  # call BVP CFEM solver
    CFEM.Error()                    # compute L2 and H1 error
    L2Error[i] = CFEM.l2Err
    H1Error[i] = CFEM.h1Err
    h[i] = CFEM.hmax
    # for i in range(0, CFEM.nnodes):
    #     print("{0:.3f}\t{1:.4e}".format(CFEM.xnod[i], CFEM.soln[i]))

# CFEM.Plot()
CFEM.Spatial_Convergence(L2Error, H1Error, h, False)
sys.exit()

# # pltCnt = FunT.ConvergencePlots(h,L2Error,H1Error,pltCnt)

#====================#
#    DFEM BVP TEST   #
#====================#
DG_L2Error = ones(len(NumOfElem))
DG_H1Error = ones(len(NumOfElem))
# for i,ne in enumerate(NumOfElem):
    # grid.DFEMGrid1D(0.0, pi, ne, 1)



#     testGrid = SetUpGrid.DFEMGrid1D([0.0,pi],[ne],1,delta)
#     FESoln = BVP.DFEM_BVP_1D(testGrid,T,source,cond)
#     # BVP.plot(testGrid,FESoln)
#     # sys.exit()
#     l2Err, h1Err = BVP.Error(testGrid,FESoln,T,Tp)
#     DG_L2Error[i] = l2Err
#     DG_H1Error[i] = h1Err
#     h[i] = testGrid.hmax

# pltCnt = FunT.ConvergencePlots(h,DG_L2Error,DG_H1Error,pltCnt)
