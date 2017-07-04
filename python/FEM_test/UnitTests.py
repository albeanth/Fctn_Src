import sys
import math
import numpy as np
import matplotlib.pyplot as plt
sys.path.append('../')
import SetUpGrid
import QuadParams
import shape
import IVP_Test as IVP
import BVP_Test as BVP
try:
    from colorama import Fore, Style, init
    init(autoreset=True)
    Yellow = Fore.YELLOW; Red = Fore.RED; Green = Fore.GREEN; Cyan = Fore.CYAN; Magenta = Fore.MAGENTA
    StyDim = Style.DIM
except ImportError:
    print('\nYou should get colorama. It\'s pretty sweet.\n')
    Yellow = ''; Red = ''; Green = ''; Cyan = ''; Magenta = ''
    StyDim = ''

def ConvergencePlots(h,L2Error,H1Error,pltCnt):
    # print('hmax = '+str(h))
    # print('L2Error = '+str(L2Error))
    # print('H1Error = '+str(H1Error))

    plt.figure(pltCnt); pltCnt+=1
    plt.loglog(h,h**2,'-x',linewidth=2,color='blue',label='h^2')
    plt.loglog(h,L2Error,'-x',linewidth=2,color='red',label='DFEM')
    plt.xlabel('h')
    plt.ylabel('L2Error')
    plt.grid(True)
    plt.legend(loc='upper left')

    plt.figure(pltCnt); pltCnt+=1
    plt.loglog(h,h,'x-',linewidth=2,color='blue',label='h')
    plt.loglog(h,H1Error,'x-',linewidth=2,color='red',label='DFEM')
    plt.xlabel('h')
    plt.ylabel('H1Error')
    plt.grid(True)
    plt.legend(loc='upper left')

    # plt.show()
    return(pltCnt)


pltCnt = 1 # counter for plots
delta = 1E-6 # spacing of endpoints from cell edges within an element
NumOfElem = [4,8,16,32,64,128,256]
h = np.zeros(len(NumOfElem))
L2Error = np.ones(len(NumOfElem))
H1Error = np.ones(len(NumOfElem))

######################
##   CFEM BVP TEST  ##
#####################
# T = lambda x: x**2 #math.sin(x)
# Tp = lambda x: 2*x #math.cos(x)
# source = lambda x: -4*x # cond(x)*math.sin(x) - math.cos(x)
# cond = lambda x: x
# for i,ne in enumerate(NumOfElem):
#     testGrid = SetUpGrid.CFEMGrid1D([0.0,math.pi],[ne],1)
#     FESoln = BVP.CFEM_BVP_1D(testGrid,T,source,cond)
#     # BVP.plot(testGrid,FESoln)
#     # sys.exit()
#     l2Err, h1Err = BVP.Error(testGrid,FESoln,T,Tp)
#     L2Error[i] = l2Err
#     H1Error[i] = h1Err
#     h[i] = testGrid.hmax
#
# pltCnt = ConvergencePlots(h,L2Error,H1Error,pltCnt)

#####################
#   DFEM BVP TEST  ##
#####################
T = lambda x: x**2 #math.sin(x)
Tp = lambda x: 2*x #math.cos(x)
source = lambda x: -4*x #cond(x)*math.sin(x) - math.cos(x)
cond = lambda x: x
for i,ne in enumerate(NumOfElem):
    testGrid = SetUpGrid.DFEMGrid1D([0.0,math.pi],[ne],1,delta)
    FESoln = BVP.DFEM_BVP_1D(testGrid,T,source,cond)
    BVP.plot(testGrid,FESoln)
    sys.exit()
    l2Err, h1Err = BVP.Error(testGrid,FESoln,T,Tp)
    L2Error[i] = l2Err
    H1Error[i] = h1Err
    h[i] = testGrid.hmax

pltCnt = ConvergencePlots(h,L2Error,H1Error,pltCnt)

######################
##   DFEM IVP TEST  ##
######################
# lmbda = 1
# T = lambda t: 1/math.exp(t)
# Tp = lambda t: -1/math.exp(t)
# source = lambda t: -1/math.exp(t) * (1 + lmbda)
# for i,ne in enumerate(NumOfElem):
#     testGrid = SetUpGrid.DFEMGrid1D([0,3],[ne],1,delta)
#     FESoln = IVP.IVP_1D(testGrid,T,source,lmbda)
#     # IVP.plot(testGrid,FESoln)
#     # sys.exit()
#     l2Err, h1Err = IVP.Error(FESoln,testGrid,T,Tp)
#     L2Error[i] = l2Err
#     H1Error[i] = h1Err
#     h[i] = testGrid.hmax
#
# pltCnt = ConvergencePlots(h,L2Error,H1Error,pltCnt)

plt.show()
