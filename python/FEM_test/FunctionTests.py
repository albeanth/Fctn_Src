import sys
import math
import numpy as np

def FuncEval(order,el,nod,sol,psi,dpsi,dx):
    # evaluate functions
    uhval = 0.0; duhval = 0.0
    for k in range(0,order[el]):
        mynum = nod[el,k] # (nod=global numbering of nodes). mynum = node,k, of element,el
        uhval = uhval + sol[mynum]*psi[k]
        duhval = duhval + sol[mynum]*dpsi[k]/dx

    return(uhval,duhval)

def L2_Norm(grid,sol):
    nw,xw,w = QuadParams.QP(grid.maxord)
    ## set mesh parameters to spatial grid
    nels = grid.nels
    order = grid.order
    nod = grid.nod
    xnod = grid.xnod

    tmp = 0.0
    for el in range(0,nels): # for each element in all elements
        xL = xnod[nod[el,0]] # left endpoint on element
        xR = xnod[nod[el,order[el]-1]]  # right endpoint on element
        dx = (xR-xL)/2.  # Jacobian of transformation

        for l in range(0,nw):
            psi,dpsi = shape.shape(xw[l],order[el])      # calculations on ref.element

            uhval,duhval = FuncEval(order,el,nod,sol,psi,dpsi,dx)

            #complete integration over element
            tmp = tmp + uhval *w[l]*dx

    return(tmp)

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
