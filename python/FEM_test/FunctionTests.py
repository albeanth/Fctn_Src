import sys
import math
import numpy as np
import QuadParams
import shape
import matplotlib.pyplot as plt

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
    print('hmax = '+str(h))

    plt.figure(pltCnt); pltCnt+=1
    plt.loglog(h,h**2,'-x',linewidth=2,color='blue',label='h^2')
    plt.loglog(h,L2Error,'-x',linewidth=2,color='red',label='DFEM IVP')
    plt.xlabel('h')
    plt.ylabel('L2Error')
    plt.grid(True)
    plt.legend(loc='upper left')

    plt.figure(pltCnt); pltCnt+=1
    plt.loglog(h,h,'x-',linewidth=2,color='blue',label='h')
    plt.loglog(h,H1Error,'x-',linewidth=2,color='red',label='DFEM IVP')
    plt.xlabel('h')
    plt.ylabel('H1Error')
    plt.grid(True)
    plt.legend(loc='upper left')

    # plt.show()
    return(pltCnt)


def Spatial_Convergence(L2Error, H1Error, h, flag=False):
    """
    computes the convergence of the L2/H1Error as a function of mesh refinement, h
    INPUT:
        L2Error = vector, computed as a function of h
        H1Error = vector, computed as a function of h
        h       = vector, mesh refinement
    OUTPUT:
        if flag == False
            table of convergence as a function of mesh refinement
        if flag == True
            table + plot as a function of mesh refinement
    """
    if len(h) == 1:
        print("Only one mesh discretization tested, mesh convergence plots not applicable.")
    else:
        # compute convergence and show table
        print("dx\t   L2Conv   H1Conv")
        print("{0:.3e}  --\t    --".format(h[0]))
        for i in range(1,len(L2Error)-1):
            conv_l2 = ( math.log(L2Error[i-1]) - math.log(L2Error[i]) ) / ( math.log(h[i-1]) - math.log(h[i]) )
            conv_h1 = ( math.log(H1Error[i-1]) - math.log(H1Error[i]) ) / ( math.log(h[i-1]) - math.log(h[i]) )
            print("{0:.3e}  {1:.3f}    {2:.3f}".format(h[i+1], conv_l2, conv_h1))
        if flag == True:
            plt.figure(2)
            # plt.subplot(121)
            plt.loglog(h, h**4/1000, 'x-', linewidth=2, label='O(-4)')
            plt.loglog(h, L2Error, 'x-', linewidth=2, label='L2Error')
            plt.loglog(h, h**2/10, 'x-', linewidth=2, label='O(-2)')
            plt.loglog(h, H1Error, 'x-', linewidth=2, label='H1Error')
            plt.legend(loc='upper left')
            plt.grid(True, which='both', axis='y', alpha=0.5)
            # plt.subplot(122)
            # plt.legend(loc='upper left')
            # plt.grid(True, which='both', axis='y', alpha=0.5)
            plt.show()