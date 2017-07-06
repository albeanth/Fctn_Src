import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import QuadParams
import shape
try:
    from colorama import Fore, Style, init
    init(autoreset=True)
    Yellow = Fore.YELLOW; Red = Fore.RED; Green = Fore.GREEN; Cyan = Fore.CYAN; Magenta = Fore.MAGENTA
    StyDim = Style.DIM
except ImportError:
    print('\nYou should get colorama. It\'s pretty sweet.\n')
    Yellow = ''; Red = ''; Green = ''; Cyan = ''; Magenta = ''
    StyDim = ''

def IVP_1D(time,T,source,lmbda):
    '''
    This function solves a 1D initial value problem using discontinuous Galerkin finite elements.

    INPUT:
        nw, xw, w: number of integration points, nodes, and weights
        nels: number of elements
        ord: lost of polynomial order for corresponsing elements
        xnod, nod, nnodes: list of node coordinates, global numbering of nodes, and total number of nodes
    OUTPUT:
        sol: is the numerical solution of the corresponding FEM coefficients

    A. Alberti 06/2017
    '''

    nw,xw,w = QuadParams.QP(time.maxord)
    ## set mesh parameters to spatial grid
    nels = time.nels
    order = time.order
    nod = time.nod
    xnod = time.xnod
    edge = time.edge
    xedge = time.xedge
    nnodes = time.nnodes

    ## set up solution vector.
    sol = np.zeros(nnodes)

    for el in range(0,nels): # for each element...
        tL = xnod[nod[el,0]] # ...leftmost node
        tR = xnod[nod[el,order[el]-1]]  # ...rightmost node
        dt = (tR-tL)/2.  # Jacobian of transformation for nodes
        tEdgeL = xedge[edge[el,0]] # ...leftmost edge
        tEdgeR = xedge[edge[el,1]] # ...rightmost edge
        dtEdge = (tEdgeR-tEdgeL)/2.  # Jacobian of transformation for edges

        # compute element stiffness matrix and load vector
        m = np.zeros((order[el],order[el]))  # element mass matrix
        k = np.zeros((order[el],order[el]))  # element stiffness matrix
        f = np.zeros(order[el])       # element load vector
        for l in range(0,nw):
            t = tL + (1 + xw[l])*dt      # t runs in true element, xw runs in reference element
            tEdge = tEdgeL + (1 + xw[l])*dtEdge

            psi,dpsi = shape.shape(xw[l],order[el])      # calculations on ref.element
            fval = source(tEdge)

            f = f + fval * psi *w[l]*dtEdge
            k = k + np.outer(dpsi,psi)/dt * w[l]*dt
            m = m + np.outer(psi,psi) * w[l]*dt

        if el == 0: # first element, so
            f[0] = f[0] + T(time.bounds[0])
        else:
            f[0] = f[0] + sol[nod[el-1,order[el]-1]]

        # # uncomment to see the local (stiffness+mass) matrix and local load vector
        # print(Cyan+'mass mat');print(m)
        # print(Cyan+'stiffness mat');print(k)
        # print(Cyan+'RHS');print(f)
        # sys.exit()

        mat = -k + lmbda*m
        mat[order[el]-1,order[el]-1] = mat[order[el]-1,order[el]-1] + 1

        sol[nod[el,0]:nod[el,0]+order[el]] = np.linalg.solve(mat,f)

    return(sol)

def Error(sol,time,T,Tp):

    nw,xw,w = QuadParams.QP(time.maxord)
    ## set mesh parameters to spatial grid
    nels = time.nels
    order = time.order
    nod = time.nod
    xnod = time.xnod

    l2Err = 0.0; h1Err=0.0

    for el in range(0,nels): # for each element in all elements
        tL = xnod[nod[el,0]] # left endpoint on element
        tR = xnod[nod[el,order[el]-1]]  # right endpoint on element
        dt = (tR-tL)/2.  # Jacobian of transformation

        for l in range(0,nw):
            t = tL + (1 + xw[l])*dt                # x runs in true element, xw runs in reference element
            psi,dpsi = shape.shape(xw[l],order[el])      # calculations on ref.element

            uval = T(t)
            duval = Tp(t)

            uhval = 0.0; duhval = 0.0 #reset uh and duh evaluations over each element
            for k in range(0,order[el]):
                mynum = nod[el,k] # (nod=global numbering of nodes). mynum = node,k, of element,el
                uhval = uhval + sol[mynum]*psi[k]
                duhval = duhval + sol[mynum]*dpsi[k]/dt

            #complete integration over element
            l2Err = l2Err + (uval-uhval)**2 *w[l]*dt
            h1Err = h1Err + (duval-duhval)**2 *w[l]*dt

    l2Err = math.sqrt(l2Err)
    h1Err = math.sqrt(l2Err+h1Err)

    return(l2Err,h1Err)

def plot(time,sol,lmbda):

    T = lambda t: np.exp(-lmbda*t)

    plt.figure(1)
    xplot = np.linspace(time.bounds[0],time.bounds[1],time.nels*10+1)
    plt.plot(xplot,T(xplot), linewidth = 2, color = 'blue')

    nw,xw,w = QuadParams.QP(time.maxord)
    ## set mesh parameters to spatial grid
    nels = time.nels
    order = time.order
    nod = time.nod
    xnod = time.xnod

    for el in range(0,nels):
        xL = xnod[nod[el,0]] # left endpoint on element
        xR = xnod[nod[el,order[el]-1]]  # right endpoint on element
        dx = (xR-xL)/2.  # Jacobian of transformation

        xpp = np.linspace(xL,xR,10);                # set up 10 points in the true element for plotting
        xww =(xpp-xL)/dx-1;                         # xww runs in reference element
        ypp = 0*xpp; # initializes a zero array of len(xpp)

        for j in range(0,len(xpp)):
            psi,dpsi = shape.shape(xww[j],order[el])
            uhval = 0
            for k in range(0,order[el]):
                mynum = nod[el,k] # (nod=global numbering of nodes). mynum = node,k, of element,el
                uhval = uhval + sol[mynum]*psi[k]
            ypp[j]=uhval

        plt.plot(xpp,ypp, linewidth = 2, color = 'red')

    plt.title('Exact solution (blue) and IVP DFEM solution (red)')
    plt.xlabel('Time')
    plt.grid(True)

    plt.show()
