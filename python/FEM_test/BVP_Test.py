import sys
import numpy as np
import math
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

def CFEM_BVP_1D(space,T,source,cond):
    '''
    This function solves a two-point BVP using continuous Galerkin finite elements.
      - either first or second order FE
      - Numerical integration is used in computation of stiffness matrix and rhs

    INPUT:
        nw, xw, w: number of integration points, nodes, and weights
        nels: number of elements
        ord: lost of polynomial order for corresponsing elements
        xnod, nod, nnodes: list of node coordinates, global numbering of nodes, and total number of nodes
    OUTPUT:
        sol: is the numerical solution of the corresponding FEM coefficients

    A. Alberti 06/2017
    Referenced from M. Peszynska MTH659 FEM Course
    '''

    nw,xw,w = QuadParams.QP(space.maxord)

    ## set mesh parameters to spatial grid
    nels = space.nels
    order = space.order
    nod = space.nod
    xnod = space.xnod
    nnodes = space.nnodes

    ## set up matrix and rhs of linear system
    stiff = np.zeros((nnodes,nnodes)); rhsf = np.zeros(nnodes)#rhsf = np.zeros((nnodes,1));
    for el in range(0,nels): # for each element in all elements
        xL = xnod[nod[el,0]] # left endpoint on element
        xR = xnod[nod[el,order[el]-1]]  # right endpoint on element
        dx = (xR-xL)/2.  # Jacobian of transformation

        # compute element stiffness matrix and load vector
        k = np.zeros((order[el],order[el]))  # element stiffness matrix
        f = np.zeros(order[el])     # element load vector

        for l in range(0,nw):
            x = xL + (1 + xw[l])*dx                # x runs in true element, xw runs in reference element
            # print('x -> '+str(x))
            psi,dpsi = shape.shape(xw[l],order[el])      # calculations on ref.element
            kval = cond(x)
            fval = source(x)
            f = f + fval * psi *w[l]*dx
            k = k + kval*np.outer(dpsi, dpsi)/dx/dx * w[l]*dx

        # # uncomment to see the local (stiffness+mass) matrix and local load vector
        # print(k)
        # print(f)
        # sys.exit()

        # add the computed element stiffness matrix and load vector to the global matrix and vector
        for idx in range(0,order[el]):
            rhsf[nod[el,idx]] = rhsf[nod[el,idx]] + f[idx]
            for idy in range(0,order[el]):
                stiff[nod[el,idx],nod[el,idy]] = stiff[nod[el,idx],nod[el,idy]] + k[idx][idy]

    # impose Dirichlet boundary conditions eliminate known values from the system
    sol = np.zeros(nnodes)
    sol[0] =  T(space.bounds[0])  #space.bounds[0]    # BC
    sol[-1] = T(space.bounds[1])  #space.bounds[1]      # BC
    mat = stiff
    # rhsf = np.subtract(rhsf , np.dot(mat,sol))

    # SOLVE! (for coefficients of finite elements, still not the \emph{actual} solution)
    sol[1:-1] = np.linalg.solve(mat[1:-1,1:-1],rhsf[1:-1])

    return(sol)

def CFEM_Func_1D(space,fun):
    '''
    This function solves a two-point BVP using continuous Galerkin finite elements.
      - either first or second order FE
      - Numerical integration is used in computation of stiffness matrix and rhs

    INPUT:
        nw, xw, w: number of integration points, nodes, and weights
        nels: number of elements
        ord: lost of polynomial order for corresponsing elements
        xnod, nod, nnodes: list of node coordinates, global numbering of nodes, and total number of nodes
    OUTPUT:
        sol: is the numerical solution of the corresponding FEM coefficients

    A. Alberti 05/2018
    Referenced from M. Peszynska MTH659 FEM Course
    '''
    nw,xw,w = QuadParams.QP(space.maxord)
    ## set mesh parameters to spatial grid
    nels = space.nels
    order = space.order
    nod = space.nod
    xnod = space.xnod
    nnodes = space.nnodes

    ## set up matrix and rhs of linear system
    mass = np.zeros((nnodes,nnodes)); rhsf = np.zeros(nnodes)#rhsf = np.zeros((nnodes,1));
    for el in range(0,nels): # for each element in all elements
        xL = xnod[nod[el,0]] # left endpoint on element
        xR = xnod[nod[el,order[el]-1]]  # right endpoint on element
        dx = (xR-xL)/2.  # Jacobian of transformation
        # compute element mass matrix and load vector
        m = np.zeros((order[el],order[el]))
        f = np.zeros(order[el])     # element load vector

        for l in range(0,nw):
            x = xL + (1 + xw[l])*dx                # x runs in true element, xw runs in reference element
            psi,dpsi = shape.shape(xw[l],order[el])      # calculations on ref.element
            f += fun(x) * psi *w[l]*dx
            m += np.outer(psi,psi) * w[l]*dx

        # add the computed element stiffness matrix and load vector to the global matrix and vector
        for idx in range(0,order[el]):
            rhsf[nod[el,idx]] += f[idx]
            for idy in range(0,order[el]):
                mass[nod[el,idx],nod[el,idy]] += m[idx][idy]

    # impose Dirichlet boundary conditions eliminate known values from the system
    sol = np.zeros(nnodes)
    sol[0] = fun(space.bounds[0])  # BC
    sol[-1] = fun(space.bounds[1])    # BC
    mat = mass
    # rhsf = np.subtract(rhsf , np.dot(mat,sol))

    # SOLVE! (for coefficients of finite elements, still not the \emph{actual} solution)
    sol[1:-1] = np.linalg.solve(mat[1:-1,1:-1],rhsf[1:-1])

    return(sol)

def DFEM_BVP_1D(space,T,source,cond):
    '''
    This function solves a two-point BVP using discontinuous Galerkin finite elements.

    INPUT:
        nw, xw, w: number of integration points, nodes, and weights
        nels: number of elements
        ord: lost of polynomial order for corresponsing elements
        xnod, nod, nnodes: list of node coordinates, global numbering of nodes, and total number of nodes
    OUTPUT:
        sol: is the numerical solution of the corresponding FEM coefficients

    A. Alberti 07/2017
    '''

    nw,xw,w = QuadParams.QP(space.maxord)
    ## set mesh parameters to spatial grid
    nels = space.nels
    order = space.order
    nod = space.nod
    xnod = space.xnod
    edge = space.edge
    xedge = space.xedge
    nnodes = space.nnodes

    ## set up matrix and rhs of linear system
    stiff = np.zeros((nnodes,nnodes)); rhsf = np.zeros(nnodes)#rhsf = np.zeros((nnodes,1));
    # impose dirichlet boudnary conditions on solution
    sol = np.zeros(nnodes)
    sol[0] =  T(space.bounds[0]) # space.bounds[0] BC
    sol[-1] = T(space.bounds[1])   #space.bounds[1]   # BC

    for el in range(0,nels): # for each element...
        xL = xnod[nod[el,0]] # ...leftmost node
        xR = xnod[nod[el,order[el]-1]]  # ...rightmost node
        dx = (xR-xL)/2.  # Jacobian of transformation for nodes
        xEdgeL = xedge[edge[el,0]] # ...leftmost edge
        xEdgeR = xedge[edge[el,1]] # ...rightmost edge
        dxEdge = (xEdgeR-xEdgeL)/2.  # Jacobian of transformation for edges

        # compute element stiffness matrix and load vector
        k = np.zeros((order[el],order[el]))  # element stiffness matrix
        f = np.zeros(order[el])       # element load vector
        lmbda = 1
        for l in range(0,nw):
            x = xL + (1 + xw[l])*dx      # t runs in true element, xw runs in reference element
            xEdge = xEdgeL + (1 + xw[l])*dxEdge

            psi,dpsi = shape.shape(xw[l],order[el])      # calculations on ref.element

            kval = cond(x)
            fval = source(x)
            f = f + fval * psi *w[l]*dxEdge
            k = k + kval*np.outer(dpsi, dpsi)/dx/dx * w[l]*dx

        if el == 0: # first element
            f[0] = f[0] - T(space.bounds[0])
            k[order[el]-1,order[el]-1] = k[order[el]-1,order[el]-1] - 1.0
        elif el==nels-1: # last element
            f[1] = f[1] + T(space.bounds[1])
        else:
            k[order[el]-1,order[el]-1] = k[order[el]-1,order[el]-1] - 1.0

        # print(Cyan+'k');print(mat)
        # print(Cyan+'RHS');print(f)

        # add the computed element stiffness matrix and load vector to the global matrix and vector
        for idx in range(0,order[el]):
            rhsf[nod[el,idx]] = rhsf[nod[el,idx]] + f[idx]
            for idy in range(0,order[el]):
                stiff[nod[el,idx],nod[el,idy]] = stiff[nod[el,idx],nod[el,idy]] + k[idx][idy]
        if el!=nels-1:
            stiff[nod[el,idx]+1,nod[el,idy]] = 1.0


    # # uncomment to see the global stiffness matrix and load vector
    # print(Cyan+'s stiffness mat');print(stiff)
    # print(Cyan+'RHS');print(rhsf)
    # sys.exit()

    rhsf = np.subtract(rhsf , np.dot(stiff,sol))

    sol[1:-1] = np.linalg.solve(stiff[1:-1,1:-1],rhsf[1:-1])

    return(sol)

def Error(space,sol,T,Tp):

    nw,xw,w = QuadParams.QP(space.maxord)
    ## set mesh parameters to spatial grid
    nels = space.nels
    order = space.order
    nod = space.nod
    xnod = space.xnod

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

def Linfty(space,sol,T):
    nw,xw,w = QuadParams.QP(space.maxord)
    ## set mesh parameters to spatial grid
    nels = space.nels
    order = space.order
    nod = space.nod
    xnod = space.xnod
    MaxErr = 0.0
    for el in range(0,nels): # for each element in all elements
        tL = xnod[nod[el,0]] # left endpoint on element
        tR = xnod[nod[el,order[el]-1]]  # right endpoint on element
        dt = (tR-tL)/2.  # Jacobian of transformation
        for l in range(0,nw):
            t = tL + (1 + xw[l])*dt                # x runs in true element, xw runs in reference element
            psi,dpsi = shape.shape(xw[l],order[el])      # calculations on ref.element
            uval = T(t)
            uhval = 0.0
            for k in range(0,order[el]):
                mynum = nod[el,k] # (nod=global numbering of nodes). mynum = node,k, of element,el
                uhval = uhval + sol[mynum]*psi[k]
            if abs(uval-uhval) > MaxErr:
                MaxErr = abs(uval-uhval)

    return(MaxErr)

def plot(space,sol):

    T = lambda x: x**2#np.sin(x)

    plt.figure(1)
    xplot = np.linspace(space.bounds[0],space.bounds[1],space.nels*10+1)
    plt.plot(xplot,T(xplot), linewidth = 2, color = 'blue')

    nw,xw,w = QuadParams.QP(space.maxord)
    ## set mesh parameters to spatial grid
    nels = space.nels
    order = space.order
    nod = space.nod
    xnod = space.xnod

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

        plt.plot(xpp,ypp,linewidth = 2, color = 'red')

    plt.title('Exact solution (blue) and BVP FEM solution (red)')
    plt.xlabel('Space')
    plt.grid(True)

    plt.show()
