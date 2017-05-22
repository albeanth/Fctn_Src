import sys
import math
import numpy as np
try:
    from colorama import Fore, Style, init
    init(autoreset=True)
    Yellow = Fore.YELLOW; Red = Fore.RED; Green = Fore.GREEN; Cyan = Fore.CYAN; Magenta = Fore.MAGENTA
    StyDim = Style.DIM
except ImportError:
    print('\nYou should get colorama. It\'s pretty sweet.\n')
    Yellow = ''; Red = ''; Green = ''; Cyan = ''; Magenta = ''
    StyDim = ''

def FEM_1D_BVP(a, b, xnel, myorder):
    '''
    This function solves a two-point BVP using continuous Galerkin finite elements.
      - either first or second order FE
      - Dirichlet or Neumann BCs
      - Numerical integration is used in computation of stiffness matrix and rhs

    For the following problem:
        

    INPUT:
        <a,b> = endpoints of calculation domain
        <xnel> =  the number or location of elements
            -> if one number (scalar) then xnel = number of elements,
            -> if a vector, then it contains position of first node of each element
        <myorder> =  degree of polynomials per element
            -> if one number (scalar) then it is uniform over all elements,
            -> if a vector, then it must be as long as the vector of elements
    OUTPUT:
        <sol> is the computed numerical solution at nodes
        <xnod> position of all nodes

    A. Alberti 05/2017
    Referenced from M. Peszynska MTH659 FEM Course
    '''

    if len(xnel) == 1:  #uniform grid
        nels = xnel[0];
        hel = (b-a)/nels; # spacing of elements (uniform)
        xel = np.arange(a,b,hel)
    else:
        nels = len(xnel);
        xel=xnel;
        hel=np.zeros(nels-1);
        for k in np.arange(0,nels-1):
            hel[k] = xnel[k+1]-xnel[k] # spacing of elements (non-uniform)

    # set up the order of elements = 1 + degree of polynomial = number of degrees of freedom
    # e.g. for linears, use two degrees of freedom for each element
    order = np.zeros(nels) + myorder + 1
    maxord = max(order)
    # number  of nodes (including boundary values)
    nnodes = sum(order-1)+1

    # derive the global indexing of nodes:
    # nod(i,1) is the global number of j'th node in element i
    nod = np.zeros((nels,int(maxord)))-1;
    myel = np.zeros((int(nnodes),2));
    n = 0
    for k in np.arange(0,nels): # loop over number of elements
        for j in np.arange(0,int(order[k])):      # loop over order
            nod[k,j] = n
            if j == 0:
                myel[n,1] = k+1
            elif j == int(order[k]-1):
                myel[n,0] = k+1
            else:
                myel[n,0] = k+1
                myel[n,1] = k+1
            if j != int(order[k])-1:
                n+=1

    # uncomment to see the way the nodes and elements are numbered
    # print(nod)
    # print(myel)
    # sys.exit()

    # xnod (i=1..nnodes): coordinates of node i
    xnod = np.zeros((int(nnodes),1));
    for k in np.arange(0,nels-1):
        h = xel[k+1]-xel[k]   # h is the size of the element,
        hi = h/(order[k]-1)   # hi is the size of subdivision
        for j in np.arange(0,int(order[k])):
            xnod[int(nod[k,j])] = xel[k] + hi*(j)
    k = nels-1
    h = b-xel[k]
    hi=h/(order[k]-1)
    for j in np.arange(0,int(order[k])):
        xnod[int(nod[k,j])] = xel[k] + hi*(j)

    try:
        hmax = max(hel)
    except TypeError:
        hmax = hel

    # set up:
    #   - numerical integration for stiffness matrix
    #   - quadrature parameters on the reference element (-1,1)
    #   - number of integration points nw, nodes xw, and weights w
    xw = np.zeros(int(maxord))
    w = np.zeros(int(maxord))
    if maxord == 1:  # exact for linears
        nw = 1
        xw[0] = 0.
        w[0] = 2.
    elif maxord == 2:  # exact for cubics
        nw = 2
        xw[0] = -1/math.sqrt(3.)
        xw[1] = -xw[0]
        w[0] = 1.
        w[1] = 1.
    elif maxord == 3:  # exact for polynomials of degree 5
      nw = 3;
      xw[0] = -math.sqrt(3./5.)
      xw[1]=0.
      xw[2] =- xw[0];
      w[0] = 5./9.
      w[1]=8./9.
      w[2]=w[0];
    else:
        print(Red+'Error: no code implemented for the desired order')
        sys.exit()

    ## set up matrix and rhs of linear system
    mat = np.zeros((int(nnodes),int(nnodes)));  rhsf = np.zeros((int(nnodes),1));

    for el in np.arange(0,nels): # for each element in all elements
        x1 = xnod[int(nod[el,0])] # left endpoint on element
        x2 = xnod[int(nod[el,int(order[el])-1])]  # right endpoint on element
        dx = (x2-x1)/2.  # Jacobian of transformation

        # # uncomment to see the geometrical information about the nodes and elements
        # print('Element {0:g} with {1:g} nodes'.format(el,int(order[el])));
        # for k in np.arange(0,int(order[el])):
        #     print(xnod[int(nod[el,k])])

        # # compute element stiffness matrix and load vector
        m = np.zeros((int(order[el]),int(order[el])))  # element stiffness matrix
        f = np.zeros((int(order[el]),1))       # element load vector

        for l in np.arange(0,nw):
            x = x1 + (1 + xw[l])*dx                # x runs in true element, xw runs in reference element

            psi,dpsi = shape(xw[l],int(order[el]))      # calculations on ref.element
            fval = 2.0
            aval = 1.0

            f = f + fval * psi * w[l]*dx
            m = m + aval * (dpsi*np.transpose(dpsi))/dx/dx * w[l]*dx

        # # uncomment to see the local stiffness matrix and local load vector
        # print(m)
        # print(f)

        # add the computed element stiffness matrix and load vector to the global matrix and vector
        for idx in np.arange(0,int(order[el])):
            rhsf[int(nod[el,idx])] = rhsf[int(nod[el,idx])] + f[idx]
            for idy in np.arange(0,int(order[el])):
                mat[int(nod[el,idx]),int(nod[el,idy])] = mat[int(nod[el,idx]),int(nod[el,idy])] + m[idx][idy]

    # impose Dirichlet boundary conditions eliminate known values from the system
    sol = np.zeros((int(nnodes),1))
    sol[0] = 1.0        # BC
    sol[-1] = 1.0       # BC
    rhsf = np.subtract(rhsf , np.dot(mat,sol))

    # SOLVE! (for coefficients of finite elements, still not the \emph{actual} solution)
    sol[1:-1] = np.linalg.solve(mat[1:-1,1:-1],rhsf[1:-1])

    return(sol,xnod)
#
#
def shape(x,n):
    '''
    shape function on reference element (-1,1)
    n = 2: linear
    n = 3: quadratic
    '''
    try:
        y = np.zeros((n,len(x)))
        dy = np.zeros((n,len(x)))
    except TypeError:
        y = np.zeros((n,1))
        dy = np.zeros((n,1))
    if n == 2:
        y[0,:] = .5*(1-x)
        y[1,:] = .5*(1+x)
        dy[0,:] = -.5
        dy[1,:] = .5
    elif n == 3:
        y[0,:] = (np.subtract(np.square(x),x))/2
        y[1,:] = 1-np.square(x)
        y[2,:] = (np.add(np.square(x),x))/2
        dy[0,:] = x-1/2
        dy[1,:] = -2*x
        dy[2,:] = x+1/2
    else:
        print(Red+'Error: Implementation for n='+str(n)+' not provided!')
        sys.exit()
    return(y,dy)
