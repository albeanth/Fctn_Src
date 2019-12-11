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

def CFEMGrid1D(bounds,xnel,myorder):
    a = bounds[0]
    b = bounds[1]
    '''
    grid1D sets up the discretized finite element domain over the interval [a,b]

    INPUT:
        <a,b> = endpoints of calculation domain
        <xnel> =  the number or location of elements
            -> if one number (scalar) then xnel = number of elements,
            -> if a vector, then it contains position of first node of each element
        <myorder> =  degree of polynomials per element
            -> if one number (scalar) then it is uniform over all elements,
            -> if a vector, then it must be as long as the vector of elements

    OUTPUT:
        Class object, "mesh". See class description below for attributes.

    '''

    class mesh:
        '''
        This class defines the characteristics of a 1D mesh. The following attributes are considered:
        1) bnds   = bounds of mesh
        2) nels   = number of elements
        3) maxord = maximum order of elements
        4) order  = list of polynomial order for corresponsing elements
        5) nod    = list of global numbering of nodes
        6) xnod   = list of node coordinates
        7) nnodes = notal number of nodes
        8) hmax   = maximum cell size (dx)
        '''
        def __init__(self, bounds,nels,order,maxord,nod,xnod,nnodes,hmax):
            self.bounds = bounds
            self.nels = nels
            self.maxord = maxord
            self.order = order
            self.nod = nod
            self.xnod = xnod
            self.nnodes = nnodes
            self.hmax = hmax

    if len(xnel) == 1:  #uniform grid
        nels = xnel[0]
        hel = (b-a)/nels # spacing of elements (uniform)
        xel = np.arange(a,b,hel)
    else:
        nels = len(xnel)
        xel=xnel
        hel=np.zeros(nels-1)
        for k in np.arange(0,nels-1):
            hel[k] = xnel[k+1]-xnel[k] # spacing of elements (non-uniform)

    # set up the order of elements = 1 + degree of polynomial = number of degrees of freedom
    # e.g. for linears, use two degrees of freedom for each element
    order = np.zeros((nels,), dtype=np.int) + myorder + 1
    maxord = max(order)
    # number  of nodes (including boundary values)
    nnodes = sum(order-1)+1

    # derive the global indexing of nodes:
    # nod(i,1) is the global number of j'th node in element i
    nod = np.zeros((nels,maxord), dtype=np.int)-1; # global node numbering
    n = 0
    for k in range(0,nels): # loop over number of elements
        for j in range(0,order[k]):      # loop over order
            nod[k,j] = n
            if j != order[k]-1:
                n+=1

    # # uncomment to see the way the nodes are numbered
    # print(nod)
    # sys.exit()

    # xnod , i=1..nnodes
    #  -> coordinates of node i
    xnod = np.zeros((nnodes,1));
    for k in range(0,nels-1):
        h = xel[k+1]-xel[k]   # h is the size of the element,
        hi = h/(order[k]-1)   # hi is the size of subdivision
        for j in range(0,order[k]):
            xnod[nod[k,j]] = xel[k] + hi*j
    k = nels-1
    h = b-xel[k]
    hi=h/(order[k]-1)
    for j in range(0,order[k]):
        xnod[nod[k,j]] = xel[k] + hi*j

    try:
        hmax = max(hel)
    except TypeError:
        hmax = hel

    return(mesh(bounds,nels,order,maxord,nod,xnod,nnodes,hmax))

def DFEMGrid1D(bounds,xnel,myorder,delta):
    a = bounds[0]
    b = bounds[1]
    '''
    grid1D sets up the discretized finite element domain over the interval [a,b]

    INPUT:
        <a,b> = endpoints of calculation domain
        <xnel> =  the number or location of elements
            -> if one number (scalar) then xnel = number of elements,
            -> if a vector, then it contains position of first node of each element
        <myorder> =  degree of polynomials per element
            -> if one number (scalar) then it is uniform over all elements,
            -> if a vector, then it must be as long as the vector of elements

    OUTPUT:
        Class object, "mesh". See class description below for attributes.

    '''

    class mesh:
        '''
        This class defines the characteristics of a 1D mesh. The following attributes are considered:
        1)  bounds = bounds of mesh
        2)  nels   = number of elements
        3)  order  = list of polynomial order for corresponsing elements
        4)  maxord = maximum order of elements
        5)  nod    = list of global numbering of nodes
        6)  xnod   = list of node coordinates
        7)  edge   = list of global numbering of cell edges
        8)  xedge  = list of cell edge coordinatees
        9)  nnodes = notal number of nodes
        10) hmax   = maximum cell size (dx)
        '''
        def __init__(self, bounds,nels,order,maxord,nod,xnod,edge,xedge,nnodes,hmax):
            self.bounds = bounds
            self.nels = nels
            self.order = order
            self.maxord = maxord
            self.nod = nod
            self.xnod = xnod
            self.edge = edge
            self.xedge = xedge
            self.nnodes = nnodes
            self.hmax = hmax

    if len(xnel) == 1:  #uniform grid
        nels = xnel[0]
        hel = (b-a)/nels # spacing of elements (uniform)
        xel = np.arange(a,b,hel)
    else:
        nels = len(xnel)
        xel=xnel
        hel=np.zeros(nels-1)
        for k in np.arange(0,nels-1):
            hel[k] = xnel[k+1]-xnel[k] # spacing of elements (non-uniform)

    # set up the order of elements = 1 + degree of polynomial = number of degrees of freedom
    # e.g. for linears, use two degrees of freedom for each element
    order = np.zeros((nels,), dtype=np.int) + myorder + 1
    maxord = max(order)
    nnodes = sum(order)      # number  of nodes (including boundary values)
    nedges = nels+1  # number of edges (including boundary values)

    # derive the global indexing of nodes:
    # nod(i,1) is the global number of j'th node in element i
    nod = np.zeros((nels,maxord), dtype=np.int)-1; # global node numbering
    edge = np.zeros((nels,2), dtype=np.int)-1; # global node numbering
    n = 0; e=0
    for k in range(0,nels): # for each element...
        for j in range(0,order[k]):   # ...get the order of that element
            nod[k,j] = n
            n+=1
            # print(j,k)
        edge[k,0]=e
        edge[k,1]=e+1
        e+=1

    # # uncomment to see the way the nodes and cell edges are numbered
    # print(Magenta+'nod\n'+str(nod))
    # print(Magenta+'edge\n'+str(edge))
    # sys.exit()

    # xedge, i=1..nels
    #  -> coordinates of edge i
    xedge = np.zeros(nedges);
    for k in range(0,nels-1):
        h = xel[k+1]-xel[k]   # h is the size of the element,
        for j in range(0,2): # only two edges regardless of polynomial order
            xedge[edge[k,j]] = xel[k] + h*j
    k = nels-1
    h = b-xel[k]
    for j in range(0,2):
        xedge[edge[k,j]] = xel[k] + h*j

    # xnod , i=1..nnodes
    #  -> coordinates of node i
    xnod = np.zeros(nnodes)
    for k in range(0,nels-1):
        h = (xel[k+1]-xel[k])-(2*delta)   # (size of the element) - (delta values of each side of element),
        hi = h/(order[k]-1)   # spacing between each interior node (based on polynomial order of that element)
        for j in range(0,order[k]):
            if j==0:
                xnod[nod[k,j]] = xel[k]+delta
            else:
                xnod[nod[k,j]] = xnod[nod[k,j-1]] + hi
    k = nels-1
    h = b-xel[k]-(2*delta)
    hi=h/(order[k]-1)
    for j in range(0,order[k]):
        if j==0:
            xnod[nod[k,j]] = xel[k]+delta
        else:
            xnod[nod[k,j]] = xnod[nod[k,j-1]] + hi

    try:
        hmax = max(hel)
    except TypeError:
        hmax = hel

    return(mesh(bounds,nels,order,maxord,nod,xnod,edge,xedge,nnodes,hmax))
