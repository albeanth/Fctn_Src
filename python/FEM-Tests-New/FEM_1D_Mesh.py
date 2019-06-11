import sys
from numpy import zeros, arange, int, finfo

class mesh:
    """
    - This class defines the characteristics of a 1D mesh. Two options:
      1)  a continuous finite element (CFEM) mesh
      2)  a discontinuous finite element (DFEM) mesh
    """
    def __init__(self):
        """Default constructor"""

    def CFEMGrid1D(self, a, b, xnel, myorder):
        """
        sets up a discretized continuous finite element domain over the interval [a,b]
        INPUT:
            [a,b]   = bounds of calculation domain
            xnel    =  the number or location of elements
                    -> if one number (scalar) then xnel = number of elements,
                    -> if a vector, then it contains position of first node of each element
            myorder =  degree of polynomials per element
                    -> if one number (scalar) then it is uniform over all elements,
                    -> if a vector, then it must be as long as the vector of elements
        OUTPUT:
            mesh object with necessary attributes
        """
        self.bounds = [a,b]
        if isinstance(xnel,int) == True:  # uniform grid
            self.nels = xnel
            hel = (b-a)/self.nels  # spacing of elements (uniform)
            xel = arange(a, b, hel)
        else:
            self.nels = len(xnel)
            xel = xnel
            hel = zeros(self.nels-1)
            for k in arange(0, self.nels-1):
                hel[k] = xnel[k+1]-xnel[k]  # spacing of elements (non-uniform)

        # set up the order of elements = 1 + degree of polynomial = number of degrees of freedom
        # e.g. for linears, use two degrees of freedom for each element
        self.order = zeros((self.nels,), dtype=int) + myorder + 1
        self.maxord = max(self.order)
        # number  of nodes (including boundary values)
        self.nnodes = sum(self.order-1)+1

        # derive the global indexing of nodes:
        # nod(i,1) is the global number of j'th node in element i
        self.nod = zeros((self.nels, self.maxord), dtype=int)-1  # global node numbering
        n = 0
        for k in range(0, self.nels):  # loop over number of elements
            for j in range(0, self.order[k]):      # loop over order
                self.nod[k, j] = n
                if j != self.order[k]-1:
                    n += 1

        # # uncomment to see the way the nodes are numbered
        # print(self.nod)
        # sys.exit()

        # xnod , i=1..nnodes
        #  -> coordinates of node i
        self.xnod = zeros(self.nnodes)
        for k in range(0, self.nels-1):
            h = xel[k+1]-xel[k]   # h is the size of the element,
            hi = h/(self.order[k]-1)   # hi is the size of subdivision
            for j in range(0, self.order[k]):
                self.xnod[self.nod[k, j]] = xel[k] + hi*j
        k = self.nels-1
        h = b-xel[k]
        hi = h/(self.order[k]-1)
        for j in range(0, self.order[k]):
            self.xnod[self.nod[k, j]] = xel[k] + hi*j

        try:
            self.hmax = max(hel)
        except TypeError:
            self.hmax = hel

    def DFEMGrid1D(self, a, b, xnel, myorder, delta=finfo(float).resolution):
        """
        sets up a discretized discontinuous finite element domain over the interval [a,b]
        using a default discontinuity spacing of the order of floating point precision
        for the machine calling this function
        INPUT:
            [a,b]   = bounds of calculation domain
            xnel    =  the number or location of elements
                    -> if one number (scalar) then xnel = number of elements,
                    -> if a vector, then it contains position of first node of each element
            myorder =  degree of polynomials per element
                    -> if one number (scalar) then it is uniform over all elements,
                    -> if a vector, then it must be as long as the vector of elements
            delta   = 
        OUTPUT:
            mesh object with necessary attributes
        """
        self.bounds = [a, b]
        if isinstance(xnel, int) == True:  # uniform grid
            self.nels = xnel
            hel = (b-a)/self.nels  # spacing of elements (uniform)
            xel = arange(a, b, hel)
        else:
            self.nels = len(xnel)
            xel = xnel
            hel = zeros(self.nels-1)
            for k in arange(0, self.nels-1):
                hel[k] = xnel[k+1]-xnel[k]  # spacing of elements (non-uniform)

        # set up the order of elements = 1 + degree of polynomial = number of degrees of freedom
        # e.g. for linears, use two degrees of freedom for each element
        self.order = zeros((self.nels,), dtype=int) + myorder + 1
        self.maxord = max(self.order)
        self.nnodes = sum(self.order)      # number  of nodes (including boundary values)
        nedges = self.nels+1  # number of edges (including boundary values)

        # derive the global indexing of nodes:
        # nod(i,1) is the global number of j'th node in element i
        self.nod = zeros((self.nels, self.maxord), dtype=int)-1  # global node numbering
        self.edge = zeros((self.nels, 2), dtype=int)-1  # global node numbering
        n = 0
        e = 0
        for k in range(0, self.nels):  # for each element...
            for j in range(0, self.order[k]):   # ...get the order of that element
                self.nod[k, j] = n
                n += 1
                # print(j,k)
            self.edge[k, 0] = e
            self.edge[k, 1] = e+1
            e += 1

        # # uncomment to see the way the nodes and cell edges are numbered
        # print(Magenta+'nod\n'+str(self.nod))
        # print(Magenta+'edge\n'+str(self.edge))
        # sys.exit()

        # xedge, i=1..nels
        #  -> coordinates of edge i
        self.xedge = zeros(nedges)
        for k in range(0, self.nels-1):
            h = xel[k+1]-xel[k]   # h is the size of the element,
            for j in range(0, 2):  # only two edges regardless of polynomial order
                self.xedge[self.edge[k, j]] = xel[k] + h*j
        k = self.nels-1
        h = b-xel[k]
        for j in range(0, 2):
            self.xedge[self.edge[k, j]] = xel[k] + h*j

        # xnod , i=1..nnodes
        #  -> coordinates of node i
        self.xnod = zeros(self.nnodes)
        for k in range(0, self.nels-1):
            # (size of the element) - (delta values of each side of element),
            h = (xel[k+1]-xel[k])-(2*delta)
            # spacing between each interior node (based on polynomial order of that element)
            hi = h/(self.order[k]-1)
            for j in range(0, self.order[k]):
                if j == 0:
                    self.xnod[self.nod[k, j]] = xel[k]+delta
                else:
                    self.xnod[self.nod[k, j]] = self.xnod[self.nod[k, j-1]] + hi
        k = self.nels-1
        h = b-xel[k]-(2*delta)
        hi = h/(self.order[k]-1)
        for j in range(0, self.order[k]):
            if j == 0:
                self.xnod[self.nod[k, j]] = xel[k]+delta
            else:
                self.xnod[self.nod[k, j]] = self.xnod[self.nod[k, j-1]] + hi

        try:
            self.hmax = max(hel)
        except TypeError:
            self.hmax = hel