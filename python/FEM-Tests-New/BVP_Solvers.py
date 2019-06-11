# import python packages
import sys
from numpy import zeros, outer, linalg, linspace
from math import sqrt
import matplotlib.pyplot as plt

# import user defined functions and classes
from QuadParams import QP
from ShapeFuncs import shape
from FEM_1D_Mesh import mesh
from TestFunction import func

class BVP_Solvers(mesh, func):
    """
    solves a two-point boundary value problem (BVP) using continuous
    or discontinuous Galerkin finite elements

    A. Alberti 06/2017 (code structure updated 06/2019)
    Referenced from M. Peszynska MTH659 FEM Course
    """

    def __init__(self):
        """ default constructor """
    
    def CFEM_1D(self):
        """
        Solves a continuous Galerkin formulation
        INPUT:
            grid - mesh object containing 1D CFEM grid
            T    - TestFunction object containing MMS function info
        OUTPUT:
            sol - vector of FE coefficients for order of solution
        """
        nw, xw, w = QP(self.maxord)

        ## set up matrix and rhs of linear system
        stiff = zeros((self.nnodes, self.nnodes))
        rhsf = zeros( self.nnodes )
        for el in range(0, self.nels):
            xL = self.xnod[self.nod[el, 0]]
            xR = self.xnod[self.nod[el, self.order[el]-1]]
            dx = (xR - xL)/2.0
            # compute element stiffness matrix and rhs (load) vector
            k = zeros((self.order[el], self.order[el]))  # element stiffness matrix
            f = zeros(self.order[el])     # element load vector

            for l in range(0, nw):
                # x runs in true element, xw runs in reference element
                x = xL + (1.0 + xw[l])*dx
                # calculations on ref element
                psi, dpsi = shape(xw[l], self.order[el])
                kval = self.k(x)
                fval = self.f(x)
                f += fval * psi * w[l]*dx
                k += kval*outer(dpsi, dpsi)/dx/dx * w[l]*dx

            # # uncomment to see the local (stiffness+mass) matrix and local load vector
            # print(k)
            # print(f)
            # sys.exit()

            for idx in range(0, self.order[el]):
                rhsf[self.nod[el, idx]] += f[idx]
                for idy in range(0, self.order[el]):
                    stiff[self.nod[el, idx], self.nod[el, idy]] += k[idx][idy]

        # impose Dirichlet boundary conditions eliminate known values from the system
        self.soln = zeros(self.nnodes)
        self.soln[0] = self.u(self.bounds[0])
        self.soln[-1] = self.u(self.bounds[1])
        mat = stiff

        # rhsf = np.subtract(rhsf , np.dot(mat,sol)) # used for non-zero Dirichlet BCs

        # SOLVE! (for coefficients of finite elements, still not the \emph{actual} solution)
        self.soln[1:-1] = linalg.solve(mat[1:-1, 1:-1], rhsf[1:-1])

    def Error(self):
        """
        computes the error from an MMS type BVP problem
        """
        nw, xw, w = QP(self.maxord)

        self.l2Err = 0.0
        self.h1Err = 0.0
        for el in range(0, self.nels):
            xL = self.xnod[self.nod[el, 0]]
            xR = self.xnod[self.nod[el, self.order[el]-1]]
            dx = (xR - xL)/2.0
            for l in range(0, nw):
                x = xL + (1.0 + xw[l])*dx
                psi, dpsi = shape(xw[l], self.order[el])
                uval = self.u(x)
                duval = self.up(x)
                uhval = 0.0
                duhval = 0.0
                for k in range(0, self.order[el]):
                    mynum = self.nod[el, k]
                    uhval += self.soln[mynum]*psi[k]
                    duhval += self.soln[mynum]*dpsi[k]/dx
                
                # complete integration over element
                self.l2Err += (uval-uhval)**2 *w[l]*dx
                self.h1Err += (duval-duhval)**2 *w[l]*dx

        self.l2Err = sqrt(self.l2Err)
        self.h1Err = sqrt(self.l2Err + self.h1Err)

    def Plot(self, NumPts=10):
        """
        plots CFEM solution to BVP (red) with analytic MMS solution (blue)
        """
        plt.figure(1)
        xplot = linspace(self.bounds[0], self.bounds[1], self.nels*NumPts+1)
        plt.plot(xplot, self.u(xplot), linewidth = 2, color = 'blue')

        nw, xw, w = QP(self.maxord)
        for el in range(0, self.nels):
            xL = self.xnod[self.nod[el, 0]]
            xR = self.xnod[self.nod[el, self.order[el]-1]]
            dx = (xR - xL)/2.0
            xpp = linspace(xL, xR, NumPts)
            xww = (xpp-xL)/dx - 1.0
            ypp = zeros(NumPts)
            for j in range(0, NumPts):
                psi, dpsi = shape(xww[j], self.order[el])
                uhval = 0.0
                for k in range(0, self.order[el]):
                    mynum = self.nod[el, k]
                    uhval += self.soln[mynum]*psi[k]
                ypp[j] = uhval

            plt.plot(xpp,ypp, linewidth=2, color = 'red')
        
        plt.title('Exact solution (blue) and BVP FEM solution (red)')
        plt.xlabel('Space')
        plt.grid(True)
        plt.show()



