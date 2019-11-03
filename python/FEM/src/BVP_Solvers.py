# import python packages
import sys
from numpy import zeros, outer, linalg, linspace, add, subtract, dot, array
from math import sqrt, log
import matplotlib.pyplot as plt
try:
    from colorama import Fore, Style, init
    init(autoreset=True)
    Yellow = Fore.YELLOW; Red = Fore.RED; Green = Fore.GREEN; Cyan = Fore.CYAN; Magenta = Fore.MAGENTA
    StyDim = Style.DIM
except ImportError:
    print('\nYou should get colorama. It\'s pretty sweet.\n')
    Yellow = ''; Red = ''; Green = ''; Cyan = ''; Magenta = ''
    StyDim = ''

# import user defined functions and classes
from src.QuadParams import QP
from src.ShapeFuncs import shape
from src.FEM_1D_Mesh import mesh
from src.TestFunction import func

class BVP_Solvers(mesh, func):
    """
    solves a two-point boundary value problem (BVP) using continuous
    or discontinuous Galerkin finite elements

    A. Alberti 06/2017 (code structure updated 06/2019)
    Referenced from M. Peszynska MTH659 FEM Course
    """

    def __init__(self, choice, material):
        func.__init__(self, choice, material)
    
    def CFEM_1D(self):
        """
        Solves a two-point BVP using Galerkin finite elements
        - automatically switches between continuous or discontinuous 
          based on mesh it has access to
        """
        nw, xw, w = QP(self.maxord)

        ## set up matrix and rhs of linear system
        self.stiff = zeros((self.nnodes, self.nnodes))
        self.mass = zeros((self.nnodes, self.nnodes))
        self.rhsf = zeros( self.nnodes )
        for el in range(0, self.nels):
            xL = self.xnod[self.nod[el, 0]]
            xR = self.xnod[self.nod[el, self.order[el]-1]]
            dx = (xR - xL)/2.0
            # compute element stiffness matrix and rhs (load) vector
            k = zeros((self.order[el], self.order[el]))  # element stiffness matrix
            m = zeros((self.order[el], self.order[el]))  # element mass matrix
            f = zeros(self.order[el])     # element load vector

            for l in range(0, nw):
                # x runs in true element, xw runs in reference element
                x = xL + (1.0 + xw[l])*dx
                # calculations on ref element
                psi, dpsi = shape(xw[l], self.order[el])
                Dval = self.D(x)
                SigAbs = self.SigA(x)
                fval = self.f(x)
                f += fval * psi * w[l]*dx
                k += Dval*outer(dpsi, dpsi)/dx/dx * w[l]*dx
                m += SigAbs*outer(psi, psi) * w[l]*dx

            # # uncomment to see the local (stiffness+mass) matrix and local load vector
            # print(k)
            # print(m)
            # print(f)
            # sys.exit()

            for idx in range(0, self.order[el]):
                self.rhsf[self.nod[el, idx]] += f[idx]
                for idy in range(0, self.order[el]):
                    self.stiff[self.nod[el, idx], self.nod[el, idy]] += k[idx][idy]
                    self.mass[self.nod[el, idx], self.nod[el, idy]] += m[idx][idy]

        self.soln = zeros(self.nnodes)
        mat = add(self.stiff, self.mass)
        # SOLVE! (for coefficients of finite elements, still not the \emph{actual} solution)
        if ((self.selection == 0) or (self.selection == 1) or (self.selection == 2) or (self.selection == 3) or (self.selection == 5)):
            # if Dirichlet boundary conditions eliminate known values from the system
            self.soln[0] = self.u(self.bounds[0])
            self.soln[-1] = self.u(self.bounds[-1])
            subtract(self.rhsf, dot(mat, self.soln), out=self.rhsf) # used for non-zero Dirichlet BCs
            self.soln[1:-1] = linalg.solve(mat[1:-1, 1:-1], self.rhsf[1:-1])
        elif self.selection == 6:
            # if Dirichlet boundary conditions eliminate known values from the system
            self.soln[0] = self.u(self.bounds[0])
            subtract(self.rhsf, dot(mat, self.soln), out=self.rhsf) # used for non-zero Dirichlet BCs
            self.soln[1:] = linalg.solve(mat[1:, 1:], self.rhsf[1:])
        else:
            # if Reflecting BCs
            self.soln = linalg.solve(mat, self.rhsf)

    def DFEM_1D(self):
        """
        Solves a two-point BVP using Galerkin finite elements
        - automatically switches between continuous or discontinuous 
          based on mesh it has access to
        """
        nw, xw, w = QP(self.maxord)

        ## set up matrix and rhs of linear system
        self.stiff = zeros((self.nnodes, self.nnodes))
        self.mass = zeros((self.nnodes, self.nnodes))
        self.curr = zeros((self.nnodes, self.nnodes)) # "current" matrix used in DFEM (default to zeros for CFEM)
        self.rhsf = zeros( self.nnodes )
        for el in range(0, self.nels):
            xL = self.xnod[self.nod[el, 0]]
            xR = self.xnod[self.nod[el, self.order[el]-1]]
            dx = (xR - xL)/2.0
            # compute element stiffness matrix and rhs (load) vector
            k = zeros((self.order[el], self.order[el]))  # element stiffness matrix
            m = zeros((self.order[el], self.order[el]))  # element mass matrix
            f = zeros(self.order[el])     # element load vector

            for l in range(0, nw):
                # x runs in true element, xw runs in reference element
                x = xL + (1.0 + xw[l])*dx
                # calculations on ref element
                psi, dpsi = shape(xw[l], self.order[el])
                Dval = self.D(x)
                SigAbs = self.SigA(x)
                fval = self.f(x)
                f += fval * psi * w[l]*dx
                k += Dval*outer(dpsi, dpsi)/dx/dx * w[l]*dx
                m += SigAbs*outer(psi, psi) * w[l]*dx

            # # uncomment to see the local (stiffness+mass) matrix and local load vector
            # print(k)
            # print(m)
            # print(f)
            # sys.exit()

            # neutron current, J, treatments for DFEM
            # current set up, only valid for linear FEs 
            #  - can make higher order FEs by generalizing 
            if el == 0:
                self.Curr_iplus(el)
            elif el == self.nels-1:
                self.Curr_iminus(el)
            else:
                self.Curr_iminus(el)
                self.Curr_iplus(el)

            for idx in range(0, self.order[el]):
                self.rhsf[self.nod[el, idx]] += f[idx]
                for idy in range(0, self.order[el]):
                    self.stiff[self.nod[el, idx], self.nod[el, idy]] += k[idx][idy]
                    self.mass[self.nod[el, idx], self.nod[el, idy]] += m[idx][idy]

        # print("Writing Out Global Matrices")
        # cmat = open('CurrentMatrix.csv','w'); print("  current matrix")
        # mmat = open('MassMatrix.csv','w'); print("  mass matrix")
        # kmat = open('StiffMatrix.csv','w'); print("  stiffness matrix")
        # for idx in range(0, self.nnodes):
        #     for idy in range(0, self.nnodes):
        #         cmat.write("{0:.5e},".format(self.curr[idx, idy]))
        #         mmat.write("{0:.5e},".format(self.mass[idx, idy]))
        #         kmat.write("{0:.5e},".format(self.stiff[idx, idy]))
        #     cmat.write("\n")
        #     mmat.write("\n")
        #     kmat.write("\n")
        # cmat.close()
        # mmat.close()
        # kmat.close()
        # print("Global RHS")
        # for i in range(0,self.nnodes):
        #     print("{0:.3e} {1:.6e}".format(self.xnod[i], self.rhsf[i]))
        # sys.exit()
        self.soln = zeros(self.nnodes)
        mat = add(self.curr, add(self.stiff, self.mass))
        # print("Writing Out Global Matrix")
        # gmat = open('GlobalMatrix.csv','w'); print("  global matrix")
        # for idx in range(0, self.nnodes):
        #     for idy in range(0, self.nnodes):
        #         gmat.write("{0:.5e},".format(mat[idx, idy]))
        #     gmat.write("\n")
        # gmat.close()
        # sys.exit()
        if ((self.selection == 0) or (self.selection == 1) or (self.selection == 2) or (self.selection == 3) or (self.selection == 5)):
            self.soln[0] = self.u(self.bounds[0])
            self.soln[-1] = self.u(self.bounds[-1])
            subtract(self.rhsf, dot(mat, self.soln), out=self.rhsf) # used for non-zero Dirichlet BCs
            self.soln[1:-1] = linalg.solve(mat[1:-1, 1:-1], self.rhsf[1:-1])
        elif self.selection == 6:
            # if Dirichlet boundary conditions eliminate known values from the system
            self.soln[0] = self.u(self.bounds[0])
            subtract(self.rhsf, dot(mat, self.soln), out=self.rhsf) # used for non-zero Dirichlet BCs
            self.soln[1:] = linalg.solve(mat[1:, 1:], self.rhsf[1:])
        else:
        self.soln = linalg.solve(mat, self.rhsf)

    def mass_lump(self, m):
        [row, col] = m.shape
        for i in range(0,row):
            tmp = 0.0
            for j in range(0,col):
                tmp += m[i][j]
            m[i][i] = tmp
        return m
    
    def Curr_iplus(self, el):
        """
        current treatment for J_{i+1/2}
        """
        # el + 1
        xL = self.xnod[self.nod[el+1, 0]]
        xR = self.xnod[self.nod[el+1, 1]]
        dx = xR-xL
        self.curr[self.nod[el, 1], self.nod[el+1, 0]] = -1.0/4.0 + self.D(xL)/(2.0*dx) # \phi_{i+1,L} 
        self.curr[self.nod[el, 1], self.nod[el+1, 1]] = -self.D(xL)/(2.0*dx)           # \phi_{i+1,R}
        # el
        xL = self.xnod[self.nod[el, 0]]
        xR = self.xnod[self.nod[el, 1]]
        dx = xR-xL
        self.curr[self.nod[el, 1], self.nod[el, 0]] =  self.D(xR)/(2.0*dx)           # \phi_{i,L}
        self.curr[self.nod[el, 1], self.nod[el, 1]] =  1.0/4.0 - self.D(xR)/(2.0*dx) # \phi_{i,R}

    def Curr_iminus(self, el):
        """
        current treatment for J_{i-1/2}
        """
        # el - 1
        xL = self.xnod[self.nod[el-1, 0]]
        xR = self.xnod[self.nod[el-1, 1]]
        dx = xR-xL
        self.curr[self.nod[el, 0], self.nod[el-1, 0]] = -self.D(xR)/(2.0*dx)           # \phi_{i-1,L}
        self.curr[self.nod[el, 0], self.nod[el-1, 1]] = -1.0/4.0 + self.D(xR)/(2.0*dx) # \phi_{i-1,R}
        # el
        xL = self.xnod[self.nod[el, 0]]
        xR = self.xnod[self.nod[el, 1]]
        dx = xR-xL
        self.curr[self.nod[el, 0], self.nod[el, 0]] = 1.0/4.0 - self.D(xL)/(2.0*dx) # \phi_{i,L}
        self.curr[self.nod[el, 0], self.nod[el, 1]] = self.D(xL)/(2.0*dx)           # \phi_{i,R}

    def L2Error(self):
        """
        computes the L2 error from an MMS type BVP problem
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

    def LinfError(self):
        """
        computes the L_infty error from an MMS type BVP problem
        """
        nw, xw, w = QP(self.maxord)
        self.linf_ErrVal = 0.0
        self.linf_Loc = None
        ## set mesh parameters to spatial grid
        for el in range(0, self.nels):  # for each element in all elements
            xL = self.xnod[self.nod[el, 0]]
            xR = self.xnod[self.nod[el, self.order[el]-1]]
            dx = (xR - xL)/2.0
            for l in range(0, nw):
                x = xL + (1.0 + xw[l])*dx
                psi, dpsi = shape(xw[l], self.order[el])
                uval = self.u(x)
                uhval = 0.0
                for k in range(0, self.order[el]):
                    mynum = self.nod[el, k]
                    uhval += self.soln[mynum]*psi[k]
                if abs(uval-uhval) > self.linf_ErrVal:
                    self.linf_ErrVal = abs(uval-uhval)
                    self.linf_Loc = x
        
        if self.linf_Loc == None:
            print("Max error not found! That doesn't make sense...")
            sys.exit()

    
    def Plot(self, NumPts=10, string="CFEM"):
        """
        plots CFEM solution to BVP (red) with analytic MMS solution (blue)
        """
        plt.figure(1)
        # xplot = linspace(self.bounds[0], self.bounds[1], 100)
        # plt.plot(xplot, self.u(xplot), linewidth = 2, color = 'blue')

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

            # if string == "DFEM":    
            #     plt.plot(xpp,ypp, ".-",linewidth=1, color = 'blue')
            # else:
            #     plt.plot(xpp,ypp, ".-",linewidth=1, color = 'red')
            if self.nels == 16:
                plt.plot(xpp, ypp, ".-", linewidth=1, color='red')
            elif self.nels == 32:
                plt.plot(xpp,ypp, ".-",linewidth=1, color = 'blue')
            elif self.nels == 64:
                plt.plot(xpp,ypp, ".-",linewidth=1, color = 'green')
        
        # plt.title("Exact solution (blue) and BVP {0} solution (red)".format(string))
        plt.xlabel('Space')
        plt.grid(True)
        # if string == "DFEM":
            # plt.title("DFEM (blue) and CFEM (red)")
        #     # plt.title("CFEM (red)")
        if self.nels == 64:
            plt.title("#Elem = 16 (red) and #Elem = 64 (green)")
            plt.show()

    def Spatial_Convergence(self, L2Error, H1Error, LocalError, h, flag=False):
        """
        computes the convergence of the L2/H1Error as a function of mesh refinement, h
        INPUT:
            L2Error = vector, computed as a function of h
            H1Error = vector, computed as a function of h
            LocalError = vector, computed as a function of h
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
            print("dx\t   L2Conv   H1Conv   LocalConv")
            print("{0:.3e}  --\t    --".format(h[0]))
            for i in range(0,len(L2Error)-1):
                try:
                    conv_l2 = ( log(L2Error[i]) - log(L2Error[i+1]) ) / ( log(h[i]) - log(h[i+1]) )
                except ValueError:
                    conv_l2 == 0.0
                try:
                    conv_h1 = ( log(H1Error[i]) - log(H1Error[i+1]) ) / ( log(h[i]) - log(h[i+1]) )
                except ValueError:
                    conv_h1 = 0.0
                try:
                    conv_local = ( log(LocalError[i]) - log(LocalError[i+1]) ) / ( log(h[i]) - log(h[i+1]) )
                except ValueError:
                    conv_local = 0.0
                print("{0:.3e}  {1:.3f}    {2:.3f}    {3:.3f}".format(h[i+1], conv_l2, conv_h1, conv_local))
            if flag == True:
                plt.figure(2)
                # plt.subplot(121)
                plt.loglog(h, h**2/10, 'x-', linewidth=2, label='O(-2)', alpha=0.35)
                plt.loglog(h, L2Error, 'x-', linewidth=2, label='L2Error')
                plt.loglog(h, LocalError, 'x-', linewidth=2, label='LInfError')
                plt.loglog(h, h**3/10, 'x-', linewidth=2, label='O(-3)', alpha=0.35)
                plt.xlabel('dx')
                plt.legend(loc='best')
                plt.grid(True, which='both', axis='y', alpha=0.5)
                plt.show()



