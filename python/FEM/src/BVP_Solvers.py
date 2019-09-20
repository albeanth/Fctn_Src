# import python packages
import sys
from numpy import zeros, outer, linalg, linspace, add, subtract, dot, array
from math import sqrt, log
import matplotlib.pyplot as plt

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
    
    def General_1D(self):
        """
        Solves a two-point BVP using Galerkin finite elements
        - automatically switches between continuous or discontinuous 
          based on mesh it has access to
        """
        nw, xw, w = QP(self.maxord)

        ## set up matrix and rhs of linear system
        stiff = zeros((self.nnodes, self.nnodes))
        mass = zeros((self.nnodes, self.nnodes))
        curr = zeros((self.nnodes, self.nnodes)) # "current" matrix used in DFEM (default to zeros for CFEM)
        rhsf = zeros( self.nnodes )
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

            if (self.nnodes == self.nels*2):
                # neutron current, J, treatments for DFEM
                # current set up, only valid for linear FEs
                #  - can make higher order FEs by generalizing
                if el == 0:
                    curr = self.Curr_iplus(curr, el)
                elif el == self.nels-1:
                    curr = self.Curr_iminus(curr, el)
                else:
                    curr = self.Curr_iminus(curr, el)
                    curr = self.Curr_iplus(curr, el)

            for idx in range(0, self.order[el]):
                rhsf[self.nod[el, idx]] += f[idx]
                for idy in range(0, self.order[el]):
                    stiff[self.nod[el, idx], self.nod[el, idy]] += k[idx][idy]
                    mass[self.nod[el, idx], self.nod[el, idy]] += m[idx][idy]

        # print("Global Current Matrix")
        # for idx in range(0, self.nnodes):
        #     for idy in range(0, self.nnodes):
        #         print("{0: .3f} ".format(curr[idx, idy]), end='')
        #     print("")
        # sys.exit()
        self.soln = zeros(self.nnodes)
        mat = add(curr,add(stiff, mass))
        # SOLVE! (for coefficients of finite elements, still not the \emph{actual} solution)
        if self.selection == 2:
            # if Dirichlet boundary conditions eliminate known values from the system
            self.soln[0] = self.u(self.bounds[0])
            subtract(rhsf, dot(mat, self.soln), out=rhsf) # used for non-zero Dirichlet BCs
            self.soln[1:] = linalg.solve(mat[1:, 1:], rhsf[1:])
        elif self.selection == 3:
            self.soln[0] = self.u(self.bounds[0])
            self.soln[-1] = self.u(self.bounds[1])
            subtract(rhsf, dot(mat, self.soln), out=rhsf) # used for non-zero Dirichlet BCs
            # self.soln[1:-1] = linalg.solve(mat[1:-1, 1:-1], rhsf[1:-1])
        else:
            # if Reflecting BCs
            self.soln = linalg.solve(mat, rhsf)


    def CFEM_1D(self):
        """
        Solves a two-point BVP using continuous Galerkin finite elements
        """
        nw, xw, w = QP(self.maxord)

        ## set up matrix and rhs of linear system
        stiff = zeros((self.nnodes, self.nnodes))
        mass = zeros((self.nnodes, self.nnodes))
        rhsf = zeros( self.nnodes )
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
                rhsf[self.nod[el, idx]] += f[idx]
                for idy in range(0, self.order[el]):
                    stiff[self.nod[el, idx], self.nod[el, idy]] += k[idx][idy]
                    mass[self.nod[el, idx], self.nod[el, idy]] += m[idx][idy]

        self.soln = zeros(self.nnodes)
        mat = add(stiff, mass)
        # if Dirichlet boundary conditions eliminate known values from the system
        self.soln[0] = self.u(self.bounds[0])
        # self.soln[-1] = self.u(self.bounds[1])
        subtract(rhsf, dot(mat, self.soln), out=rhsf) # used for non-zero Dirichlet BCs
        # SOLVE! (for coefficients of finite elements, still not the \emph{actual} solution)
        self.soln[1:] = linalg.solve(mat[1:, 1:], rhsf[1:])

        # if Reflecting BCs
        # self.soln = linalg.solve(mat, rhsf)

    def DFEM_1D(self):
        """
        Solves a two-point BVP using discontinuous Galerkin finite elements
        """
        nw, xw, w = QP(self.maxord)

        ## set up matrix and rhs of linear system
        stiff = zeros((self.nnodes, self.nnodes))
        mass = zeros((self.nnodes, self.nnodes))
        curr = zeros((self.nnodes, self.nnodes))
        rhsf = zeros(self.nnodes)
        for el in range(0, self.nels):
            xL = self.xnod[self.nod[el, 0]]
            xR = self.xnod[self.nod[el, self.order[el]-1]]
            dx = (xR - xL)/2.0
            # compute element stiffness matrix and rhs (load) vector
            k = zeros((self.order[el], self.order[el]))  # element stiffness matrix
            m = zeros((self.order[el], self.order[el]))  # element stiffness matrix
            f = zeros(self.order[el])     # element load vector
            # tmp_m = zeros((self.order[el], self.order[el]))  # element mass matrix
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
                # tmp_m += outer(psi, psi) * w[l]*dx

            # uncomment to see the local (stiffness+mass) matrix and local load vector
            # print(k)
            # print(m)
            # print(f)
            # if el == 0:
            #     b = dot(tmp_m, array([self.f(xL), self.f(xR)]))
            # else:
            #     xL = self.xnod[self.nod[el-1, self.order[el-1]-1]]
            #     b = dot(tmp_m, array([self.f(xL), self.f(xR)]))
            
            # print(tmp_m)
            # print(array([self.f(xL), self.f(xR)]))
            
            # print("{0:.4e}, {1:.4e} -> [{2:.4e}, {3:.4e}]".format(xL, xR, self.f(xL), self.f(xR)))
            # sys.exit()

            # neutron current, J, treatments
            # current set up, only valid for linear FEs 
            #  - can make higher order FEs by generalizing 
            if el == 0:
                curr = self.Curr_iplus(curr, el)
            elif el == self.nels-1:
                curr = self.Curr_iminus(curr, el)
            else:  # if el > 0 and el < self.nels-1: #
                curr = self.Curr_iminus(curr, el)
                curr = self.Curr_iplus(curr, el)


            # add the computed element stiffness matrix and load vector to the global matrix and vector
            for idx in range(0,self.order[el]):
                rhsf[self.nod[el, idx]] += f[idx]
                for idy in range(0, self.order[el]):
                    stiff[self.nod[el, idx], self.nod[el, idy]] += k[idx][idy]
                    mass[self.nod[el, idx], self.nod[el, idy]] += m[idx][idy]

        # uncomment to see the global stiffness matrix and load vector
        # print("Global Stiffness Matrix")
        # for idx in range(0, self.nnodes):
        #     for idy in range(0, self.nnodes):
        #         print("{0: .3f} ".format(stiff[idx,idy]),end='')
        #     print("")
        # print("Global Mass Matrix")
        # for idx in range(0, self.nnodes):
        #     for idy in range(0, self.nnodes):
        #         print("{0: .3f} ".format(mass[idx,idy]),end='')
        #     print("")
        # print("Global Current Matrix")
        # for idx in range(0, self.nnodes):
        #     for idy in range(0, self.nnodes):
        #         print("{0: .3f} ".format(curr[idx,idy]),end='')
        #     print("")
        # print("RHSF")
        # for idx in range(0, self.nnodes):
        #     print("{0: .3f} {1: .3f}".format(self.xnod[idx],rhsf[idx]))
        # for idx in range(1, self.nnodes-1, 2):
        #     rhsf[idx] = (rhsf[idx]+rhsf[idx+1])/2.0
        #     rhsf[idx+1] = rhsf[idx]
        # for idx in range(0, self.nnodes):
        #     print("{0: .3f} {1: .3f}".format(self.xnod[idx], rhsf[idx]))
        # sys.exit()


        self.soln = zeros(self.nnodes)
        mat = add(curr,add(stiff,mass))
        # print("Global Combined Matrix")
        # for idx in range(0, self.nnodes):
        #     for idy in range(0, self.nnodes):
        #         print("{0: .3f} ".format(mat[idx, idy]), end='')
        #     print("")
        # sys.exit()

        # self.soln[0] = self.u(self.bounds[0])
        # self.soln[-1] = self.u(self.bounds[1])
        # subtract(rhsf, dot(mat, self.soln), out=rhsf) # used for non-zero Dirichlet BCs
        # SOLVE! (for coefficients of finite elements, still not the \emph{actual} solution)
        # self.soln[1:] = linalg.solve(mat[1:, 1:], rhsf[1:])

        self.soln = linalg.solve(mat, rhsf)

    def mass_lump(self, m):
        [row, col] = m.shape
        for i in range(0,row):
            tmp = 0.0
            for j in range(0,col):
                tmp += m[i][j]
            m[i][i] = tmp
        return m
    
    def Curr_iplus(self, curr, el):
        """
        current treatment for J_{i+1/2}
        """
        # el + 1
        xL = self.xnod[self.nod[el+1, 0]]
        xR = self.xnod[self.nod[el+1, 1]]
        dx = xR-xL
        curr[self.nod[el, 1], self.nod[el+1, 0]] = -( 1.0/4.0 - self.D(xL)/(2.0*dx) )# \phi_{i+1,L}
        curr[self.nod[el, 1], self.nod[el+1, 1]] = -( self.D(xL)/(2.0*dx)           )# \phi_{i+1,R}
        # el
        xL = self.xnod[self.nod[el, 0]]
        xR = self.xnod[self.nod[el, 1]]
        dx = xR-xL
        curr[self.nod[el, 1], self.nod[el, 0]] =  self.D(xR)/(2.0*dx)           # \phi_{i,L}
        curr[self.nod[el, 1], self.nod[el, 1]] =  1.0/4.0 - self.D(xR)/(2.0*dx) # \phi_{i,R}
        return curr

    def Curr_iminus(self, curr, el):
        """
        current treatment for J_{i-1/2}
        """
        # el - 1
        xL = self.xnod[self.nod[el-1, 0]]
        xR = self.xnod[self.nod[el-1, 1]]
        dx = xR-xL
        curr[self.nod[el, 0], self.nod[el-1, 0]] = -( self.D(xR)/(2.0*dx)           )# \phi_{i-1,L}
        curr[self.nod[el, 0], self.nod[el-1, 1]] = -( 1.0/4.0 - self.D(xR)/(2.0*dx) )# \phi_{i-1,R}
        # el
        xL = self.xnod[self.nod[el, 0]]
        xR = self.xnod[self.nod[el, 1]]
        dx = xR-xL
        curr[self.nod[el, 0], self.nod[el, 0]] = 1.0/4.0 - self.D(xL)/(2.0*dx) # \phi_{i,L}
        curr[self.nod[el, 0], self.nod[el, 1]] = self.D(xL)/(2.0*dx)           # \phi_{i,R}
        return curr

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

    def Spatial_Convergence(self, L2Error, H1Error, h, flag=False):
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
                conv_l2 = ( log(L2Error[i-1]) - log(L2Error[i]) ) / ( log(h[i-1]) - log(h[i]) )
                conv_h1 = ( log(H1Error[i-1]) - log(H1Error[i]) ) / ( log(h[i-1]) - log(h[i]) )
                print("{0:.3e}  {1:.3f}    {2:.3f}".format(h[i+1], conv_l2, conv_h1))
            if flag == True:
                plt.figure(2)
                # plt.subplot(121)
                plt.loglog(h, h**2, 'x-', linewidth=2, label='O(-2)')
                plt.loglog(h, L2Error, 'x-', linewidth=2, label='L2Error')
                plt.loglog(h, h, 'x-', linewidth=2, label='O(-1)')
                plt.loglog(h, H1Error, 'x-', linewidth=2, label='H1Error')
                plt.legend(loc='upper left')
                plt.grid(True, which='both', axis='y', alpha=0.5)
                # plt.subplot(122)
                # plt.legend(loc='upper left')
                # plt.grid(True, which='both', axis='y', alpha=0.5)
                plt.show()



