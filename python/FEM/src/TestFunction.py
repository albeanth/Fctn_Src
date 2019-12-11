from sys import exit
from numpy import cos, sin, pi
class func:
    """
    defines the function, derivative, and source required for 
    a method of manufactured solutions (MMS) problem set up to
    a simple one-group neutron diffusion equation type problem
    - (d/dx * D(x) du/dx) + Sigma_a(x) u(x)= f(x)
    """
    def __init__(self, problem, material):
        """
        problem -> problem selection flag
        material -> homogeneous/heterogeneous/ringing
        """
        self.selection = problem
        self.hetero = material
        self.display_selection()

    def display_selection(self):
        """
        print selection to screen
        """
        if self.selection == 0:
            print("\nLinear Solution test, Dirichlet BCs left and right, ", end="")
        elif self.selection == 1:
            print("\nQuadratic Solution test, Dirichlet BCs left and right, ", end="")
        elif self.selection == 2:
            print("\nCubic Solution test, Dirichlet BCs left and right, ", end="")
        elif self.selection == 3:
            print("\nQuintic Solution test, Dirichlet BCs left and right, ", end="")
        elif self.selection == 4:
            print("\nCosine Solution test, Reflecting BCs ", end="")
        elif self.selection == 5:
            print("\nPiecewise Solution test, Dirichlet BCs left and right, ", end="")
        elif self.selection == 6:
            print("\nDirichlet left, reflecting right ",end="")
        elif self.selection == 7:
            print("\nNumerical ringing test ",end="\n")
        else:
            print("\nunknown problem selection!\n")
            exit(-1)
        
        if self.hetero == False:
            print("with homogeneous cross sections.",end="\n\n")
        elif self.hetero == True:
            print("with smooth heterogeneous cross sections.",end="\n\n")
        else:
            print("\nunknown heterogeneity selection!\n")
            exit(-1)

    def u(self, x):
        """
        base function definition
        """
        if self.selection == 0:
            val = 1.0 - x
        elif self.selection == 1:
            val = x**2 - x + 1.0
        elif self.selection == 2:
            val = (x-1)**3 + 1.0
        elif self.selection == 3:
            val = x**5 - 5.0*x**4 + 5*x**3 + 5*x**2 - 6*x + 5
        elif self.selection == 4:
            val = cos(x)+1.0
        elif self.selection == 5:
            if x < 1:
                val = x**2 + 1.0 # x + 1.0
            else:
                val = -(x-1.0)**3 + 1.0 # -x + 2.0
        elif self.selection == 6:
            val = cos(x-pi/2)+1.0
        return val

    def up(self, x):
        """
        derivative of base function
        """
        if self.selection == 0:
            val = -1.0
        elif self.selection == 1:
            val = 2*x - 1.0
        elif self.selection == 2:
            val = 3*(x-1.0)**2
        elif self.selection == 3:
            val = 5*x**4 - 20*x**3 + 15*x**2 + 10*x - 6
        elif self.selection == 4:
            val = -sin(x)
        elif self.selection == 5:
            if x < 1.0:
                val = 2.0*x #1.0
            else:
                val = -3.0*(x-1.0)**2 #-1.0
        elif self.selection == 6:
            val = cos(x)
        return val

    def upp(self, x):
        """
        derivative of base function
        """
        if self.selection == 0:
            val = 0.0
        elif self.selection == 1:
            val = 2.0
        elif self.selection == 2:
            val = 6*(x-1.0)
        elif self.selection == 3:
            val = 20*x**3 - 60*x**2 + 30*x + 10
        elif self.selection == 4:
            val = -cos(x)
        elif self.selection == 5:
            if x < 1.0:
                val = 2.0  # 0.0
            else:
                val = -6*(x-1.0)  # 0.0
        elif self.selection == 6:
            val = -sin(x)
        return val

    def SigA(self, x):
        """
        absorption cross section
        """
        if self.hetero == False:
            val = 1.0 # homogeneous case
        elif self.hetero == True:
            val = x+1.0
        
        if self.selection == 5:
            if x < 1.0:
                val = 1.0E4
            else:
                val = 1.0E4
        # elif self.selection == 7:
        #     if x < 2:
        #         val = 1.0E10
        #     elif x < 4:
        #         val = 1.0E10
        #     else:
        #         val = 1.0E10
        return val

    def D(self, x):
        """
        diffusion coefficient, D(x)
        """
        val = 1.0/(3.0*self.SigA(x))
        if self.selection == 5:
            if x < 1.0:
                val = 1.0E-4
            else:
                val = 1.0E-2
        return val

    def Dx(self, x):
        """
        first spatial derivative of diffusion coefficient, D'(x)
        """
        if self.hetero == False:
            val = 0.0    # homogeneous case
        else:
            val = -1.0/(3*(x+1)**2.0)  # heterogeneous case
        return val
    
    def f(self, x):
        """
        source using T and Tp that defines the MMS problem
        """
        if self.selection < 7:
            val = -(self.Dx(x)*self.up(x) + self.D(x)*self.upp(x)) + self.SigA(x)*self.u(x)
        elif self.selection == 7:
            if x < 2:
                val = 1.0E-1
            elif x < 4:
                val = 1.0E-1
            else:
                val = 1.0E-1
        return val
