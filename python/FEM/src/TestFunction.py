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
        if self.selection == 1:
            print("\nReflecting BCs ", end="")
        elif self.selection == 2:
            print("\nDirichlet left, reflecting right ",end="")
        elif self.selection == 3:
            print("\nNumerical ringing test ",end="")
        if self.hetero == 1:
            print("with homogeneous cross sections.",end="\n\n")
        elif self.hetero == 2:
            print("with smooth heterogeneous cross sections.",end="\n\n")

    def u(self, x):
        """
        base function definition
        """
        if self.selection == 1:
            val = cos(x)+1.0
        elif self.selection == 2:
            val = cos(x-pi/2)+1.0
        elif self.selection == 3:
            val = 1/(1+x**2)
        else:
            print("def u(self,x): unknown problem selection!")
            exit(-1)
        return val

    def up(self, x):
        """
        derivative of base function
        """
        if self.selection == 1:
            val = -sin(x)
        elif self.selection == 2:
            val = cos(x)
        elif self.selection == 3:
            val = (-2*x)/((1+x**2)**2)
        else:
            print("def up(self, x): unknown problem selection!")
            exit(-1)
        return val

    def upp(self, x):
        """
        derivative of base function
        """
        if self.selection == 1:
            val = -cos(x)
        elif self.selection == 2:
            val = -sin(x)
        elif self.selection == 3:
            val = 8*x**2/(1+x**2)**3 - 2/(1+x**2)**2
        else:
            print("def upp(self, x): unknown problem selection!")
            exit(-1)
        return val

    def SigA(self, x):
        """
        absorption cross section
        """
        if self.hetero == 1:
            val = 1.0 # homogeneous case
        elif self.hetero == 2:
            val = x+1.0  # heterogeneous case
        elif self.hetero == 3:
            if x < 2:
                val = 1E-4
            elif x < 4:
                val = 100.0
            else:
                val = 2.0
        else:
            print("def SigA(self, x): unknown heterogeneity selection!")
            exit(-1)
        return val

    def D(self, x):
        """
        diffusion coefficient, D(x)
        """
        if ((self.hetero == 1) or (self.hetero == 2)):
            val = 1.0/(3.0*self.SigA(x)) # heterogeneous case
        elif self.hetero == 3:
            if x < 2:
                val = 4
            elif x < 4:
                val = 1.0
            else:
                val = 1.0
        else:
            print("def D(self, x): unknown material selection!")
            exit(-1)
        return val

    def Dx(self, x):
        """
        first spatial derivative of diffusion coefficient, D'(x)
        """
        if self.hetero == 1:
            val = 0.0    # homogeneous case
        else:
            val = -1.0/(3*(x+1)**2.0)  # heterogeneous case
        return val
    
    def f(self, x):
        """
        source using T and Tp that defines the MMS problem
        """
        if ((self.selection == 1) or (self.selection == 2)):
            val = -(self.Dx(x)*self.up(x) + self.D(x)*self.upp(x)) + self.SigA(x)*self.u(x)
        elif self.selection == 3:
            if x < 2:
                val = 50.0
            elif x < 4:
                val = 50.0
            else:
                val = 0.0
            return val
        else:
            print("def f(self, x): unknown problem selection!")
            exit(-1)
        return val
