from numpy import cos, sin, pi
class func:
    """
    defines the function, derivative, and source required for 
    a method of manufactured solutions (MMS) problem set up to
    a simple one-group neutron diffusion equation type problem
    - (d/dx * D(x) du/dx) + Sigma_a(x) u(x)= f(x)
    """
    def __init__(self):
        """standard constructor"""

    def u(self, x):
        """
        base function definition
        """
        # return cos(x)+1.0
        return cos(x-pi/2)+1.0
        # return 1/(1+x**2)

    def up(self, x):
        """
        derivative of base function
        """
        # return -sin(x)
        return cos(x)
        # return (-2*x)/((1+x**2)**2)

    def upp(self, x):
        """
        derivative of base function
        """
        # return -cos(x)
        return -sin(x)
        # return 8*x**2/(1+x**2)**3 - 2/(1+x**2)**2

    def SigA(self, x):
        """
        absorption cross section
        """
        return x+1.0 # heterogeneous case
        # return 1.0 # homogeneous case
        # if x < 2:
        #     val = 1E-4
        # elif x < 4:
        #     val = 100.0
        # else:
        #     val = 2.0
        # return val

    def D(self, x):
        """
        diffusion coefficient, D(x)
        """
        return 1.0/(3.0*self.SigA(x))
        # if x < 2:
        #     val = 4
        # elif x < 4:
        #     val = 1.0
        # else:
        #     val = 1.0
        # return val

    def Dx(self, x):
        """
        first spatial derivative of diffusion coefficient, D'(x)
        """
        return -1.0/(3*(x+1)**2.0)  # heterogeneous case
        # return 0.0    # homogeneous case
    
    def f(self, x):
        """
        source using T and Tp that defines the MMS problem
        """
        return -(self.Dx(x)*self.up(x) + self.D(x)*self.upp(x)) + self.SigA(x)*self.u(x)
        # return 10.0
        # return x+1.0
        # if x < 2:
        #     val = 50.0
        # elif x < 4:
        #     val = 50.0
        # else:
        #     val = 0.0
        # return val