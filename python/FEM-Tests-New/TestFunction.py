from numpy import cos, sin
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
        return cos(x)+1.0

    def up(self, x):
        """
        derivative of base function
        """
        return -sin(x)

    def upp(self, x):
        """
        derivative of base function
        """
        # return 2.0
        return -cos(x)

    def SigA(self, x):
        """
        absorption cross section
        """
        return x+1.0

    def D(self, x):
        """
        diffusion coefficient, D(x)
        """
        return 1.0/(3.0*self.SigA(x))

    def Dx(self, x):
        """
        first spatial derivative of diffusion coefficient, D'(x)
        """
        return -1.0/(3*(x+1)**2.0)
    
    def f(self, x):
        """
        source using T and Tp that defines the MMS problem
        """
        return -(self.Dx(x)*self.up(x) + self.D(x)*self.upp(x)) + self.SigA(x)*self.u(x)