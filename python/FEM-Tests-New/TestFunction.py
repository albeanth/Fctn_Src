from numpy import cos, sin
class func:
    """
    defines the function, derivative, and source required for 
    a method of manufactured solutions (MMS) problem set up to
    a simple heat equation type problem
    - (d/dx * K(x) du/dx) = f(x)
    """
    def __init__(self):
        """standard constructor"""

    def u(self, x):
        """
        base function definition
        """
        # return x**2.0
        return sin(x)

    def up(self, x):
        """
        derivative of base function
        """
        # return 2.0*x
        return cos(x)

    def upp(self, x):
        """
        derivative of base function
        """
        # return 2.0
        return -sin(x)

    def k(self, x):
        """
        thermal conductivity, k(x)
        """
        return x

    def kp(self, x):
        """
        thermal conductivity, K(x)
        """
        return 1.0
    
    def f(self, x):
        """
        source using T and Tp that defines the MMS problem
        """
        return -(self.kp(x)*self.up(x) + self.k(x)*self.upp(x))
    
