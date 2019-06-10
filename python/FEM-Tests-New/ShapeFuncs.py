import sys
from numpy import zeros
try:
    from colorama import Fore, Style, init
    init(autoreset=True)
    Yellow = Fore.YELLOW; Red = Fore.RED; Green = Fore.GREEN; Cyan = Fore.CYAN; Magenta = Fore.MAGENTA
    StyDim = Style.DIM
except ImportError:
    print('\nYou should get colorama. It\'s pretty sweet.\n')
    Yellow = ''; Red = ''; Green = ''; Cyan = ''; Magenta = ''
    StyDim = ''


def shape(x,n):
    """
    evalutes shape function defined 
    on reference element (-1,1) at x
    OUTPUT:
      y  = base function
      dy = first derivative
    """
    y = zeros(n)
    dy = zeros(n)
    if n == 2: # linear
        y[0] = 0.5*(1.0 - x)
        y[1] = 0.5*(1.0 + x)
        dy[0] = -0.5
        dy[1] =  0.5
    elif n == 3: # quadratic
        y[0] = (x**2 - x) / 2.0
        y[1] = 1.0 - x**2.0
        y[2] = (x**2 + x)/2.0
        dy[0] = x - 1.0/2.0
        dy[1] = -2.0 * x
        dy[2] = x + 1.0/2.0
    elif n == 4: # cubic
        y[0] =  1.0/16.0 * (-9.0 *(x**3) + 9.0*(x**2) + x - 1.0)
        y[1] =  1.0/16.0 * ( 27.0*(x**3) - 9.0*(x**2) - 27.0*x + 9.0)
        y[2] =  1.0/16.0 * (-27.0*(x**3) - 9.0*(x**2) + 27.0*x + 9.0)
        y[3] =  1.0/16.0 * ( 9.0 *(x**3) + 9.0*(x**2) - x - 1.0)
        dy[0] = 1.0/16.0 * (-27.0*(x**2) + 18.0*x + 1.0)
        dy[1] = 1.0/16.0 * ( 81.0*(x**2) - 18.0*x - 27.0)
        dy[2] = 1.0/16.0 * (-81.0*(x**2) - 18.0*x + 27.0)
        dy[3] = 1.0/16.0 * ( 27.0*(x**2) + 18.0*x - 1.0)
    else:
        print(Red+'Error: Implementation for n='+str(n)+' not provided!')
        sys.exit()
    
    return(y,dy)
