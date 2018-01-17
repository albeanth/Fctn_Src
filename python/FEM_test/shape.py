import sys
import numpy as np
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
    '''
    shape function on reference element (-1,1)
    n = 2: linear
    n = 3: quadratic
    '''
    y = np.zeros(n)
    dy = np.zeros(n)
    if n == 2:
        y[0] = .5*(1-x)
        y[1] = .5*(1+x)
        dy[0] = -.5
        dy[1] = .5
    elif n == 3:
        y[0] = (x**2-x)/2
        y[1] = 1-x**2
        y[2] = (x**2+x)/2
        dy[0] = x-1/2
        dy[1] = -2*x
        dy[2] = x+1/2
    elif n == 4:
        y[0] = 1/16*(-9*(x**3) + 9*(x**2) + x - 1)
        y[1] = 1/16*(27*(x**3) - 9*(x**2) - 27*x + 9)
        y[2] = 1/16*(-27*(x**3) - 9*(x**2) + 27*x + 9)
        y[3] = 1/16*(9*(x**3) + 9*(x**2) - x - 1)
        dy[0] = 1/16*(-27*(x**2) + 18*x + 1)
        dy[1] = 1/16*(81*(x**2) - 18*x - 27)
        dy[2] = 1/16*(-81*(x**2) - 18*x + 27)
        dy[3] = 1/16*(27*(x**2) + 18*x - 1)
    else:
        print(Red+'Error: Implementation for n='+str(n)+' not provided!')
        sys.exit()
    return(y,dy)
