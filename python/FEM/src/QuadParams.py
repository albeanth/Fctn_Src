import sys
from math import sqrt
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

def QP(maxord):
    """
    quadrature set for a given order defined on the reference element [-1,1]
    OUTPUT:
        nw = number of integration points
        xw = nodes locations
        w = weights
    """
    xw = zeros(maxord)
    w = zeros(maxord)
    if maxord == 1:  # exact for linears
        nw = 1
        xw[0] = 0.0
        w[0] = 2.0
    elif maxord == 2:  # exact for cubics
        nw = 2
        xw[0] = -1.0/sqrt(3.0)
        xw[1] = -xw[0]
        w[0] = 1.0
        w[1] = 1.0
    elif maxord == 3:  # exact for polynomials of degree 5
        nw = 3
        xw[0] = -sqrt(3.0/5.0)
        xw[1] = 0.0
        xw[2] = -xw[0]
        w[0] = 5.0/9.0
        w[1] = 8.0/9.0
        w[2] = w[0]
    elif maxord == 4:
        nw = 4
        xw[0] = -sqrt(3.0/7.0 + 2.0/7.0 * sqrt(6.0/5.0))
        xw[1] = -sqrt(3.0/7.0 - 2.0/7.0 * sqrt(6.0/5.0))
        xw[2] = -xw[1]
        xw[3] = -xw[0]
        w[0] = (18.0 - sqrt(30.0)) / 36.0
        w[1] = (18.0 + sqrt(30.0)) / 36.0
        w[2] = w[1]
        w[3] = w[0]
    elif maxord == 5:
        nw = 5
        xw[0] = -1.0/3.0 *sqrt(5.0 + 2.0 * sqrt(10.0/7.0))
        xw[1] = -1.0/3.0 *sqrt(5.0 - 2.0 * sqrt(10.0/7.0))
        xw[2] = 0.0
        xw[3] = -xw[1]
        xw[4] = -xw[0]
        w[0] = (322.0 - 13.0 * sqrt(70)) / 900.0
        w[1] = (322.0 + 13.0 * sqrt(70)) / 900.0
        w[2] = 128.0/225.0
        w[3] = w[1]
        w[4] = w[0]
    else:
        print(Red+'Error: no code implemented for the desired order')
        sys.exit()
    
    return(nw,xw,w)
