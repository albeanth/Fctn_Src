import sys
import math
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

def QP(maxord):
    '''
    set up:
        - numerical integration for stiffness matrix
        - quadrature parameters on the reference element (-1,1)
        - number of integration points nw, nodes xw, and weights w
     '''
    xw = np.zeros(maxord)
    w = np.zeros(maxord)
    if maxord == 1:  # exact for linears
        nw = 1
        xw[0] = 0.
        w[0] = 2.
    elif maxord == 2:  # exact for cubics
        nw = 2
        xw[0] = -1/math.sqrt(3.)
        xw[1] = -xw[0]
        w[0] = 1.
        w[1] = 1.
    elif maxord == 3:  # exact for polynomials of degree 5
        nw = 3;
        xw[0] = -math.sqrt(3./5.)
        xw[1] = 0.
        xw[2] = -xw[0];
        w[0] = 5./9.
        w[1] = 8./9.
        w[2] = w[0];
    elif maxord == 4:
        nw = 4;
        xw[0] = -math.sqrt(3./7. + 2./7.*math.sqrt(6./5.))
        xw[1] = -math.sqrt(3./7. - 2./7.*math.sqrt(6./5.))
        xw[2] = -xw[1]
        xw[3] = -xw[0]
        w[0] = (18.-math.sqrt(30.))/36.
        w[1] = (18.+math.sqrt(30.))/36.
        w[2] = w[1]
        w[3] = w[0]
    elif maxord == 5:
        nw = 5;
        xw[0] = -1./3.*math.sqrt(5.+2.*math.sqrt(10./7.))
        xw[1] = -1./3.*math.sqrt(5.-2.*math.sqrt(10./7.))
        xw[2] = 0
        xw[3] = -xw[1]
        xw[4] = -xw[0]
        w[0] = (322-13*math.sqrt(70))/900
        w[1] = (322+13*math.sqrt(70))/900
        w[2] = 128/225
        w[3] = w[1]
        w[4] = w[0]
    else:
        print(Red+'Error: no code implemented for the desired order')
        sys.exit()
    return(nw,xw,w)
