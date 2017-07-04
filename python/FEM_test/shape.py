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
        y[0] = (np.subtract(np.square(x),x))/2
        y[1] = 1-np.square(x)
        y[2] = (np.add(np.square(x),x))/2
        dy[0] = x-1/2
        dy[1] = -2*x
        dy[2] = x+1/2
    else:
        print(Red+'Error: Implementation for n='+str(n)+' not provided!')
        sys.exit()
    return(y,dy)

# def shape(x,n):
#     '''
#     shape function on reference element (-1,1)
#     n = 2: linear
#     n = 3: quadratic
#     '''
#     try:
#         y = np.zeros((n,len(x)))
#         dy = np.zeros((n,len(x)))
#     except TypeError:
#         y = np.zeros((n,1))
#         dy = np.zeros((n,1))
#     if n == 2:
#         y[0,:] = .5*(1-x)
#         y[1,:] = .5*(1+x)
#         dy[0,:] = -.5
#         dy[1,:] = .5
#     elif n == 3:
#         y[0,:] = (np.subtract(np.square(x),x))/2
#         y[1,:] = 1-np.square(x)
#         y[2,:] = (np.add(np.square(x),x))/2
#         dy[0,:] = x-1/2
#         dy[1,:] = -2*x
#         dy[2,:] = x+1/2
#     else:
#         print(Red+'Error: Implementation for n='+str(n)+' not provided!')
#         sys.exit()
#     return(y,dy)
