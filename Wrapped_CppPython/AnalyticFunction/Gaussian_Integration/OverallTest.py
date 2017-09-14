import GaussInt
import math
try:
    from colorama import Fore, Style, init
    init(autoreset=True)
    Yellow = Fore.YELLOW; Red = Fore.RED; Green = Fore.GREEN; Cyan = Fore.CYAN; Magenta = Fore.MAGENTA
    StyDim = Style.DIM
except ImportError:
    print('\nYou should get colorama. It\'s pretty sweet.\n')
    Yellow = ''; Red = ''; Green = ''; Cyan = ''; Magenta = ''
    StyDim = ''

## 1D TEST
Test1D = GaussInt.GaussianIntegration()
print(Cyan+'1D Integral');print(Test1D.GaussInt_1D(0,math.pi,100))
## 2D TEST
Test2D = GaussInt.GaussianIntegration()
print(Cyan+'2D Integral');print(Test2D.GaussInt_2D(0,math.pi,100, 0,math.pi,100))
## 3D TEST
Test3D = GaussInt.GaussianIntegration()
print(Cyan+'3D Integral');print(Test3D.GaussInt_3D(0,math.pi,100, 0,math.pi,100, 0,math.pi,100))
