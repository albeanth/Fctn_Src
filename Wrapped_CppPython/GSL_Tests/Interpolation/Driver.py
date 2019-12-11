import sys
from numpy import array
from math import pi
import InterpTest

class Test():
    def __init__(self, size, method):
        self.size = size
        self.method = method

PyObj = InterpTest.InterpClass_Test()
TBnds = InterpTest.LineDouble()
TBnds = array([0,pi])
size = 10

method = "cspline"
dummy = Test(size,method)
ScatterData = PyObj.set_xy(size, TBnds)
MySpline = PyObj.get_spline(dummy)
MyData = PyObj.eval_spline(MySpline, ScatterData, size, 1)
IntVal = PyObj.integrate_spline(MySpline, ScatterData, size)

# print('#Sample Data')
# for i in range(0,len(ScatterData.x)):
#     print('{0:.4f}, {1:.6f}'.format(ScatterData.x[i],ScatterData.y[i]))
print('\n\n#Interp. Data')
for i in range(0,len(MyData.x)):
    print('{0:.4f}, {1:.6f}, {2:.6f}'.format(MyData.x[i],MyData.y[i],MyData.y_p[i]))
print("\n\nIntVal = "+str(IntVal))
