import sys
from numpy import array
import ODEInt

PyObj = ODEInt.ODEClass_Test()
TBnds = ODEInt.LineDouble()
TBnds = array([0.0,20.0])
dim = 1
RelTol = 1e-6
AbsTol = 0.0
mu = 0.5

MySoln = PyObj.compute_ODE(TBnds, dim, RelTol, AbsTol, mu)
