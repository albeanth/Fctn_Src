import GaussInt2D
import math

myGaussian = GaussInt2D.GaussianIntegration()
myGaussian.SetupGrid(0,math.pi,100,0,math.pi,100)
test_integral = myGaussian.GaussInt_2D()

print(test_integral)
