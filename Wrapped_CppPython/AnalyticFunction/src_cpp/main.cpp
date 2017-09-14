#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include "GaussianIntegration.h"
using namespace std;

int main()
{
  GaussianIntegration mytest;
  double test_integral;
  test_integral = mytest.GaussInt_1D(0,M_PI,100);
  // test_integral = mytest.GaussInt_2D(0,M_PI,100, 0,M_PI,100);
  // test_integral = mytest.GaussInt_3D(0,M_PI,100, 0,M_PI,100, 0,M_PI,100);
  cout << "\n" << test_integral << endl;
}
