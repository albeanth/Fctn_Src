#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include "GaussianIntegration.h"
using namespace std;

int main()
{
  GaussianIntegration mytest;
  mytest.SetupGrid(0,M_PI,100, 0,M_PI,100);

  double test_integral;
  test_integral = mytest.GaussInt_2D();
  cout << "\n" << test_integral << endl;
}
