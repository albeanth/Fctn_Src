#include "Class_Test.h"
#include <iostream>
#include <vector>

int main(){
  int dim = 1;
  double RelTol = 1.0e-6;
  double AbsTol = 0.0;
  Params mus;
  mus.mu0 = 0.5;
  mus.mu1 = 0.5;
  std::vector<double> Bounds = {0.0,1.0};
  ODEClass_Test mytest; // initialize object for computing ODE
  Solution mySoln; // initialize struct for solution
  mySoln = mytest.compute_ODE(Bounds, dim, RelTol, AbsTol);//, mus);
  return 0;
}
