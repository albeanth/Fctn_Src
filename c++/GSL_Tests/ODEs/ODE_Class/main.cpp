#include "Class_Test.hpp"
#include <iostream>
#include <vector>

int main(){
  // --------------------------------------------
  int dim = 1;
  double RelTol = 0.0;
  double AbsTol = 1.0E-6;
  std::vector<double> Bounds = {0.0,1.0};
  // --------------------------------------------
  ODEClass_Test mytest(AbsTol, RelTol, Bounds); // initialize object for computing ODE
  Solution mySoln;                              // initialize struct for solution
  // --------------------------------------------
  mySoln = mytest.compute_ODE(dim);
  mySoln = mytest.compute_Sparse_ODE(dim);

  return 0;
}
