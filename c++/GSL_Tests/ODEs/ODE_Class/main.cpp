#include "Class_Test.hpp"
#include <iostream>
#include <vector>

int main(){
  // --------------------------------------------
  int dim = 1;
  double RelTol = 1E-6;
  double AbsTol = 1.0E-9;
  std::vector<double> Bounds = {0.0,1.0};
  // --------------------------------------------
  ODEClass_Test mytest(AbsTol, RelTol, Bounds, dim); // initialize object for computing ODE
  // --------------------------------------------
  mytest.compute_ODE();
  // mytest.compute_Sparse_ODE(dim);

  std::vector<double> error;
  mytest.compute_error(error);
  printf("     t (sec)\t\t    GSL Soln\t\t  Relative Error\n");
  for (int i = 0; i < mytest.soln.t.size(); i++){
    printf("%.12e\t%.12e\t%.12e\n", mytest.soln.t[i], mytest.soln.y[i], error[i]);
  }

  return 0;
}
