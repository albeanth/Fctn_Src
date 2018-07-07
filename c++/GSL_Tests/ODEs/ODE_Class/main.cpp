#include "Class_Test.h"
#include <iostream>
#include <vector>

int main(){
  int dim = 1;
  double RelTol = 1.0e-6;
  double AbsTol = 0.0;
  // double mu = 0.5;
  Params mus = {0.35,0.15};
  std::vector<double> Bounds;
  Bounds.push_back(0.0);
  Bounds.push_back(20.0);
  ODEClass_Test mytest; // initialize object for computing ODE
  Solution mySoln; // initialize struct for solution
  mySoln = mytest.compute_ODE(Bounds, dim, RelTol, AbsTol, mus);
  for(int i=0; i<mySoln.t.size(); i++){
    printf("%.5e %.5e\n",mySoln.t[i], mySoln.y[i]);
  }
  return 0;
}
