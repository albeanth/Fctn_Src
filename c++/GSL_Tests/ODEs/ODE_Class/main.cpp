#include "Class_Test.h"
#include <iostream>
#include <vector>

/* NOTE: THE PROBLEM WITH GNU SCIENTIFIC LIBRARY IS IT IS WRITTEN IN C!! NOT C++....
         TO GET AROUND THIS, LOOK INTO BOOST. IT IS WRITTEN IN C++ */

int main(){
  int dim = 1;
  double RelTol = 1.0e-6;
  double AbsTol = 0.0;
  double mu = 0.5;
  // double T0 = 0.0;
  // double TEnd = 1.0;
  std::vector<double> Bounds;
  Bounds.push_back(0.0);
  Bounds.push_back(20.0);
  ODEClass_Test mytest; // initialize object for computing ODE
  Solution mySoln; // initialize struct for solution
  mySoln = mytest.compute_ODE(Bounds, dim, RelTol, AbsTol, mu);
  for(int i=0; i<mySoln.t.size(); i++){
    printf("%.5e %.5e\n",mySoln.t[i], mySoln.y[i]);
  }
  return 0;
}
