#ifndef Mesh_InterfaceXT
#define Mesh_InterfaceXT_h 1

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
/*
This interface allows for python generated 1D FEM meshes to be read into c++ via swig.
*/

// Struct for defining quadrature parameters
struct QuadParams{
  int nw;
  std::vector<double> xw;
  std::vector<double> w;
};

struct ShapeFunction{
  std::vector<double> psi;
  std::vector<double> dpsi;
};

class GaussianIntegration{

  public:
    double Source_Integrate(int NUMTHREADS,
            std::vector<double> X, double Tval, std::vector<double> X_tmp,
            int nels, std::vector<int> order, std::vector< std::vector< int > > nod, std::vector< double > xnod, int maxord);
  private:
    QuadParams getQPs(int maxord, QuadParams qps);
    ShapeFunction getShapeFuns(double x, int n, ShapeFunction shape);

};

#endif
