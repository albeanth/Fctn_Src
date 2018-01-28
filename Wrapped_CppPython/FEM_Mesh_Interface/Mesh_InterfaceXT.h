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
    double Source_Integrate(
      int selector,int NUMTHREADS, char flag, std::vector<double> B_tmp, double a,
      int nelsB, std::vector<int> orderB, std::vector< std::vector< int > > nodB, std::vector< double > xnodB, int maxordB);
    std::vector<double> Error_Integrate2D(
      int selector,int NUMTHREADS, std::vector< std::vector< double > > A, std::vector< std::vector< double > > B,
      int nelsA, std::vector<int> orderA, std::vector< std::vector< int > > nodA, std::vector< double > xnodA, int maxordA,
      int nelsB, std::vector<int> orderB, std::vector< std::vector< int > > nodB, std::vector< double > xnodB, int maxordB);
  private:
    double MMS_Source(int selector, double x, double t);
    double phi_fun(int selector, double x, double t);
    double phi_px(int selector, double x, double t);
    double phi_pxx(int selector, double x, double t);
    double phi_pt(int selector, double x, double t);
    double D(double x);
    double SigAbs(int selector, double x);
    double Xend(int selector);
    const double v=1.5;
    QuadParams getQPs(int maxord, QuadParams qps);
    ShapeFunction getShapeFuns(double x, int n, ShapeFunction shape);

};

#endif
