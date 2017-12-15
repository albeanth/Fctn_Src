#ifndef FEM_Integration
#define FEM_Integration_h 1

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
    double Error_Integrate1D( std::vector<double> FEMSoln, int nelsA, std::vector<int> orderA, std::vector<std::vector<int> > nodA, std::vector<double> xnodA, int maxordA);
    double FEM_Func_Integrate_1D( std::vector<double> FEMSoln, int nelsA, std::vector<int> orderA, std::vector<std::vector<int> > nodA, std::vector<double> xnodA, int maxordA);
    double FEM_Func_Integrate_2D_Serial( std::vector<double> FEMSoln,
                  int nelsA, std::vector<int> orderA, std::vector<std::vector<int> > nodA, std::vector<double> xnodA, int maxordA,
                  int nelsB, std::vector<int> orderB, std::vector<std::vector<int> > nodB, std::vector<double> xnodB, int maxordB);
    double FEM_Func_Integrate_2D_Parallel(int NTHREADS, std::vector<double> FEMSoln,
                  int nelsA, std::vector<int> orderA, std::vector<std::vector<int> > nodA, std::vector<double> xnodA, int maxordA,
                  int nelsB, std::vector<int> orderB, std::vector<std::vector<int> > nodB, std::vector<double> xnodB, int maxordB);
  private:
    double ExactFun(double x);
    double TestFun1D(double x);
    double TestFun2D(double x, double y);
    QuadParams getQPs(int maxord, QuadParams qps);
    ShapeFunction getShapeFuns(double x, int n, ShapeFunction shape);

};

#endif
