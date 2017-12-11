#ifndef Mesh_InterfaceXYT
#define Mesh_InterfaceXYT_h 1

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <string.h>
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
    double Source_Integrate(char flag, std::vector<double> B_tmp, std::vector<double> C_tmp, double a,
                        int nelsB, std::vector<int> orderB, std::vector< std::vector< int > > nodB, std::vector< double > xnodB, int maxordB,
                        int nelsC, std::vector<int> orderC, std::vector< std::vector< int > > nodC, std::vector< double > xnodC, int maxordC);
    std::vector<double> Error_Integrate3D(
      std::vector< std::vector< double > > X, std::vector< std::vector< double > > T, std::vector< std::vector< double > > E,
      int nelsA, std::vector<int> orderA, std::vector< std::vector< int > > nodA, std::vector< double > xnodA, int maxordA,
      int nelsB, std::vector<int> orderB, std::vector< std::vector< int > > nodB, std::vector< double > xnodB, int maxordB,
      int nelsC, std::vector<int> orderC, std::vector< std::vector< int > > nodC, std::vector< double > xnodC, int maxordC);
  private:
    double MMS_Source(double a, double b, double c);
    double phi_fun(double x, double y, double t);
    double phi_px(double x, double y, double t);
    double phi_pxx(double x, double y, double t);
    double phi_py(double x, double y, double t);
    double phi_pyy(double x, double y, double t);
    double phi_pt(double x, double y, double t);
    double D(double x, double y);
    double D_px(double y);
    double D_py(double x);
    double SigAbs(double x, double y);
    double v=1.5;
    double Xbnd = M_PI;
    double Ybnd = M_PI;
    // double Xbnd = sqrt(2.0)*sqrt(M_PI);
    // double Ybnd = sqrt(2.0)*sqrt(M_PI);
    QuadParams getQPs(int maxord, QuadParams qps);
    ShapeFunction getShapeFuns(double x, int n, ShapeFunction shape);

};

#endif
