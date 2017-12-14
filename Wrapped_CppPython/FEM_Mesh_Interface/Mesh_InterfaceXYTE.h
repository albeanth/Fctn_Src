#ifndef Mesh_InterfaceXYTE
#define Mesh_InterfaceXYTE_h 1

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
    double Source_Integrate(
            char flag, std::vector<double> B_tmp, std::vector<double> C_tmp, std::vector<double> D_tmp, double a,
            int nelsB, std::vector<int> orderB, std::vector< std::vector< int > > nodB, std::vector< double > xnodB, int maxordB,
            int nelsC, std::vector<int> orderC, std::vector< std::vector< int > > nodC, std::vector< double > xnodC, int maxordC,
            int nelsD, std::vector<int> orderD, std::vector< std::vector< int > > nodD, std::vector< double > xnodD, int maxordD);

    std::vector<double> Error_Integrate4D(
            std::vector< std::vector< double > > A, std::vector< std::vector< double > > B, std::vector< std::vector< double > > C, std::vector< std::vector< double > > D,
            int nelsA, std::vector<int> orderA, std::vector< std::vector< int > > nodA, std::vector< double > xnodA, int maxordA,
            int nelsB, std::vector<int> orderB, std::vector< std::vector< int > > nodB, std::vector< double > xnodB, int maxordB,
            int nelsC, std::vector<int> orderC, std::vector< std::vector< int > > nodC, std::vector< double > xnodC, int maxordC,
            int nelsD, std::vector<int> orderD, std::vector< std::vector< int > > nodD, std::vector< double > xnodD, int maxordD);

    // functions to get functions for basis function calculations
    std::vector< std::vector< double > > UpdateH(
            double a, double b1, double c1, std::vector<double> B_tmp, std::vector<double> C_tmp,
            int nelsB, std::vector<int> orderB, std::vector< std::vector< int > > nodB, std::vector< double > xnodB, int maxordB,
            int nelsC, std::vector<int> orderC, std::vector< std::vector< int > > nodC, std::vector< double > xnodC, int maxordC);
    std::vector< std::vector< double > > UpdateH_Enr(
            double a, double b2, double c2, std::vector<double> B_tmp, std::vector<double> C_tmp, std::vector< std::vector< double > > B, std::vector< std::vector< double > > C,
            int nelsB, std::vector<int> orderB, std::vector< std::vector< int > > nodB, std::vector< double > xnodB, int maxordB,
            int nelsC, std::vector<int> orderC, std::vector< std::vector< int > > nodC, std::vector< double > xnodC, int maxordC);
    std::vector< std::vector< double > > UpdateG(
            std::vector<double> A_tmp, double x1, int nelsA, std::vector<int> orderA, std::vector< std::vector< int > > nodA, std::vector< double > xnodA, int maxordA,
            std::vector<double> B_tmp, double y1, int nelsB, std::vector<int> orderB, std::vector< std::vector< int > > nodB, std::vector< double > xnodB, int maxordB,
            std::vector<double> C_tmp, double c1, int nelsC, std::vector<int> orderC, std::vector< std::vector< int > > nodC, std::vector< double > xnodC, int maxordC);


  private:
    // functions used for to get reference solution
    double MMS_Source(double x, double y, double t, double E);
    double phi_fun(double x, double y, double t, double E);
    double phi_px(double x, double y, double t, double E);
    double phi_pxx(double x, double y, double t, double E);
    double phi_py(double x, double y, double t, double E);
    double phi_pyy(double x, double y, double t, double E);
    double phi_pt(double x, double y, double E);

    // functions describing cross sections and problem domain
    double Q(double E);
    double D(double x, double y, double E);
    double D_px(double y, double E);
    double D_py(double x, double E);
    double SigT(double x, double y, double E);
    double V(double E);
    double Xbnd = M_PI;
    double Ybnd = M_PI;

    // structs used for integrations
    QuadParams getQPs(int maxord, QuadParams qps);
    ShapeFunction getShapeFuns(double x, int n, ShapeFunction shape);

};

#endif
