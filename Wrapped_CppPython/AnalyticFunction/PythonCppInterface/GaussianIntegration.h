#ifndef GaussianIntegration
#define GaussianIntegration_h 1

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>

// Struct for defining quadrature parameters
struct QuadParams{
  int nw;
  std::vector<double> xw;
  std::vector<double> w;
};

class GaussianIntegration{

  public:
    double GaussInt_1D(double xMin, double xMax, int xNode);
    double GaussInt_2D_Serial(double xMin, double xMax, int xNode,
                       double yMin, double yMax, int yNode);
    double GaussInt_2D_Parallel(int NTHREADS, double xMin, double xMax, int xNode,
                                double yMin, double yMax, int yNode);
    double GaussInt_3D_Serial(double xMin, double xMax, int xNode,
                       double yMin, double yMax, int yNode,
                       double zMin, double zMax, int zNode);
    double GaussInt_3D_Parallel(int NTHREADS, double xMin, double xMax, int xNode,
                       double yMin, double yMax, int yNode,
                       double zMin, double zMax, int zNode);

  private:
    double TestFun1D(double x);
    double TestFun2D(double x, double y);
    double TestFun3D(double x, double y, double z);
    std::vector<double> linspace(double a, double b, int NumElem);
    QuadParams getQPs(int maxord, QuadParams qps);
};

#endif
