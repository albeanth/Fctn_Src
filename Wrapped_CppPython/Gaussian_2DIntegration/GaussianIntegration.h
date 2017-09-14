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
    double GaussInt_2D();
    void SetupGrid(double xMin, double xMax, int xNum, double yMin, double yMax, int yNum);

  private:
    std::vector<double> linspace(double a, double b, int NumElem);
    QuadParams getQPs(int maxord, QuadParams qps);
    double TestFun(double x, double y);
    std::vector<double> xVec;
    std::vector<double> yVec;
};

#endif
