#include "GaussianIntegration.h"
#include <stdlib.h>

using namespace std;

double GaussianIntegration::TestFun(double x, double y){
  return cos(x*y);
}

void GaussianIntegration::SetupGrid(double xMin, double xMax, int xNum, double yMin, double yMax, int yNum){
  xVec = linspace(xMin, xMax,xNum);
  yVec = linspace(yMin, yMax, yNum);
}

double GaussianIntegration::GaussInt_2D ( ){
  std::vector<double> xRange;
  std::vector<double> yRange;
  xRange = xVec;
  yRange = yVec;
  QuadParams qps;
  qps = getQPs(2, qps);

  double u_h;
  u_h = 0.0;
  double xL,xR,dx,x;
  double yL,yR,dy,y;
  for (unsigned int idx=0; idx<xRange.size()-1; idx++){
    xL = xRange[idx];
    xR = xRange[idx+1];
    dx = (xR-xL)/2.;

    for (int l1=0; l1<qps.nw; l1++){
      x = xL + (1 + qps.xw[l1])*dx;

      for (unsigned int idy=0; idy<yRange.size()-1; idy++){
        yL = yRange[idy];
        yR = yRange[idy+1];
        dy = (yR-yL)/2.;

        for (int l2=0; l2<qps.nw; l2++){
          y = yL + (1 + qps.xw[l2])*dy;
          //complete integrations
          u_h += TestFun(x,y) *qps.w[l1]*dx*qps.w[l2]*dy;
        }
      }
    }
  }
  return u_h;
}

vector<double>  GaussianIntegration::linspace(double a, double b, int NumElem){
  std::vector<double> space;
  double dx;
  dx = (b-a)/NumElem;
  for (double i=a; i<=b; i+=dx){
    space.push_back(i);
  }
 return space;
}

QuadParams GaussianIntegration::getQPs(int maxord, QuadParams qps){
  int nw;
  std::vector<double> xw;
  std::vector<double> w;

  if (maxord == 1) {
    nw = 1;
    xw.push_back(0.0);
    w.push_back(2.0);
  }
  else if (maxord == 2){
    nw = 2;
    xw.push_back(-1./sqrt(3.));
    xw.push_back(-xw[0]);
    w.push_back(1.);
    w.push_back(1.);
  }
  else if (maxord == 3){
    nw = 3;
    xw.push_back(sqrt(3./5.));
    xw.push_back(0.);
    xw.push_back(-xw[0]);
    w.push_back(5./9.);
    w.push_back(8./9.);
    w.push_back(w[0]);
  }
  else if (maxord == 4){
    nw = 4;
    xw.push_back(sqrt(3./7. + 2./7.*sqrt(6./5.)));
    xw.push_back(sqrt(3./7. - 2./7.*sqrt(6./5.)));
    xw.push_back(-xw[1]);
    xw.push_back(-xw[0]);
    w.push_back((18.-sqrt(30.))/36.);
    w.push_back((18.+sqrt(30.))/36.);
    w.push_back(w[1]);
    w.push_back(w[0]);
  }
  else if (maxord == 5){
    nw = 5;
    xw.push_back(-1./3.*sqrt(5.+2.*sqrt(10./7.)));
    xw.push_back(-1./3.*sqrt(5.-2.*sqrt(10./7.)));
    xw.push_back(0.0);
    xw.push_back(-xw[1]);
    xw.push_back(-xw[0]);
    w.push_back((322.-13.*sqrt(70.))/900.);
    w.push_back((322.+13.*sqrt(70.))/900.);
    w.push_back(128./225.);
    w.push_back(w[1]);
    w.push_back(w[0]);
  }
  else{
    cout << "Incorrect quadrature order." << endl;
    cout << "Order = " << maxord << " does not exist." << endl;
    exit(1);
  }
  qps.nw = nw;
  qps.xw = xw;
  qps.w = w;
  return qps;
}
