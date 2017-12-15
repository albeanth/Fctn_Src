#include "GaussianIntegration.h"
#include <stdlib.h>
#include <omp.h>

using namespace std;

double GaussianIntegration::TestFun1D(double x){
  return cos(x);
}
double GaussianIntegration::GaussInt_1D (double xMin, double xMax, int xNode){
  std::vector<double> xRange;
  xRange = linspace(xMin, xMax, xNode);
  QuadParams qps;
  qps = getQPs(2, qps);

  double u_h;
  u_h = 0.0;
  double xL,xR,dx,x;
  for (unsigned int idx=0; idx<xRange.size()-1; idx++){
    xL = xRange[idx];
    xR = xRange[idx+1];
    dx = (xR-xL)/2.;

    for (int l1=0; l1<qps.nw; l1++){
      x = xL + (1 + qps.xw[l1])*dx;
      //complete integrations
      u_h += TestFun1D(x) *qps.w[l1]*dx;
      }
    }
  return u_h;
}

double GaussianIntegration::TestFun2D(double x, double y){
  return pow(y,(x+0.5))*sin(x);
}
double GaussianIntegration::GaussInt_2D_Serial (double xMin, double xMax, int xNode, double yMin, double yMax, int yNode){
  std::vector<double> xRange;
  std::vector<double> yRange;
  xRange = linspace(xMin, xMax, xNode);
  yRange = linspace(yMin, yMax, yNode);
  QuadParams qps;
  qps = getQPs(3, qps);

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
          u_h += TestFun2D(x,y) *qps.w[l1]*dx*qps.w[l2]*dy;
        }
      }
    }
  }
  return u_h;
}

double GaussianIntegration::GaussInt_2D_Parallel (int NTHREADS, double xMin, double xMax, int xNode, double yMin, double yMax, int yNode){

  std::vector<double> xRange;
  std::vector<double> yRange;
  xRange = linspace(xMin, xMax, xNode);
  yRange = linspace(yMin, yMax, yNode);
  QuadParams qps;
  qps = getQPs(3, qps);

  double u_h;
  u_h = 0.0;
  double xL,xR,dx,x;
  double yL,yR,dy,y;

  #pragma omp parallel for default(none),shared(xRange,yRange,qps),private(xL,xR,dx,x,yL,yR,dy,y),reduction(+:u_h)
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
          u_h += TestFun2D(x,y) *qps.w[l1]*dx*qps.w[l2]*dy;
        }
      }
    }
  }
  return u_h;
}

double GaussianIntegration::TestFun3D(double x, double y, double z){
  return cos(x*y*z);
}
double GaussianIntegration::GaussInt_3D_Serial (double xMin, double xMax, int xNode, double yMin, double yMax, int yNode, double zMin, double zMax, int zNode){
  std::vector<double> xRange;
  std::vector<double> yRange;
  std::vector<double> zRange;
  xRange = linspace(xMin, xMax, xNode);
  yRange = linspace(yMin, yMax, yNode);
  zRange = linspace(zMin, zMax, zNode);
  QuadParams qps;
  qps = getQPs(4, qps);

  double u_h;
  u_h = 0.0;
  double xL,xR,dx,x;
  double yL,yR,dy,y;
  double zL,zR,dz,z;
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

          for (unsigned int idz=0; idz<zRange.size()-1; idz++){
            zL = zRange[idz];
            zR = zRange[idz+1];
            dz = (zR-zL)/2;
            for (int l3=0; l3<qps.nw; l3++){
              z = zL + (1 + qps.xw[l3])*dz;

              //complete integrations
              u_h += TestFun3D(x,y,z) *qps.w[l1]*dx *qps.w[l2]*dy *qps.w[l3]*dz;
            }
          }
        }
      }
    }
  }
  return u_h;
}

double GaussianIntegration::GaussInt_3D_Parallel (int NTHREADS, double xMin, double xMax, int xNode, double yMin, double yMax, int yNode, double zMin, double zMax, int zNode){
  
  omp_set_num_threads(NTHREADS);

  std::vector<double> xRange;
  std::vector<double> yRange;
  std::vector<double> zRange;
  xRange = linspace(xMin, xMax, xNode);
  yRange = linspace(yMin, yMax, yNode);
  zRange = linspace(zMin, zMax, zNode);
  QuadParams qps;
  qps = getQPs(4, qps);

  double u_h;
  u_h = 0.0;
  double xL,xR,dx,x;
  double yL,yR,dy,y;
  double zL,zR,dz,z;

#pragma omp parallel for default(none),shared(xRange,yRange,zRange,qps),private(xL,xR,dx,x,yL,yR,dy,y,zL,zR,dz,z),reduction(+:u_h)
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

          for (unsigned int idz=0; idz<zRange.size()-1; idz++){
            zL = zRange[idz];
            zR = zRange[idz+1];
            dz = (zR-zL)/2;
            for (int l3=0; l3<qps.nw; l3++){
              z = zL + (1 + qps.xw[l3])*dz;

              //complete integrations
              u_h += TestFun3D(x,y,z) *qps.w[l1]*dx *qps.w[l2]*dy *qps.w[l3]*dz;
            }
          }
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
