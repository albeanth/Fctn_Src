#include "Mesh_InterfaceXT.h"
#include <stdlib.h>
#include <omp.h>

using namespace std;

double GaussianIntegration::get_AbsX_NL(int NUMTHREADS,
                    vector<double> X, double Tval, <double> X_tmp,
                    int nels, vector<int> order, vector< vector< int > > nod, vector< double > xnod, int maxord){

  omp_set_num_threads(NUMTHREADS);

  QuadParams qps;
  qps = getQPs(maxord, qps);
  double xL; double xR; double dx; double x; int mynum; ShapeFunction shape;

  double x_val; double xi_val;
  double EnrSum; double C1;
  int TotNumEnr;
  TotNumEnr = X.size();

  vector< vector< double > > xAbs;
  double xAbs1=0.0; double xAbs2=0.0;
  vector<double> xAbs3(TotNumEnr);
  vector<double> xAbs4(TotNumEnr);

  #pragma omp parallel for default(none),\
  shared(qps, nels,order,nod,xnod, X,Tval,X_tmp),	\
  private(xL,xR,dx,x,mynum,shape, x_val,xi_val,EnrSum,C1), \
  reduction(+:xAbs1,xAbs2,xAbs3,xAbs4)
  for (int el=0; el<nels; el++){
    xL = xnod[nod[el][0]];
    xR = xnod[nod[el][order[el]-1]];
    dx = (xR-xL)/2.0;
    for (int l2=0; l2<qps.nw; l2++){
      x = xL + (1.0+qps.xw[l2])*dx;
      shape = getShapeFuns(qps.xw[l2], order[el], shape);
      x_val = 0.0;
      for (int k2 = 0; k2<order[el]; k2++){
        mynum = nod[el][k2];
        x_val += X_tmp[mynum]*shape.psi[k2];
      }
      // calculate C1
      EnrSum = 0.0;
      for (int ide = 0; ide<TotNumEnr; ide++){
        xi_val = 0.0;
        for (int k2 = 0; k2<order[el]; k2++){
          mynum = nod[el][k2];
          xi_val += X[ide][mynum]*shape.psi[k2];
        }
        EnrSum += xi_val*Tval
      }
      C1 = EnrSum
      // complete integration over element
      xAbs1 += (pow(x_val,2.0)*pow(C1,2.0)) * qps.w[l2]*dx;
      xAbs2 += (pow(x_val,2.0)*C1)          * qps.w[l2]*dx;
      for (int ide = 0; ide<TotNumEnr; ide++){
        xi_val = 0.0;
        for (int k2 = 0; k2<order[el]; k2++){
          mynum = nod[el][k2];
          xi_val += X[ide][mynum]*shape.psi[k2];
        }
        xAbs3[enr] += (x_val*xi_val*C1)          * qps.w[l2]*dx;
        xAbs4[enr] += (x_val*xi_val*pow(C1,2.0)) * qps.w[l2]*dx;
      }
    }
  }
  xAbs.push_back(xAbs1); xAbs.push_back(xAbs2);
  xAbs.push_back(xAbs3); xAbs.push_back(xAbs4);
  return xAbs;
}

QuadParams GaussianIntegration::getQPs(int maxord, QuadParams qps){
  int nw;
  vector<double> xw;
  vector<double> w;

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

ShapeFunction GaussianIntegration::getShapeFuns(double x, int n, ShapeFunction shape){
  /*
  shape function on reference element (-1,1)
  n = 2: linear
  n = 3: quadratic
  */
  std::vector<double> y;
  std::vector<double> dy;
  if (n==2){
    y.push_back(0.5*(1.0-x));
    y.push_back(0.5*(1.0+x));
    dy.push_back(-0.5);
    dy.push_back(0.5);
  }
  else if (n==3){
    y.push_back((pow(x,2)-x)/2.0);
    y.push_back(1-pow(x,2));
    y.push_back((pow(x,2)+x)/2.0);
    dy.push_back(x-1./2.);
    dy.push_back(-2.0*x);
    dy.push_back(x+1./2.);
  }
  else if (n==4){
    y.push_back( 1./16.*(-9.*pow(x,3) + 9.*pow(x,2) + x - 1.) );
    y.push_back( 1./16.*(27.*pow(x,3) - 9.*pow(x,2) - 27.*x + 9.) );
    y.push_back( 1./16.*(-27.*pow(x,3) - 9.*pow(x,2) + 27.*x + 9.) );
    y.push_back( 1./16.*(9.*pow(x,3) + 9.*pow(x,2) - x - 1.) );
    dy.push_back( 1./16.*(-27.*pow(x,2) + 18.*x + 1.) );
    dy.push_back( 1./16.*(81.*pow(x,2) - 18.*x - 27.) );
    dy.push_back( 1./16.*(-81.*pow(x,2) - 18.*x + 27.) );
    dy.push_back( 1./16.*(27.*pow(x,2) + 18.*x - 1.) );
  }
  else{
    cout << "Order = " << n << " shape function does not exist." << endl;
    exit(1);
  }
  shape.psi = y;
  shape.dpsi = dy;
  return shape;
}
