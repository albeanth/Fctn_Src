#include "FEM_Integration.h"
#include <omp.h>
using namespace std;

double GaussianIntegration::Error_Integrate1D(vector<double> FEMSoln, int nelsA, vector<int> orderA, vector<vector<int> > nodA, vector<double> xnodA, int maxordA){
  /*
  This function is to test an error comparison for a FE solution and a known analytic solution.
  */
  QuadParams qpsA;
  qpsA = getQPs(maxordA, qpsA);
  double aL; double aR; double da; double a; int mynumX; ShapeFunction Xshape;

  double tmp_uh;
  double uh_int; double diff;

  uh_int = 0.0; diff = 0.0;
  for(int Ael=0; Ael<nelsA; Ael++){ // for each element in dimension A, get l/r bounds, and transformation
    aL = xnodA[nodA[Ael][0]];
    aR = xnodA[nodA[Ael][orderA[Ael]-1]];
    da = (aR-aL)/2.0;
    for (int l1=0; l1<qpsA.nw; l1++){ // over the number of weights in the quadrature order
      a = aL + (1.0+qpsA.xw[l1])*da; // get the coordinate in the real mesh
      Xshape = getShapeFuns(qpsA.xw[l1], orderA[Ael], Xshape);
      tmp_uh = 0.0;
      for (int k1 = 0; k1<orderA[Ael]; k1++){ // over the nodes in the element based on mesh order
        mynumX = nodA[Ael][k1];
        tmp_uh = tmp_uh + FEMSoln[mynumX] *Xshape.psi[k1];
      }
      // complete integration over element
      uh_int = uh_int + tmp_uh *qpsA.w[l1]*da;
      diff = diff + pow( ExactFun(a)-tmp_uh, 2) *qpsA.w[l1]*da;
    }
  }
  printf("    L2 diff = %.4f, u_h Int = %.4f\n",diff,uh_int);
  return diff;
}

double GaussianIntegration::FEM_Func_Integrate_1D(vector<double> FEMSoln, int nelsA, vector<int> orderA, vector<vector<int> > nodA, vector<double> xnodA, int maxordA){
  /*
  this fucntion is for testing the integration of an integral over one FEM fucntion and a known analytic function.
  \int FEM*Function

  This is useful for the known analytic fucntions that are the coefficients in the BVP and IVP solves.
      - Think X_old*D(x) etc.
  */
  QuadParams qpsA;
  qpsA = getQPs(maxordA, qpsA);
  double aL; double aR; double da; double a; int mynumX; ShapeFunction Xshape;

  double tmp_uh;
  double uh_int;

  uh_int = 0.0;
  for(int Ael=0; Ael<nelsA; Ael++){ // for each element in dimension A, get l/r bounds, and transformation
    aL = xnodA[nodA[Ael][0]];
    aR = xnodA[nodA[Ael][orderA[Ael]-1]];
    da = (aR-aL)/2.0;
    for (int l1=0; l1<qpsA.nw; l1++){ // over the number of weights in the quadrature order
      a = aL + (1.0+qpsA.xw[l1])*da; // get the coordinate in the real mesh
      Xshape = getShapeFuns(qpsA.xw[l1], orderA[Ael], Xshape);
      tmp_uh = 0.0;
      for (int k1 = 0; k1<orderA[Ael]; k1++){ // over the nodes in the element based on mesh order
        mynumX = nodA[Ael][k1];
        tmp_uh = tmp_uh + FEMSoln[mynumX] *Xshape.psi[k1];
      }
      // complete integration over element
      uh_int = uh_int + pow(tmp_uh,3.0)*TestFun1D(a) *qpsA.w[l1]*da;
    }
  }
  return uh_int;
}

double GaussianIntegration::FEM_Func_Integrate_2D_Serial(vector<double> FEMSoln,
  int nelsA, vector<int> orderA, vector<vector<int> > nodA, vector<double> xnodA, int maxordA,
  int nelsB, vector<int> orderB, vector<vector<int> > nodB, vector<double> xnodB, int maxordB){

  /*
  this fucntion is for testing the PGD source. it completes an integral over two FEM fucntions and a known analytic function.
  \int\int FEM*FEM*Function

  can also be used for finding r coefficients that are integrated over x and y space.
  */

  QuadParams qpsA;
  qpsA = getQPs(maxordA, qpsA);
  double aL; double aR; double da; double a; int mynumA; ShapeFunction Ashape;
  QuadParams qpsB;
  qpsB = getQPs(maxordB, qpsB);
  double bL; double bR; double db; double b; int mynumB; ShapeFunction Bshape;

  double tmp_uhA; double tmp_uhB; double uh_int;
  uh_int = 0.0;
  for(int Ael=0; Ael<nelsA; Ael++){ // for each element in dimension A, get l/r bounds, and transformation
    aL = xnodA[nodA[Ael][0]];
    aR = xnodA[nodA[Ael][orderA[Ael]-1]];
    da = (aR-aL)/2.0;
    for (int l1=0; l1<qpsA.nw; l1++){ // over the number of weights in the quadrature order
      a = aL + (1.0+qpsA.xw[l1])*da; // get the coordinate in the real mesh
      Ashape = getShapeFuns(qpsA.xw[l1], orderA[Ael], Ashape);
      tmp_uhA = 0.0;
      for (int k1 = 0; k1<orderA[Ael]; k1++){ // over the nodes in the element based on mesh order
        mynumA = nodA[Ael][k1];
        tmp_uhA = tmp_uhA + FEMSoln[mynumA] *Ashape.psi[k1];
      }

      for (int Bel=0; Bel<nelsB; Bel++){
        bL = xnodB[nodB[Bel][0]];
        bR = xnodB[nodB[Bel][orderB[Bel]-1]];
        db = (bR-bL)/2.0;
        for (int l2=0; l2<qpsB.nw; l2++){
          b = bL + (1.0+qpsB.xw[l2])*db;
          Bshape = getShapeFuns(qpsB.xw[l2], orderB[Bel], Bshape);
          tmp_uhB = 0.0;
          for (int k2 = 0; k2<orderB[Bel]; k2++){
            mynumB = nodB[Bel][k2];
            tmp_uhB = tmp_uhB + FEMSoln[mynumB]*Bshape.psi[k2];
          }
          // complete integration over element
          uh_int = uh_int + tmp_uhA*tmp_uhB*TestFun2D(a,b) *qpsA.w[l1]*da*qpsB.w[l2]*db;
        }
      }
    }
  }
  return uh_int;
}

double GaussianIntegration::FEM_Func_Integrate_2D_Parallel(int NTHREADS, vector<double> FEMSoln,
  int nelsA, vector<int> orderA, vector<vector<int> > nodA, vector<double> xnodA, int maxordA,
  int nelsB, vector<int> orderB, vector<vector<int> > nodB, vector<double> xnodB, int maxordB){

  /*
  this fucntion is for testing the PGD source. it completes an integral over two FEM fucntions and a known analytic function.
  \int\int FEM*FEM*Function

  can also be used for finding r coefficients that are integrated over x and y space.
  */

  omp_set_num_threads(NTHREADS);

  QuadParams qpsA;
  qpsA = getQPs(maxordA, qpsA);
  double aL; double aR; double da; double a; int mynumA; ShapeFunction Ashape;
  QuadParams qpsB;
  qpsB = getQPs(maxordB, qpsB);
  double bL; double bR; double db; double b; int mynumB; ShapeFunction Bshape;

  double tmp_uhA; double tmp_uhB; double uh_int;
  uh_int = 0.0;

#pragma omp parallel for default(none),shared(qpsA,qpsB, FEMSoln, nelsA,orderA,nodA,xnodA,maxordA, nelsB,orderB,nodB,xnodB,maxordB),private(aL,aR,da,a,mynumA,Ashape,tmp_uhA, bL,bR,db,b,mynumB,Bshape,tmp_uhB),reduction(+:uh_int)
  for(int Ael=0; Ael<nelsA; Ael++){ // for each element in dimension A, get l/r bounds, and transformation
    aL = xnodA[nodA[Ael][0]];
    aR = xnodA[nodA[Ael][orderA[Ael]-1]];
    da = (aR-aL)/2.0;
    for (int l1=0; l1<qpsA.nw; l1++){ // over the number of weights in the quadrature order
      a = aL + (1.0+qpsA.xw[l1])*da; // get the coordinate in the real mesh
      Ashape = getShapeFuns(qpsA.xw[l1], orderA[Ael], Ashape);
      tmp_uhA = 0.0;
      for (int k1 = 0; k1<orderA[Ael]; k1++){ // over the nodes in the element based on mesh order
        mynumA = nodA[Ael][k1];
        tmp_uhA = tmp_uhA + FEMSoln[mynumA] *Ashape.psi[k1];
      }

      for (int Bel=0; Bel<nelsB; Bel++){
        bL = xnodB[nodB[Bel][0]];
        bR = xnodB[nodB[Bel][orderB[Bel]-1]];
        db = (bR-bL)/2.0;
        for (int l2=0; l2<qpsB.nw; l2++){
          b = bL + (1.0+qpsB.xw[l2])*db;
          Bshape = getShapeFuns(qpsB.xw[l2], orderB[Bel], Bshape);
          tmp_uhB = 0.0;
          for (int k2 = 0; k2<orderB[Bel]; k2++){
            mynumB = nodB[Bel][k2];
            tmp_uhB = tmp_uhB + FEMSoln[mynumB]*Bshape.psi[k2];
          }
          // complete integration over element
          uh_int = uh_int + tmp_uhA*tmp_uhB*TestFun2D(a,b) *qpsA.w[l1]*da*qpsB.w[l2]*db;
        }
      }
    }
  }
  return uh_int;
}


double GaussianIntegration::ExactFun(double x){
  return sin(x);
  // return pow(x,2.0);
  // return cos(pow(x,2.0));
}

double GaussianIntegration::TestFun1D(double x){
  return (x+pow(x,2));
}

double GaussianIntegration::TestFun2D(double x,double y){
  return x+y;
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
    y.push_back( 1./16.*(-9.*pow(x,3.) + 9.*pow(x,2.) + x - 1.) );
    y.push_back( 1./16.*(27.*pow(x,3.) - 9.*pow(x,2.) - 27.*x + 9.) );
    y.push_back( 1./16.*(-27.*pow(x,3.) - 9.*pow(x,2.) + 27.*x + 9.) );
    y.push_back( 1./16.*(9.*pow(x,3.) + 9.*pow(x,2.) - x - 1.) );
    dy.push_back( 1./16.*(-27.*pow(x,2.) + 18.*x + 1. ) );
    dy.push_back( 1./16.*(81.*pow(x,2.) - 18.*x - 27.) );
    dy.push_back( 1./16.*(-81.*pow(x,2.) - 18.*x + 27.) );
    dy.push_back( 1./16.*(27.*pow(x,2.) + 18.*x - 1. ) );
  }
  else{
    cout << "Order = " << n << " shape function does not exist." << endl;
    exit(1);
  }
  shape.psi = y;
  shape.dpsi = dy;
  return shape;
}
