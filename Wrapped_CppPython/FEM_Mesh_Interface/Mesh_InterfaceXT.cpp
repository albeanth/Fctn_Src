#include "Mesh_InterfaceXT.h"
#include <stdlib.h>
#include <omp.h>

using namespace std;

vector<double> GaussianIntegration::Error_Integrate2D(
      int NUMTHREADS, vector< vector< double > > A, vector< vector< double > > B,
      int nelsA, vector<int> orderA, vector< vector< int > > nodA, vector< double > xnodA, int maxordA,
      int nelsB, vector<int> orderB, vector< vector< int > > nodB, vector< double > xnodB, int maxordB){

  omp_set_num_threads(NUMTHREADS);

  int TotNumEnr;
  TotNumEnr = A.size();

  QuadParams qpsA;
  qpsA = getQPs(maxordA, qpsA);
  QuadParams qpsB;
  qpsB = getQPs(maxordB, qpsB);
  double aL; double aR; double da; double a; int mynumA; ShapeFunction Ashape;
  double bL; double bR; double db; double b; int mynumB; ShapeFunction Bshape;
  vector<double> l2Error(TotNumEnr);
  double pgdA; double pgdB;
  double diff; double pgd_int; double pgdEnr_int;
  double tmp;

  for (int enr=0; enr<TotNumEnr; enr++){
    diff = 0.0; pgd_int = 0.0; pgdEnr_int = 0.0;
    vector<double> tmpEnr_pgdA(enr+1); vector<double> tmpEnr_pgdB(enr+1);

   #pragma omp parallel for default(none),\
   shared(enr, qpsA,qpsB, A,B, nelsA,orderA,nodA,xnodA, nelsB,orderB,nodB,xnodB), \
   private(aL,aR,da,a,Ashape,pgdA,mynumA, bL,bR,db,b,Bshape,pgdB,mynumB, tmp), \
   firstprivate(tmpEnr_pgdA, tmpEnr_pgdB),\
   reduction(+:diff,pgd_int,pgdEnr_int)
    for (int Ael=0; Ael<nelsA; Ael++){ // for each element in dimension A, get l/r bounds, and transformation
      aL = xnodA[nodA[Ael][0]];
      aR = xnodA[nodA[Ael][orderA[Ael]-1]];
      da = (aR-aL)/2.0;
      for (int l1=0; l1<qpsA.nw; l1++){ // over the number of weights in the quadrature order
        a = aL + (1.0+qpsA.xw[l1])*da; // get the coordinate in the real mesh
        Ashape = getShapeFuns(qpsA.xw[l1], orderA[Ael], Ashape);
        pgdA = 0.0; fill(tmpEnr_pgdA.begin(),tmpEnr_pgdA.end(),0.0);
        for (int k1 = 0; k1<orderA[Ael]; k1++){ // over the nodes in the element based on mesh order
          mynumA = nodA[Ael][k1];
          pgdA = pgdA + A[enr][mynumA]*Ashape.psi[k1];
          for (int ide = enr; ide > -1; ide--){
            tmpEnr_pgdA[ide] = tmpEnr_pgdA[ide] + A[ide][mynumA]*Ashape.psi[k1];
          }
        }

        for (int Bel=0; Bel<nelsB; Bel++){
          bL = xnodB[nodB[Bel][0]];
          bR = xnodB[nodB[Bel][orderB[Bel]-1]];
          db = (bR-bL)/2.0;
          for (int l2=0; l2<qpsB.nw; l2++){
            b = bL + (1.0+qpsB.xw[l2])*db;
            Bshape = getShapeFuns(qpsB.xw[l2], orderB[Bel], Bshape);
            pgdB = 0.0; fill(tmpEnr_pgdB.begin(),tmpEnr_pgdB.end(),0.0);
            for (int k2 = 0; k2<orderB[Bel]; k2++){
              mynumB = nodB[Bel][k2];
              pgdB = pgdB + B[enr][mynumB]*Bshape.psi[k2];
              for (int ide = enr; ide > -1; ide--){
                tmpEnr_pgdB[ide] = tmpEnr_pgdB[ide] + B[ide][mynumB]*Bshape.psi[k2];
              }
            }

            // complete integration over element
            pgdEnr_int = pgdEnr_int + pgdA*pgdB  *qpsA.w[l1]*da *qpsB.w[l2]*db; //integral of step enr
            tmp = 0.0;
            for (int ide = enr; ide > -1; ide--){
              pgd_int = pgd_int + tmpEnr_pgdA[ide]*tmpEnr_pgdB[ide] *qpsA.w[l1]*da *qpsB.w[l2]*db; //total integral PGD up to step enr
              tmp = tmp + tmpEnr_pgdA[ide]*tmpEnr_pgdB[ide];
            }
            diff = diff + pow( phi_fun(a,b)-tmp, 2.0) *qpsA.w[l1]*da *qpsB.w[l2]*db;
          }
        }
      }
    }
    l2Error[enr] = sqrt(diff) ;
    printf("    %d, L2Error = %.6e, PGD Int = %.7f, Tot PGD Int = %.7f\n",enr+1,l2Error[enr],pgdEnr_int,pgd_int);
  }
  return l2Error;
}

double GaussianIntegration::Source_Integrate(int NUMTHREADS, char flag, vector<double> B_tmp, double a,
                    int nelsB, vector<int> orderB, vector< vector< int > > nodB, vector< double > xnodB, int maxordB){

  omp_set_num_threads(NUMTHREADS);

  QuadParams qpsB;
  qpsB = getQPs(maxordB, qpsB);
  double bL; double bR; double db; double b; int mynumB; ShapeFunction Bshape;

  double tmp_hvalB;
  double u_h;
  u_h = 0.0;
  #pragma omp parallel for default(none),\
  shared(qpsB, flag,a,nelsB,orderB,nodB,xnodB,B_tmp), \
  private(bL,bR,db,b,mynumB,Bshape, tmp_hvalB), \
  reduction(+:u_h)
  for (int Bel=0; Bel<nelsB; Bel++){
    bL = xnodB[nodB[Bel][0]];
    bR = xnodB[nodB[Bel][orderB[Bel]-1]];
    db = (bR-bL)/2.0;
    for (int l2=0; l2<qpsB.nw; l2++){
      b = bL + (1.0+qpsB.xw[l2])*db;
      Bshape = getShapeFuns(qpsB.xw[l2], orderB[Bel], Bshape);
      tmp_hvalB = 0.0;
      for (int k2 = 0; k2<orderB[Bel]; k2++){
        mynumB = nodB[Bel][k2];
        tmp_hvalB = tmp_hvalB + B_tmp[mynumB]*Bshape.psi[k2];
      }

      // complete integration over element
      if (flag == 'x'){
        u_h = u_h + tmp_hvalB*MMS_Source(a,b) * qpsB.w[l2]*db;
      }
      else if (flag == 't'){
        u_h = u_h + tmp_hvalB*MMS_Source(b,a) * qpsB.w[l2]*db;
      }
    }
  }
  return u_h;
}

double GaussianIntegration::MMS_Source(double x, double t){
  return 1.0/v*phi_pt(x,t) - ( 1.0*phi_px(x,t) + D(x)*phi_pxx(x,t) ) + SigAbs(x)*phi_fun(x,t);
}
double GaussianIntegration::phi_fun(double x, double t){
  //return t*cos(t*x)*sin(x);
  //return t*log(x/t)*sin(x);
  //return t*sin(x);
  // return pow(t,x)*sin(x);
}
double GaussianIntegration::phi_px(double x, double t){
  //return t*cos(x)*cos(t*x) - pow(t,2.0)*sin(x)*sin(t*x);
  //return t*cos(x)*log(x/t) + t*sin(x)/x;
  //return t*cos(x);
  // return pow(t,x)*(log(t)*sin(x)+cos(x));
}
double GaussianIntegration::phi_pxx(double x, double t){
  //return -t*cos(t*x)*sin(x) - pow(t,3.0)*cos(t*x)*sin(x) - 2*pow(t,2.0)*cos(x)*sin(t*x);
  //return 2.0*t*cos(x)/x - t*sin(x)/pow(x,2) - t*log(x/t)*sin(x);
  //return -t*sin(x);
  //return -pow(t,x)*sin(x) + pow(t,x)*pow(log(t),2)*sin(x) + 2*pow(t,x)*log(t)*cos(x);
}
double GaussianIntegration::phi_pt(double x, double t){
  //return cos(t*x)*sin(x) - t*x*sin(x)*sin(t*x);
  //return -sin(x) + log(x/t)*sin(x);
  //return sin(x);
  //return x*pow(t,(x-1))*sin(x);
}
double GaussianIntegration::D(double x){
  return (x+1.0);
}
double GaussianIntegration::SigAbs(double x){
  return ((Xend-x)+1.0);
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
