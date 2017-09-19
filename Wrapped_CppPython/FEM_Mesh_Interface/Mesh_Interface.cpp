#include "Mesh_Interface.h"
#include <stdlib.h>
#include <math.h>

using namespace std;

std::vector<double> GaussianIntegration::Source_Integrate(char flag,
                    int nelsA, std::vector<int> orderA, std::vector< std::vector< int > > nodA, std::vector< double > xnodA, int maxordA,
                    int nelsB, std::vector<int> orderB, std::vector< std::vector< int > > nodB, std::vector< double > xnodB, int maxordB,
                    int nelsC, std::vector<int> orderC, std::vector< std::vector< int > > nodC, std::vector< double > xnodC, int maxordC){
  QuadParams qpsA;
  qpsA = getQPs(maxordA, qpsA);
  QuadParams qpsB;
  qpsB = getQPs(maxordB, qpsB);
  QuadParams qpsC;
  qpsC = getQPs(maxordC, qpsC);
  std::vector<double> srcA;
  srcA.resize(nelsA);
  double aL; double aR; double da; double a;
  double bL; double bR; double db; double b;
  double cL; double cR; double dc; double c;
  for(int Ael=0; Ael<nelsA; Ael++){
    aL = xnodA[nodA[Ael][0]];
    aR = xnodA[nodA[Ael][orderA[Ael]-1]];
    da = (aR-aL)/2.0;
    for (int l1=0; l1<qpsA.nw; l1++){
      a = aL + (1.0+qpsA.xw[l1])*da;

      for (int Bel=0; Bel<nelsB; Bel++){
        bL = xnodB[nodB[Bel][0]];
        bR = xnodB[nodB[Bel][orderB[Bel]-1]];
        db = (bR-bL)/2.0;
        for (int l2=0; l2<qpsB.nw; l2++){
          b = bL + (1.0+qpsB.xw[l2])*db;

          for(int Cel=0; Cel<nelsC; Cel++){
            cL = xnodC[nodC[Cel][0]];
            cR = xnodC[nodC[Cel][orderC[Cel]-1]];
            dc = (cR-cL)/2.0;
            // printf("%.3f, %.3f -> %.3f\n", cL,cR,dc);
            for (int l3=0; l1<qpsC.nw; l1++){
              c = cL + (1.0+qpsC.xw[l3])*dc;

              // complete integration over element
              if (flag == 'x'){
                srcA[Ael+l1] = srcA[Ael+l1] + MMS_Source(a,b,c) * qpsB.w[l2]*db *qpsC.w[l3]*dc;
              }
              else if (flag == 'y'){
                srcA[Ael+l1] = srcA[Ael+l1] + MMS_Source(b,a,c) * qpsB.w[l2]*db *qpsC.w[l3]*dc;
              }
              else if (flag == 't'){
                srcA[Ael+l1] = srcA[Ael+l1] + MMS_Source(b,c,a) * qpsB.w[l2]*db *qpsC.w[l3]*dc;
              }
              else{
                cout << "Unknown spatial dimension flag in spatial source. \n Quitting." << endl;
                exit(1);
              }
            }
          }
        }
      }
    }
  }
  return srcA;
}

double GaussianIntegration::MMS_Source(double x, double y, double t){
  return 1/v*phi_pt(x,y) - (( 1.0*phi_px(x,y,t) + D(x,y)*phi_pxx(x,y,t) ) + ( 1.0*phi_py(x,y,t) + D(x,y)*phi_pyy(x,y,t)) ) + SigAbs(x,y)*phi_fun(x,y,t);
}
double GaussianIntegration::phi_fun(double x, double y, double t){
  return t*cos(x)*cos(y);
}
double GaussianIntegration::phi_px(double x, double y, double t){
  return -t*sin(x)*cos(y);
}
double GaussianIntegration::phi_pxx(double x, double y, double t){
  return -t*cos(x)*cos(y);
}
double GaussianIntegration::phi_py(double x, double y, double t){
  return -t*cos(x)*sin(y);
}
double GaussianIntegration::phi_pyy(double x, double y, double t){
  return -t*cos(x)*cos(y);
}
double GaussianIntegration::phi_pt(double x, double y){
  return cos(x)*cos(y);
}
double GaussianIntegration::D(double x, double y){
  return x+y+1;
}
double GaussianIntegration::SigAbs(double x, double y){
  return x+y+1;
}

void GaussianIntegration::Get2DInfo(int nels, std::vector<int> order, std::vector< std::vector< int > > nod, std::vector< double > xnod, int maxord){
  cout << "print nels" << endl << nels <<endl;

  cout << "print order" << endl;
  for (unsigned long i=0; i<order.size(); i++){
    printf("%d\n", order[i]);
  }

  cout << "print nod" << endl;
  for (unsigned long i=0; i<nod.size(); i++){
    for (unsigned long j = 0; j<nod[i].size(); j++){
      printf("%d, ",nod[i][j]);
    }
    printf("\n");
  }

  cout << "print xnod" << endl;
  for (unsigned long i=0; i<xnod.size(); i++){
    printf("%f\n", xnod[i]);
  }

  cout << "print maxord" << endl << maxord << endl;
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
