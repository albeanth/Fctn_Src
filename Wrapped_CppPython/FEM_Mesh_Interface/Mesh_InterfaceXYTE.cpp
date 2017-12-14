#include "Mesh_InterfaceXYTE.h"
#include <stdlib.h>
#include <math.h>
#include <iostream>

using namespace std;

vector<double> GaussianIntegration::Error_Integrate4D(
      vector< vector< double > > A, vector< vector< double > > B, vector< vector< double > > C, vector< vector< double > > D,
      int nelsA, vector<int> orderA, vector< vector< int > > nodA, vector< double > xnodA, int maxordA,
      int nelsB, vector<int> orderB, vector< vector< int > > nodB, vector< double > xnodB, int maxordB,
      int nelsC, vector<int> orderC, vector< vector< int > > nodC, vector< double > xnodC, int maxordC,
      int nelsD, vector<int> orderD, vector< vector< int > > nodD, vector< double > xnodD, int maxordD){

  int TotNumEnr;
  TotNumEnr = A.size();

  QuadParams qpsA;
  qpsA = getQPs(maxordA, qpsA);
  QuadParams qpsB;
  qpsB = getQPs(maxordB, qpsB);
  QuadParams qpsC;
  qpsC = getQPs(maxordC, qpsC);
  QuadParams qpsD;
  qpsD = getQPs(maxordD, qpsD);
  double aL; double aR; double da; double a; int mynumA; ShapeFunction Ashape;
  double bL; double bR; double db; double b; int mynumB; ShapeFunction Bshape;
  double cL; double cR; double dc; double c; int mynumC; ShapeFunction Cshape;
  double dL; double dR; double dd; double d; int mynumD; ShapeFunction Dshape;
  vector<double> l2Error(TotNumEnr);
  double pgdA; double pgdB; double pgdC; double pgdD;
  double diff; double pgd_int; double pgdEnr_int;
  double tmp;
  printf("    on enrichment number...\n");
  for (int enr=0; enr<TotNumEnr; enr++){
    diff = 0.0; pgd_int = 0.0; pgdEnr_int = 0.0;
    vector<double> tmpEnr_pgdA(enr+1); vector<double> tmpEnr_pgdB(enr+1); vector<double> tmpEnr_pgdC(enr+1); vector<double> tmpEnr_pgdD(enr+1);

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

            for (int Cel=0; Cel<nelsC; Cel++){
              cL = xnodC[nodC[Cel][0]];
              cR = xnodC[nodC[Cel][orderC[Cel]-1]];
              dc = (cR-cL)/2.0;
              for (int l3=0; l3<qpsC.nw; l3++){
                c = cL + (1.0+qpsC.xw[l3])*dc;
                Cshape = getShapeFuns(qpsC.xw[l3], orderC[Cel], Cshape);
                pgdC = 0.0; fill(tmpEnr_pgdC.begin(),tmpEnr_pgdC.end(),0.0);
                for (int k3 = 0; k3<orderC[Cel]; k3++){
                  mynumC = nodC[Cel][k3];
                  pgdC = pgdC + C[enr][mynumC]*Cshape.psi[k3];
                  for (int ide = enr; ide > -1; ide--){
                    tmpEnr_pgdC[ide] = tmpEnr_pgdC[ide] + C[ide][mynumC]*Cshape.psi[k3];
                  }
                }

                for (int Del=0; Del<nelsD; Del++){
                  dL = xnodD[nodD[Del][0]];
                  dR = xnodD[nodD[Del][orderD[Del]-1]];
                  dd = (dR-dL)/2.0;
                  for (int l4=0; l4<qpsD.nw; l4++){
                    d = dL + (1.0+qpsD.xw[l4])*dd;
                    Dshape = getShapeFuns(qpsD.xw[l4], orderD[Del], Dshape);
                    pgdD = 0.0; fill(tmpEnr_pgdD.begin(),tmpEnr_pgdD.end(),0.0);
                    for (int k4 = 0; k4<orderD[Del]; k4++){
                      mynumD = nodD[Del][k4];
                      pgdD = pgdD + D[enr][mynumD]*Dshape.psi[k4];
                      for (int ide = enr; ide > -1; ide--){
                        tmpEnr_pgdD[ide] = tmpEnr_pgdD[ide] + D[ide][mynumD]*Dshape.psi[k4];
                      }
                    }

                    // complete integration over element
                    pgdEnr_int = pgdEnr_int + pgdA*pgdB*pgdC*pgdD  *qpsA.w[l1]*da *qpsB.w[l2]*db *qpsC.w[l3]*dc *qpsD.w[l4]*dd; //integral of step enr
                    tmp = 0.0;
                    for (int ide = enr; ide > -1; ide--){
                      pgd_int = pgd_int + tmpEnr_pgdA[ide]*tmpEnr_pgdB[ide]*tmpEnr_pgdC[ide]*tmpEnr_pgdD[ide]  *qpsA.w[l1]*da *qpsB.w[l2]*db *qpsC.w[l3]*dc *qpsD.w[l4]*dd; //total integral PGD up to step enr
                      tmp = tmp + tmpEnr_pgdA[ide]*tmpEnr_pgdB[ide]*tmpEnr_pgdC[ide]*tmpEnr_pgdD[ide];
                    }
                    diff = diff + pow( phi_fun(a,b,c,d)-tmp, 2.0) *qpsA.w[l1]*da *qpsB.w[l2]*db *qpsC.w[l3]*dc *qpsD.w[l4]*dd;
                  }
                }
              }
            }
          }
        }
      }
    }
    l2Error[enr] = sqrt(diff) ;
    printf("    %d, L2Error = %.6e, PGD Int = %.7f, Tot PGD Int = %.7f\n",enr+1,l2Error[enr],pgdEnr_int,pgd_int);
  }
  return l2Error;
}

double GaussianIntegration::Source_Integrate(
                char flag, vector<double> B_tmp, vector<double> C_tmp, vector<double> D_tmp, double a,
                int nelsB, vector<int> orderB, vector< vector< int > > nodB, vector< double > xnodB, int maxordB,
                int nelsC, vector<int> orderC, vector< vector< int > > nodC, vector< double > xnodC, int maxordC,
                int nelsD, vector<int> orderD, vector< vector< int > > nodD, vector< double > xnodD, int maxordD){
  QuadParams qpsB;
  qpsB = getQPs(maxordB, qpsB);
  QuadParams qpsC;
  qpsC = getQPs(maxordC, qpsC);
  QuadParams qpsD;
  qpsD = getQPs(maxordD, qpsD);
  double bL; double bR; double db; double b; int mynumB; ShapeFunction Bshape;
  double cL; double cR; double dc; double c; int mynumC; ShapeFunction Cshape;
  double dL; double dR; double dd; double d; int mynumD; ShapeFunction Dshape;

  double tmp_hvalB; double tmp_hvalC; double tmp_hvalD;
  double u_h;
  u_h = 0.0;
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

      for (int Cel=0; Cel<nelsC; Cel++){
        cL = xnodC[nodC[Cel][0]];
        cR = xnodC[nodC[Cel][orderC[Cel]-1]];
        dc = (cR-cL)/2.0;
        for (int l3=0; l3<qpsC.nw; l3++){
          c = cL + (1.0+qpsC.xw[l3])*dc;
          Cshape = getShapeFuns(qpsC.xw[l3], orderC[Cel], Cshape);
          tmp_hvalC = 0.0;
          for (int k3 = 0; k3<orderC[Cel]; k3++){
            mynumC = nodC[Cel][k3];
            tmp_hvalC = tmp_hvalC + C_tmp[mynumC]*Cshape.psi[k3];
          }

          for (int Del=0; Del<nelsD; Del++){
            dL = xnodD[nodD[Del][0]];
            dR = xnodD[nodD[Del][orderD[Del]-1]];
            dd = (dR-dL)/2.0;
            for (int l4=0; l4<qpsD.nw; l4++){
              d = dL + (1.0+qpsD.xw[l4])*dd;
              Dshape = getShapeFuns(qpsD.xw[l4], orderD[Del], Dshape);
              tmp_hvalD = 0.0;
              for (int k4 = 0; k4<orderD[Del]; k4++){
                mynumD = nodD[Del][k4];
                tmp_hvalD = tmp_hvalD + D_tmp[mynumD]*Dshape.psi[k4];
              }

              // complete integration over element
              if (flag == 'x'){
                u_h = u_h + tmp_hvalB*tmp_hvalC*MMS_Source(a,b,c,d) * qpsB.w[l2]*db *qpsC.w[l3]*dc;
              }
              else if (flag == 'y'){
                u_h = u_h + tmp_hvalB*tmp_hvalC*MMS_Source(b,a,c,d) * qpsB.w[l2]*db *qpsC.w[l3]*dc;
              }
              else if (flag == 't'){
                u_h = u_h + tmp_hvalB*tmp_hvalC*MMS_Source(b,c,a,d) * qpsB.w[l2]*db *qpsC.w[l3]*dc;
              }
              else if (flag == 'e'){
                u_h = u_h + tmp_hvalB*tmp_hvalC*MMS_Source(b,c,d,a) * qpsB.w[l2]*db *qpsC.w[l3]*dc;
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
  return u_h;
}

vector< vector< double > > GaussianIntegration::UpdateH(
                double a, double b1, double c1, vector<double> B_tmp, vector<double> C_tmp,
                int nelsB, vector<int> orderB, vector< vector< int > > nodB, vector< double > xnodB, int maxordB,
                int nelsC, vector<int> orderC, vector< vector< int > > nodC, vector< double > xnodC, int maxordC){

  QuadParams qpsB;
  qpsB = getQPs(maxordB, qpsB);
  QuadParams qpsC;
  qpsC = getQPs(maxordC, qpsC);
  double bL; double bR; double db; double b; int mynumB; ShapeFunction Bshape;
  double Bhval; double dBhval;
  double cL; double cR; double dc; double c; int mynumC; ShapeFunction Cshape;
  double Chval; double dChval;

  vector< vector< double > > h;
  double h1int; double h1L; double h1R;
  double h2int; double h2L; double h2R;
  double h3int;
  vector<double> h1; vector<double> h2; vector<double> h3;

  for (int Bel=0; Bel<nelsB; Bel++){
    bL = xnodB[nodB[Bel][0]];
    bR = xnodB[nodB[Bel][orderB[Bel]-1]];
    db = (bR-bL)/2.0;
    for (int l2=0; l2<qpsB.nw; l2++){
      b = bL + (1.0+qpsB.xw[l2])*db;
      Bshape = getShapeFuns(qpsB.xw[l2], orderB[Bel], Bshape);
      Bhval = 0.0; dBhval = 0.0;
      for (int k2 = 0; k2<orderB[Bel]; k2++){
        mynumB = nodB[Bel][k2];
        Bhval = Bhval + B_tmp[mynumB]*Bshape.psi[k2];
        dBhval = dBhval + B_tmp[mynumB]*Bshape.dpsi[k2]/db;
      }

      for (int Cel=0; Cel<nelsC; Cel++){
        cL = xnodC[nodC[Cel][0]];
        cR = xnodC[nodC[Cel][orderC[Cel]-1]];
        dc = (cR-cL)/2.0;
        for (int l3=0; l3<qpsC.nw; l3++){
          c = cL + (1.0+qpsC.xw[l3])*dc;
          Cshape = getShapeFuns(qpsC.xw[l3], orderC[Cel], Cshape);
          Chval = 0.0; dChval = 0.0;
          for (int k3 = 0; k3<orderC[Cel]; k3++){
            mynumC = nodC[Cel][k3];
            Chval = Chval + C_tmp[mynumC]*Cshape.psi[k3];
            dChval = dChval + C_tmp[mynumC]*Cshape.dpsi[k3]/dc;
          }

          // complete integration over element
          h1int = h1int + D(b,c,a)*pow(dBhval,2.0)*pow(Chval,2.0)  *qpsB.w[l2]*db *qpsC.w[l3]*dc;
          h2int = h2int + D(b,c,a)*pow(Bhval,2.0)*pow(dChval,2.0)  *qpsB.w[l2]*db *qpsC.w[l3]*dc;
          h3int = h3int + SigT(b,c,a)*pow(Bhval,2.0)*pow(Chval,2.0)*qpsB.w[l2]*db *qpsC.w[l3]*dc;
          if ((Bel==0) && (l2==0) && (Cel==0) && (l3==0)){
            h1L = Bhval*D(b,c,a)*dBhval;
            h2L = Chval*D(b,c,a)*dChval;
          }
          else if ((Bel==nelsB-1) && (l2==0) && (Cel==0) && (l3==0)){
            h1R = Bhval*D(b,c,a)*dBhval;
            h2R = Chval*D(b,c,a)*dChval;
          }
        }
      }
    }
  }
  h1.push_back( (h1R - h1L)*c1 - h1int );
  h2.push_back( (h2R - h2L)*b1 - h2int );
  h3.push_back( h3int );
  h.push_back(h1); h.push_back(h2); h.push_back(h3);
  return(h);
}

vector< vector< double > > GaussianIntegration::UpdateH_Enr(
                double a, double b2, double c2, vector<double> B_tmp, vector<double> C_tmp, vector< vector< double > > B, vector< vector< double > > C,
                int nelsB, vector<int> orderB, vector< vector< int > > nodB, vector< double > xnodB, int maxordB,
                int nelsC, vector<int> orderC, vector< vector< int > > nodC, vector< double > xnodC, int maxordC){

  QuadParams qpsB;
  qpsB = getQPs(maxordB, qpsB);
  QuadParams qpsC;
  qpsC = getQPs(maxordC, qpsC);
  double bL; double bR; double db; double b; int mynumB; ShapeFunction Bshape;
  double Bhval; double dBhval; double Bihval; double dBihval;
  double cL; double cR; double dc; double c; int mynumC; ShapeFunction Cshape;
  double Chval; double dChval; double Cihval; double dCihval;

  int TotNumEnr;
  TotNumEnr = B.size();
  vector< vector< double > > h;
  double h4int; double h4L; double h4R;
  double h5int; double h5L; double h5R;
  vector<double> h4(TotNumEnr); vector<double> h5(TotNumEnr);
  vector<double> h6(TotNumEnr); fill(h6.begin(),h6.end(),0.0);

  for (int enr=0; enr<TotNumEnr; enr++){
    h4int=0.0; h4L=0.0; h4R=0.0;
    h5int=0.0; h5L=0.0; h5R=0.0;
    for (int Bel=0; Bel<nelsB; Bel++){
      bL = xnodB[nodB[Bel][0]];
      bR = xnodB[nodB[Bel][orderB[Bel]-1]];
      db = (bR-bL)/2.0;
      for (int l2=0; l2<qpsB.nw; l2++){
        b = bL + (1.0+qpsB.xw[l2])*db;
        Bshape = getShapeFuns(qpsB.xw[l2], orderB[Bel], Bshape);
        Bhval = 0.0; dBhval = 0.0; Bihval = 0.0; dBihval = 0.0;
        for (int k2 = 0; k2<orderB[Bel]; k2++){
          mynumB = nodB[Bel][k2];
          Bhval = Bhval + B_tmp[mynumB]*Bshape.psi[k2]; dBhval = dBhval + B_tmp[mynumB]*Bshape.dpsi[k2]/db;
          Bihval = Bihval + B[enr][mynumB]*Bshape.psi[k2]; dBihval = dBihval + B[enr][mynumB]*Bshape.dpsi[k2]/db;
        }

        for (int Cel=0; Cel<nelsC; Cel++){
          cL = xnodC[nodC[Cel][0]];
          cR = xnodC[nodC[Cel][orderC[Cel]-1]];
          dc = (cR-cL)/2.0;
          for (int l3=0; l3<qpsC.nw; l3++){
            c = cL + (1.0+qpsC.xw[l3])*dc;
            Cshape = getShapeFuns(qpsC.xw[l3], orderC[Cel], Cshape);
            Chval = 0.0; dChval = 0.0; Cihval = 0.0; dCihval = 0.0;
            for (int k3 = 0; k3<orderC[Cel]; k3++){
              mynumC = nodC[Cel][k3];
              Cihval = Chval + C[enr][mynumC]*Cshape.psi[k3]; dCihval = dChval + C[enr][mynumC]*Cshape.dpsi[k3]/dc;
              Chval = Chval + C_tmp[mynumC]*Cshape.psi[k3]; dChval = dChval + C_tmp[mynumC]*Cshape.dpsi[k3]/dc;
            }

            // complete integration over element
            h4int = h4int + D(b,c,a)*dBhval*dBihval*Chval*Cihval  *qpsB.w[l2]*db *qpsC.w[l3]*dc;
            h5int = h5int + D(b,c,a)*Bhval*Bihval*dChval*dCihval  *qpsB.w[l2]*db *qpsC.w[l3]*dc;
            h6[enr] = h6[enr] + SigT(b,c,a)*Bhval*Bihval*Chval*Cihval *qpsB.w[l2]*db *qpsC.w[l3]*dc;
            if ((Bel==0) && (l2==0) && (Cel==0) && (l3==0)) {
              h4L = Bhval*D(b,c,a)*dBihval;
              h5L = Chval*D(b,c,a)*dCihval;
            }
            else if ((Bel==nelsB-1) && (l2==0) && (Cel==0) && (l3==0)){
              h4R = Bhval*D(b,c,a)*dBihval;
              h5R = Chval*D(b,c,a)*dCihval;
            }
          }
        }
      }
    }
    h4[enr] = (h4R-h4L)*c2 - h4int;
    h5[enr] = (h5R-h5L)*b2 - h5int;
  }
h.push_back(h4); h.push_back(h5); h.push_back(h6);
return h;
}

vector< vector< double > > UpdateG(
        vector<double> A_tmp, double x1, int nelsA, vector<int> orderA, vector< vector< int > > nodA, vector< double > xnodA, int maxordA,
        vector<double> B_tmp, double y1, int nelsB, vector<int> orderB, vector< vector< int > > nodB, vector< double > xnodB, int maxordB,
        vector<double> C_tmp, double c1, int nelsC, vector<int> orderC, vector< vector< int > > nodC, vector< double > xnodC, int maxordC);





double GaussianIntegration::MMS_Source(double x, double y, double t, double E){
  return 1/V(E)*phi_pt(x,y,E) - (( D_px(y,E)*phi_px(x,y,t,E) + D(x,y,E)*phi_pxx(x,y,t,E) ) + ( D_py(x,E)*phi_py(x,y,t,E) + D(x,y,E)*phi_pyy(x,y,t,E)) ) + SigT(x,y,E)*phi_fun(x,y,t,E);
}
double GaussianIntegration::Q(double E){
  return fabs(log10(E)) + 1.0;
}
double GaussianIntegration::phi_fun(double x, double y, double t, double E){
  return Q(E)*t*sin(x)*sin(y);
}
double GaussianIntegration::phi_px(double x, double y, double t, double E){
  return Q(E)*t*cos(x)*sin(y);
}
double GaussianIntegration::phi_pxx(double x, double y, double t, double E){
  return Q(E)*-t*sin(x)*sin(y);
}
double GaussianIntegration::phi_py(double x, double y, double t, double E){
  return Q(E)*t*sin(x)*cos(y);
}
double GaussianIntegration::phi_pyy(double x, double y, double t, double E){
  return Q(E)*-t*sin(x)*sin(y);
}
double GaussianIntegration::phi_pt(double x, double y, double E){
  return Q(E)*sin(x)*sin(y);
}
double GaussianIntegration::D(double x, double y, double E){
  return Q(E)*((x*y) + 1.0);
}
double GaussianIntegration::D_px(double y, double E){
  return Q(E)*y;
}
double GaussianIntegration::D_py(double x, double E){
  return Q(E)*x;
}
double GaussianIntegration::SigT(double x, double y, double E){
  return Q(E)*((Xbnd-x)*(Ybnd-y)) + 1.0;
}

double GaussianIntegration::V(double E){
  return Q(E);
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
  else{
    cout << "Order = " << n << " shape function does not exist." << endl;
    exit(1);
  }
  shape.psi = y;
  shape.dpsi = dy;
  return shape;
}
