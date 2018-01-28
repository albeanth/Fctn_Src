#include "Mesh_InterfaceXYT.h"
#include <stdlib.h>
#include <omp.h>

using namespace std;

vector<double> GaussianIntegration::Error_Integrate3D(
      int NUMTHREADS,vector< vector< double > > A, vector< vector< double > > B, vector< vector< double > > C,
      int nelsA, vector<int> orderA, vector< vector< int > > nodA, vector< double > xnodA, int maxordA,
      int nelsB, vector<int> orderB, vector< vector< int > > nodB, vector< double > xnodB, int maxordB,
      int nelsC, vector<int> orderC, vector< vector< int > > nodC, vector< double > xnodC, int maxordC){

  omp_set_num_threads(NUMTHREADS);

  int TotNumEnr;
  TotNumEnr = A.size();

  QuadParams qpsA;
  qpsA = getQPs(maxordA, qpsA);
  QuadParams qpsB;
  qpsB = getQPs(maxordB, qpsB);
  QuadParams qpsC;
  qpsC = getQPs(maxordC, qpsC);
  double aL; double aR; double da; double a; int mynumA; ShapeFunction Ashape;
  double bL; double bR; double db; double b; int mynumB; ShapeFunction Bshape;
  double cL; double cR; double dc; double c; int mynumC; ShapeFunction Cshape;
  vector<double> l2Error(TotNumEnr);
  double pgdA; double pgdB; double pgdC;
  double diff; double pgd_int; double pgdEnr_int;
  double tmp;

  for (int enr=0; enr<TotNumEnr; enr++){
    diff = 0.0; pgd_int = 0.0; pgdEnr_int = 0.0;
    vector<double> tmpEnr_pgdA(enr+1); vector<double> tmpEnr_pgdB(enr+1); vector<double> tmpEnr_pgdC(enr+1);

    #pragma omp parallel for default(none),\
    shared(enr, qpsA,qpsB,qpsC, A,B,C, nelsA,orderA,nodA,xnodA, nelsB,orderB,nodB,xnodB, nelsC,orderC,nodC,xnodC), \
    private(aL,aR,da,a,Ashape,pgdA,mynumA, bL,bR,db,b,Bshape,pgdB,mynumB, cL,cR,dc,c,Cshape,pgdC,mynumC, tmp), \
    firstprivate(tmpEnr_pgdA, tmpEnr_pgdB, tmpEnr_pgdC),\
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
                // complete integration over element
                pgdEnr_int = pgdEnr_int + pgdA*pgdB*pgdC  *qpsA.w[l1]*da *qpsB.w[l2]*db *qpsC.w[l3]*dc; //integral of step enr
                tmp = 0.0;
                for (int ide = enr; ide > -1; ide--){
                  pgd_int = pgd_int + tmpEnr_pgdA[ide]*tmpEnr_pgdB[ide]*tmpEnr_pgdC[ide]  *qpsA.w[l1]*da *qpsB.w[l2]*db *qpsC.w[l3]*dc; //total integral PGD up to step enr
                  tmp = tmp + tmpEnr_pgdA[ide]*tmpEnr_pgdB[ide]*tmpEnr_pgdC[ide];
                }
                diff = diff + pow( phi_fun(a,b,c)-tmp, 2.0) *qpsA.w[l1]*da *qpsB.w[l2]*db *qpsC.w[l3]*dc;
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

double GaussianIntegration::Source_Integrate(int NUMTHREADS, char flag, vector<double> B_tmp, vector<double> C_tmp, double a,
                    int nelsB, vector<int> orderB, vector< vector< int > > nodB, vector< double > xnodB, int maxordB,
                    int nelsC, vector<int> orderC, vector< vector< int > > nodC, vector< double > xnodC, int maxordC){

  omp_set_num_threads(NUMTHREADS);

  QuadParams qpsB;
  qpsB = getQPs(maxordB, qpsB);
  QuadParams qpsC;
  qpsC = getQPs(maxordC, qpsC);
  double bL; double bR; double db; double b; int mynumB; ShapeFunction Bshape;
  double cL; double cR; double dc; double c; int mynumC; ShapeFunction Cshape;

  double tmp_hvalB; double tmp_hvalC;
  double u_h;
  u_h = 0.0;
  #pragma omp parallel for default(none),\
  shared(qpsB,qpsC, flag,a,nelsB,orderB,nodB,xnodB,B_tmp, nelsC,orderC,nodC,xnodC,C_tmp), \
  private(bL,bR,db,b,mynumB,Bshape, cL,cR,dc,c,mynumC,Cshape, tmp_hvalB,tmp_hvalC), \
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
          // complete integration over element
          if (flag == 'x'){
            u_h = u_h + tmp_hvalB*tmp_hvalC*MMS_Source(a,b,c) * qpsB.w[l2]*db *qpsC.w[l3]*dc;
          }
          else if (flag == 'y'){
            u_h = u_h + tmp_hvalB*tmp_hvalC*MMS_Source(b,a,c) * qpsB.w[l2]*db *qpsC.w[l3]*dc;
          }
          else if (flag == 't'){
            u_h = u_h + tmp_hvalB*tmp_hvalC*MMS_Source(b,c,a) * qpsB.w[l2]*db *qpsC.w[l3]*dc;
          }
          // else{
          //   cout << "Unknown spatial dimension flag in spatial source. \n Quitting." << endl;
          //   exit(1);
          // }
        }
      }
    }
  }
  return u_h;
}

vector< vector< double > > GaussianIntegration::UpdateR_Serial(
                double b1, double c1, vector<double> B_tmp, vector<double> C_tmp,
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

  vector< vector< double > > r;
  double r1int=0.0; double r1L=0.0; double r1R=0.0;
  double r2int=0.0; double r2L=0.0; double r2R=0.0;
  double r3int=0.0;
  vector<double> r1; vector<double> r2; vector<double> r3; // "1D" vectors. just single entry vectors.

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
          r1int = r1int + D(b,c)*pow(dBhval,2.0)*pow(Chval,2.0)  *qpsB.w[l2]*db *qpsC.w[l3]*dc;
          r2int = r2int + D(b,c)*pow(Bhval,2.0)*pow(dChval,2.0)  *qpsB.w[l2]*db *qpsC.w[l3]*dc;
          r3int = r3int + SigAbs(b,c)*pow(Bhval,2.0)*pow(Chval,2.0)*qpsB.w[l2]*db *qpsC.w[l3]*dc;
          if ((Bel==0) && (l2==0) && (Cel==0) && (l3==0)){
            r1L = Bhval*D(b,c)*dBhval;
            r2L = Chval*D(b,c)*dChval;
          }
          else if ((Bel==nelsB-1) && (l2==0) && (Cel==0) && (l3==0)){
            r1R = Bhval*D(b,c)*dBhval;
            r2R = Chval*D(b,c)*dChval;
          }
        }
      }
    }
  }
  r1.push_back( (r1R - r1L)*c1 - r1int );
  r2.push_back( (r2R - r2L)*b1 - r2int );
  r3.push_back( r3int );
  r.push_back(r1); r.push_back(r2); r.push_back(r3);
  return r;
}

vector< vector< double > > GaussianIntegration::UpdateR_Enr_Serial(
                vector<double> b2, vector<double> c2, vector<double> B_tmp, vector<double> C_tmp, vector< vector< double > > B, vector< vector< double > > C,
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
  vector< vector< double > > r;
  double r4int; double r4L; double r4R;
  double r5int; double r5L; double r5R;
  vector<double> r4(TotNumEnr); 
  vector<double> r5(TotNumEnr);
  vector<double> r6(TotNumEnr); fill(r6.begin(),r6.end(),0.0);

  for (int enr=0; enr<TotNumEnr; enr++){
    r4int=0.0; r4L=0.0; r4R=0.0;
    r5int=0.0; r5L=0.0; r5R=0.0;
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
          Bhval += B_tmp[mynumB]*Bshape.psi[k2]; dBhval += B_tmp[mynumB]*Bshape.dpsi[k2]/db;
          Bihval += B[enr][mynumB]*Bshape.psi[k2]; dBihval += B[enr][mynumB]*Bshape.dpsi[k2]/db;
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
              Cihval += C[enr][mynumC]*Cshape.psi[k3]; dCihval += C[enr][mynumC]*Cshape.dpsi[k3]/dc;
              Chval += C_tmp[mynumC]*Cshape.psi[k3]; dChval += C_tmp[mynumC]*Cshape.dpsi[k3]/dc;
            }

            // complete integration over element
            r4int = r4int + D(b,c)*dBhval*dBihval*Chval*Cihval  *qpsB.w[l2]*db *qpsC.w[l3]*dc;
            r5int = r5int + D(b,c)*Bhval*Bihval*dChval*dCihval  *qpsB.w[l2]*db *qpsC.w[l3]*dc;
            r6[enr] = r6[enr] + SigAbs(b,c)*Bhval*Bihval*Chval*Cihval *qpsB.w[l2]*db *qpsC.w[l3]*dc;
            if ((Bel==0) && (l2==0) && (Cel==0) && (l3==0)) {
              r4L = Bhval*D(b,c)*dBihval;
              r5L = Chval*D(b,c)*dCihval;
            }
            else if ((Bel==nelsB-1) && (l2==0) && (Cel==0) && (l3==0)){
              r4R = Bhval*D(b,c)*dBihval;
              r5R = Chval*D(b,c)*dCihval;
            }
          }
        }
      }
    }
    r4[enr] = (r4R-r4L)*c2[enr] - r4int;
    r5[enr] = (r5R-r5L)*b2[enr] - r5int;
  }
r.push_back(r4); r.push_back(r5); r.push_back(r6);
return r;
}

vector< vector< double > > GaussianIntegration::UpdateR_Parallel(
                int NUMTHREADS, double b1, double c1, vector<double> B_tmp, vector<double> C_tmp,
                int nelsB, vector<int> orderB, vector< vector< int > > nodB, vector< double > xnodB, int maxordB,
                int nelsC, vector<int> orderC, vector< vector< int > > nodC, vector< double > xnodC, int maxordC){

  omp_set_num_threads(NUMTHREADS);

  QuadParams qpsB;
  qpsB = getQPs(maxordB, qpsB);
  QuadParams qpsC;
  qpsC = getQPs(maxordC, qpsC);
  double bL; double bR; double db; double b; int mynumB; ShapeFunction Bshape;
  double Bhval; double dBhval;
  double cL; double cR; double dc; double c; int mynumC; ShapeFunction Cshape;
  double Chval; double dChval;

  vector< vector< double > > r;
  double r1int=0.0; double r1L=0.0; double r1R=0.0;
  double r2int=0.0; double r2L=0.0; double r2R=0.0;
  double r3int=0.0;
  vector<double> r1; vector<double> r2; vector<double> r3; // "1D" vectors. just single entry vectors.
  #pragma omp parallel for default(none),\
  shared(qpsB,qpsC, nelsB,orderB,nodB,xnodB,B_tmp, nelsC,orderC,nodC,xnodC,C_tmp), \
  private(bL,bR,db,b,mynumB,Bshape,Bhval,dBhval, cL,cR,dc,c,mynumC,Cshape,Chval,dChval), \
  reduction(+:r1int,r1L,r1R,r2int,r2L,r2R,r3int)
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
          r1int = r1int + D(b,c)*pow(dBhval,2.0)*pow(Chval,2.0)  *qpsB.w[l2]*db *qpsC.w[l3]*dc;
          r2int = r2int + D(b,c)*pow(Bhval,2.0)*pow(dChval,2.0)  *qpsB.w[l2]*db *qpsC.w[l3]*dc;
          r3int = r3int + SigAbs(b,c)*pow(Bhval,2.0)*pow(Chval,2.0)*qpsB.w[l2]*db *qpsC.w[l3]*dc;
          if ((Bel==0) && (l2==0) && (Cel==0) && (l3==0)){
            r1L = Bhval*D(b,c)*dBhval;
            r2L = Chval*D(b,c)*dChval;
          }
          else if ((Bel==nelsB-1) && (l2==0) && (Cel==0) && (l3==0)){
            r1R = Bhval*D(b,c)*dBhval;
            r2R = Chval*D(b,c)*dChval;
          }
        }
      }
    }
  }
  r1.push_back( (r1R - r1L)*c1 - r1int );
  r2.push_back( (r2R - r2L)*b1 - r2int );
  r3.push_back( r3int );
  r.push_back(r1); r.push_back(r2); r.push_back(r3);
  return r;
}

vector< vector< double > > GaussianIntegration::UpdateR_Enr_Parallel(
		int NUMTHREADS, vector<double> b2, vector<double> c2, vector<double> B_tmp, vector<double> C_tmp, vector< vector< double > > B, vector< vector< double > > C,
                int nelsB, vector<int> orderB, vector< vector< int > > nodB, vector< double > xnodB, int maxordB,
                int nelsC, vector<int> orderC, vector< vector< int > > nodC, vector< double > xnodC, int maxordC){

  omp_set_num_threads(NUMTHREADS);

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
  vector< vector< double > > r;
  double r4int; double r4L; double r4R;
  double r5int; double r5L; double r5R;
  double r6tmp;
  vector<double> r4(TotNumEnr); vector<double> r5(TotNumEnr);
  vector<double> r6(TotNumEnr); fill(r6.begin(),r6.end(),0.0);

  for (int enr=0; enr<TotNumEnr; enr++){
    r4int=0.0; r4L=0.0; r4R=0.0;
    r5int=0.0; r5L=0.0; r5R=0.0;
    r6tmp = 0.0;
    #pragma omp parallel for default(none),\
    shared(enr, qpsB,qpsC, B_tmp,B,C_tmp,C, nelsB,orderB,nodB,xnodB, nelsC,orderC,nodC,xnodC), \
    private(bL,bR,db,b,Bshape,mynumB,Bhval,dBhval,Bihval,dBihval, cL,cR,dc,c,Cshape,mynumC,Chval,dChval,Cihval,dCihval), \
    reduction(+:r4int,r4L,r4R, r5int,r5L,r5R, r6tmp)
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
          Bhval += B_tmp[mynumB]*Bshape.psi[k2]; dBhval += B_tmp[mynumB]*Bshape.dpsi[k2]/db;
          Bihval += B[enr][mynumB]*Bshape.psi[k2]; dBihval += B[enr][mynumB]*Bshape.dpsi[k2]/db;
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
              Cihval += C[enr][mynumC]*Cshape.psi[k3]; dCihval += C[enr][mynumC]*Cshape.dpsi[k3]/dc;
              Chval += C_tmp[mynumC]*Cshape.psi[k3]; dChval += C_tmp[mynumC]*Cshape.dpsi[k3]/dc;
            }

            // complete integration over element
            r4int = r4int + D(b,c)*dBhval*dBihval*Chval*Cihval  *qpsB.w[l2]*db *qpsC.w[l3]*dc;
            r5int = r5int + D(b,c)*Bhval*Bihval*dChval*dCihval  *qpsB.w[l2]*db *qpsC.w[l3]*dc;
            r6tmp = r6tmp + SigAbs(b,c)*Bhval*Bihval*Chval*Cihval *qpsB.w[l2]*db *qpsC.w[l3]*dc;
            if ((Bel==0) && (l2==0) && (Cel==0) && (l3==0)) {
              r4L = Bhval*D(b,c)*dBihval;
              r5L = Chval*D(b,c)*dCihval;
            }
            else if ((Bel==nelsB-1) && (l2==0) && (Cel==0) && (l3==0)){
              r4R = Bhval*D(b,c)*dBihval;
              r5R = Chval*D(b,c)*dCihval;
            }
          }
        }
      }
    }
    r4[enr] = (r4R-r4L)*c2[enr] - r4int;
    r5[enr] = (r5R-r5L)*b2[enr] - r5int;
    r6[enr] = r6tmp;
  }
r.push_back(r4); r.push_back(r5); r.push_back(r6);
return r;
}


double GaussianIntegration::MMS_Source(double x, double y, double t){
  return 1/v*phi_pt(x,y,t) - (( D_px(y)*phi_px(x,y,t) + D(x,y)*phi_pxx(x,y,t) ) + ( D_py(x)*phi_py(x,y,t) + D(x,y)*phi_pyy(x,y,t)) ) + SigAbs(x,y)*phi_fun(x,y,t);
}
double GaussianIntegration::phi_fun(double x, double y, double t){
  //return t*sin(x)*sin(y);
  return t*sin(pow(x,2))*sin(pow(y,2));
  // return pow(t,x)*sin(x)*sin(y);
}
double GaussianIntegration::phi_px(double x, double y, double t){
  //return t*cos(x)*sin(y);
  return 2*x*t*cos(pow(x,2))*sin(pow(y,2));
  //return pow(t,x)*cos(x)*sin(y) + pow(t,x)*log(t)*sin(x)*sin(y);
}
double GaussianIntegration::phi_pxx(double x, double y, double t){
  //return -t*sin(x)*sin(y);
  return 2*t*sin(pow(y,2))*(cos(pow(x,2))-2*(pow(x,2))*sin(pow(x,2)));
  //return sin(y)*(2*pow(t,x)*cos(x)*log(t) - pow(t,x)*sin(x) + pow(t,x)*pow(log10(t),2.0)*sin(x));
}
double GaussianIntegration::phi_py(double x, double y, double t){
  //return t*sin(x)*cos(y);
  return 2*y*t*sin(pow(x,2))*cos(pow(y,2));
  //return pow(t,x)*cos(y)*sin(x);
}
double GaussianIntegration::phi_pyy(double x, double y, double t){
  //return -t*sin(x)*sin(y);
  return 2*t*sin(pow(x,2))*( cos(pow(y,2))-2*(pow(y,2))*sin(pow(y,2)) );
  //return -pow(t,x)*sin(x)*sin(y);
}
double GaussianIntegration::phi_pt(double x, double y, double t){
  //return sin(x)*sin(y);
  return sin(pow(x,2))*sin(pow(y,2));
  //return pow(t,x-1.0)*x*sin(x)*sin(y);
}
double GaussianIntegration::D(double x, double y){
  return (x*y) + 1.0;
}
double GaussianIntegration::D_px(double y){
  return y;
}
double GaussianIntegration::D_py(double x){
  return x;
}
double GaussianIntegration::SigAbs(double x, double y){
  return ((Xbnd-x)*(Ybnd-y)) + 1.0;
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
