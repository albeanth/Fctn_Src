#include "UtilityFunctions.hpp"

//==========================================================================//
//                      General Purpose operations                          //
//==========================================================================//
std::vector<double> UtilityFunctions::linspace(const double a, const double b, const int NumPts){
  std::vector<double> x;
  double dx { (b-a)/(static_cast<double>(NumPts-1)) };
  for (int i=0; i<NumPts; i++){
    x.push_back( a + dx*static_cast<double>(i) );
  }
  return x;
}

std::vector<int> UtilityFunctions::linspace(const int a, const int b){
  std::vector<int> x;
  for (int i=a; i<b; i++){
    x.push_back( i );
  }
  return x;
}

std::vector<double> UtilityFunctions::logspace(const double a, const double b, const int NumPts){
  std::vector<double> x, xout;
  x = linspace(a,b,NumPts);
  for (double val : x){
    xout.push_back(pow(10.0,val));
  }
  return xout;
}

std::vector<double> UtilityFunctions::arange(const double a, const double b, const double dx){
  std::vector<double> x;
  long int NumPts;
  NumPts = lrint((b-a)/dx + 1);
  double val = a;
  while (val<=b){
    x.push_back( val );
    val += dx;
  }
  return x;
}

std::vector<double> UtilityFunctions::random(const int NumPts) {
  std::vector<double> x;
  for (int i = 0; i < NumPts; i++) {
    x.push_back(rand());
  }
  return x;
}

std::vector<double> UtilityFunctions::ones(const int NumPts){
  std::vector<double> x;
  for (int i=0; i<NumPts; i++){
    x.push_back( 1.0 );
  }
  return x;
}

std::vector<double> UtilityFunctions::zeros(const int NumPts){
  std::vector<double> x;
  for (int i=0; i<NumPts; i++){
    x.push_back( 0.0 );
  }
  return x;
}

std::vector<double> UtilityFunctions::multiply(const std::vector<double> &v, const double x){
  // takes a scalar value, x, and multiplies vector v by x, elementwise.
  int size = v.size();
  std::vector<double> target(size,0);
  for (int i=0; i<size; i++){
    target[i] = v[i]*x;
  }
  return target;
}

std::vector<double> UtilityFunctions::multiply(const std::vector<double> &v, const std::vector<double> &u){
  // takes two vectors, v & u, and does element wise multiplication
  if (v.size() != u.size()){
    printf("\n\n%lu, %lu\n",v.size(),u.size());
    throw std::invalid_argument("\n\nMultiply -> size of vectors u and v must be the same!!\n");
  }
  int size = v.size();
  std::vector<double> target(size,0);
  for (int i=0; i<size; i++){
    target[i] = v[i]*u[i];
  }
  return target;
}

std::vector<double> UtilityFunctions::divide(const std::vector<double> &v, const double x){
  // takes a scalar value, x, and divides vector v by x, elementwise.
  int size = v.size();
  std::vector<double> target(size,0);
  for (int i=0; i<size; i++){
    target[i] = v[i]/x;
    if (isfinite(target[i]) == false){
      printf("  -> UtilityFunctions::divide(const std::vector<double> &v, const double x)\n");
      printf("     division is producing a non-finite number! Quitting.\n");
      exit(-1);
    }
  }
  return target;
}

std::vector<double> UtilityFunctions::divide(const std::vector<double> &v, const std::vector<double> &u){
  // takes two vectors, v & u, and does element wise subtraction
  if (v.size() != u.size()){
    printf("\n\n%lu, %lu\n",v.size(),u.size());
    throw std::invalid_argument("\n\nDivide -> size of vectors u and v must be the same!!\n");
  }
  int size = v.size();
  std::vector<double> target(size,0);
  for (int i=0; i<size; i++){
    if (u[i] == 0.0){
      printf("\n\n dividing by 0.0! This may lead to unexepected results. Check it out\n");
      printf("%.4e/%.4e = %.4e\n", v[i],u[i],v[i]/u[i]);
      exit(-1);
    }
    target[i] = v[i]/u[i];
  }
  return target;
}

std::vector<double> UtilityFunctions::subtract(const std::vector<double> &v, const std::vector<double> &u){
  // takes two vectors, v & u, and does element wise subtraction
  if (v.size() != u.size()){
    printf("\n\n%lu, %lu\n",v.size(),u.size());
    throw std::invalid_argument("\n\nSubtract -> size of vectors u and v must be the same!!\n");
  }
  int size = v.size();
  std::vector<double> target(size,0);
  for (int i=0; i<size; i++){
    target[i] = v[i]-u[i];
  }
  return target;
}

std::vector<double> UtilityFunctions::add(const std::vector<double> &v, const std::vector<double> &u){
  // takes two vectors, v & u, and does element wise addition
  if (v.size() != u.size()){
    printf("\n\n%lu, %lu\n",v.size(),u.size());
    throw std::invalid_argument("\nadd -> size of vectors u and v must be the same!!\n");
  }
  int size = v.size();
  std::vector<double> target(size,0);
  for (int i=0; i<size; i++){
    target[i] = v[i]+u[i];
  }
  return target;
}

std::vector<double> UtilityFunctions::add(const std::vector<double> &v, const double x){
  // adds x to each element of v
  int size = v.size();
  std::vector<double> target(size,0);
  for (int i=0; i<size; i++){
    target[i] = v[i] + x;
  }
  return target;
}


std::vector<double> UtilityFunctions::square(const std::vector<double> &v){
  // takes a vector, v, and squares each element
  int size = v.size();
  std::vector<double> target(size,0);
  for (int i=0; i<size; i++){
    target[i] = pow(v[i],2.0);
  }
  return target;
}

std::vector<double> UtilityFunctions::absolute(const std::vector<double> &v){
  // takes a vector, v, and takes the absolute value of each element
  int size = v.size();
  std::vector<double> target(size,0);
  for (int i=0; i<size; i++){
    target[i] = abs(v[i]);
  }
  return target;
}

double UtilityFunctions::sum(const std::vector<double> &v) {
  // takes a vector, v, and takes the absolute value of each element
  double target {0.0};
  for (double val : v)
    target += val;
  return target;
}

int UtilityFunctions::sum(const std::vector<int> &v) {
  // takes a vector, v, and takes the absolute value of each element
  int target{0};
  for (int val : v)
    target += val;
  return target;
}

//==========================================================================//
//                        General FEM operations                            //
//==========================================================================//
void UtilityFunctions::FEM1D_Eval(const FEMGrid &grid, ShapeFunction1D &shape,
        const double dx, const int el, const std::vector<double> &v, EvalFEM1D &mynums){
  // evaluate FEM functions in v
  double uhval {0.0}; double duhval {0.0}; int mynum;
  for (int i = 0; i < grid.order[el]; i++){
    mynum = grid.nod[el][i];
    uhval += v[mynum]*shape.psi[i];
    duhval += v[mynum]*shape.dpsi[i]/dx;
  }
  mynums.b = uhval;
  mynums.bp = duhval;
}

void UtilityFunctions::FEM1D_Eval(const FEMGrid &grid, ShapeFunction1D &shape,
        const double dx, const int el, const std::vector<std::vector<double>> &v, EvalFEM1D &mynums){
  // evaluate FEM functions in v
  int groups {static_cast<int>(v.size())};
  double uhval, duhval; int mynum;
  for (int g=0; g<groups; g++){
    uhval = 0.0; duhval = 0.0;
    for (int i = 0; i < grid.order[el]; i++){
      mynum = grid.nod[el][i];
      uhval += v[g][mynum]*shape.psi[i];
      duhval += v[g][mynum]*shape.dpsi[i]/dx;
    }
    mynums.u[g] = uhval;
    mynums.up[g] = duhval;
  }
}

void UtilityFunctions::FEM1D_Eval(const FEMGrid &grid, ShapeFunction1D &shape,
        const double dx, const int el, const std::vector<std::vector<std::vector<double>>> &v, EvalFEM1D &mynums){
  // evaluate FEM functions in v
  int groups {static_cast<int>(v.size())};
  int N;
  double uhval, duhval; int mynum;
  for (int g=0; g<groups; g++){
    N = static_cast<int>(v[g].size());
    for (int ide=0; ide<N; ide++){
      uhval = 0.0; duhval = 0.0;
      for (int i = 0; i < grid.order[el]; i++){
        mynum = grid.nod[el][i];
        uhval += v[g][ide][mynum]*shape.psi[i];
        duhval += v[g][ide][mynum]*shape.dpsi[i]/dx;
      }
      mynums.u_i[g][ide] = uhval;
      mynums.up_i[g][ide] = duhval;
    }
  }
}

double UtilityFunctions::SqL2Norm(const std::vector<double> &v, const FEMGrid &xgrid){
  QuadParams1D qps; ShapeFunction1D shape; EvalFEM1D evald;
  get1D_QPs(xgrid.maxord, qps);
  double xL; double xR; double dx;
  double tmp = 0.0;
  for (int el = 0; el < xgrid.nels; el++){
    xL = xgrid.xnod[ xgrid.nod[el][0] ];
    xR = xgrid.xnod[ xgrid.nod[el][xgrid.order[el]-1] ];
    dx = (xR-xL)/2.0;
    for (int l=0; l<qps.nw; l++){
      EvalBasis1D(qps.xw[l], xgrid.order[el], shape);
      FEM1D_Eval(xgrid,shape,dx,el,v,evald);
      // complete integration over element
      tmp += pow(evald.b,2.0)*qps.w[l]*dx;
    }
  }
  return tmp;
}

double UtilityFunctions::FD_PNorm(const std::vector<double> &v, const double p){
  double tmp = 0.0;
  for (double elem : v){
    tmp += pow(abs(elem),p);
  }
  tmp = pow(tmp,1/p);
  return tmp;
}

//==========================================================================//
//                         General PGD operations                           //
//==========================================================================//
double UtilityFunctions::FEM_Enr_Norm(const std::vector<double> &v, const FEMGrid &xgrid){
  QuadParams1D qps; ShapeFunction1D shape; EvalFEM1D evald;
  get1D_QPs(xgrid.maxord, qps);
  double xL; double xR; double dx;
  double tmp = 0.0;
  for (int el = 0; el < xgrid.nels; el++){
    xL = xgrid.xnod[ xgrid.nod[el][0] ];
    xR = xgrid.xnod[ xgrid.nod[el][xgrid.order[el]-1] ];
    dx = (xR-xL)/2.0;
    for (int l=0; l<qps.nw; l++){
      EvalBasis1D(qps.xw[l], xgrid.order[el], shape);
      FEM1D_Eval(xgrid,shape,dx,el,v,evald);
      // complete integration over element
      tmp += pow(evald.b,2.0)*qps.w[l]*dx;
    }
  }
  return tmp/(pow(xgrid.xnod.back() - xgrid.xnod.front(), 2.0));
}

double UtilityFunctions::FD_Enr_Norm(const std::vector<double> &v, const double a, const double b){
  double tmp = 0.0;
  for (double elem : v){//(int i=0; i<v.size(); i++){
    tmp += pow(elem,2.0);
  }
  return tmp/(pow(b-a,2.0));
}


//==========================================================================//
//          Quadrature information for numerical integration.               //
//==========================================================================//
void UtilityFunctions::get2D_QPs(int numnodes, QuadParams2D &qps){
  int nw;
  std::vector<std::vector<double>> xw;
  std::vector<double> w;

  std::vector<double> tmp;
  if (numnodes == 4){
    nw = 4;
    tmp.push_back(-1./sqrt(3.)); tmp.push_back(-1./sqrt(3.));
    xw.push_back(tmp); tmp.clear();
    tmp.push_back(1./sqrt(3.)); tmp.push_back(-1./sqrt(3.));
    xw.push_back(tmp); tmp.clear();
    tmp.push_back(1./sqrt(3.)); tmp.push_back(1./sqrt(3.));
    xw.push_back(tmp); tmp.clear();
    tmp.push_back(-1./sqrt(3.)); tmp.push_back(1./sqrt(3.));
    xw.push_back(tmp); tmp.clear();
    w.push_back(1.);
    w.push_back(1.);
    w.push_back(1.);
    w.push_back(1.);
  }
  else{
    printf("Incorrect quadrature order.\n");
    printf("quad order for num nodes = %d does not exist.\n",numnodes);
    exit(1);
  }
  qps.nw = nw;
  qps.xw = xw;
  qps.w = w;
}

void UtilityFunctions::get1D_QPs(int maxord, QuadParams1D &qps){
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
    printf("Incorrect quadrature order.\n");
    printf("Order = %d does not exist.\n",maxord);
    exit(1);
  }
  qps.nw = nw;
  qps.xw = xw;
  qps.w = w;
}

void UtilityFunctions::EvalBasis1D(double x, int n, ShapeFunction1D &shape){
  // shape function on reference element (-1,1)
  // FEM order (n), has (n+1) quadrature points, and integrates (2n-1) functions EXACTLY
  std::vector<double> y;
  std::vector<double> dy;
  if (n==2){ //for linear FE's (n=1). Integrates up to linear FE's exactly.
    y.push_back(0.5*(1.0-x));
    y.push_back(0.5*(1.0+x));
    dy.push_back(-0.5);
    dy.push_back(0.5);
  }
  else if (n==3){ //for quadratic FE's (n=2). Integrates up to cubic FE's exactly.
    y.push_back((pow(x,2)-x)/2.0);
    y.push_back(1-pow(x,2));
    y.push_back((pow(x,2)+x)/2.0);
    dy.push_back(x-1./2.);
    dy.push_back(-2.0*x);
    dy.push_back(x+1./2.);
  }
  else if (n==4){ //for cubic FE's (n=3). Integrates up to 5th order FE's exactly.
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
    printf("Order = %d shape function does not exist.",n);
    exit(1);
  }
  shape.psi = y;
  shape.dpsi = dy;
}

//==========================================================//
//          General Purpose PGD I/O Functions               //
//==========================================================//
void UtilityFunctions::make_dir(const char *path){
  int dir_err;
  dir_err = mkdir(path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if (dir_err == -1 && errno != 17){
    printf("\n\nError in making time directory \"%s\"\n  %s\n", path,strerror(errno));
    exit(-1);
  }
}
