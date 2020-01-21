// INCLUDE GUARD
#ifndef __Utilities_H_INCLUDED__
#define __Utilities_H_INCLUDED__

// headers for dependents (what this class needs to function)
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sys/stat.h> // used for mkdir in Write_PGDSoln
#include <errno.h>

struct QuadParams1D{
  int nw;
  std::vector<double> xw;
  std::vector<double> w;
};

struct QuadParams2D{
  int nw;
  std::vector<std::vector<double>> xw;
  std::vector<double> w;
};

struct ShapeFunction1D{
  std::vector<double> psi;
  std::vector<double> dpsi;
};

struct FEMGrid{
  // Space mesh attributes
  int nels;
  int maxord;
  std::vector<int> order;
  std::vector< std::vector< int > > nod;
  std::vector< double > xnod;
  int nnodes;
  std::vector<double> bounds;
  double hel;
};

struct EvalFEM1D{
  /* NOTE, groups can be either total num of energy or precursor groups */
  EvalFEM1D(){};
  EvalFEM1D(const int groups){
    u.resize(groups);
    up.resize(groups);
  }
  EvalFEM1D(const int groups, const int NumEnr){
    u.resize(groups);
    up.resize(groups);
    u_i.resize(groups, std::vector<double>(NumEnr));
    up_i.resize(groups, std::vector<double>(NumEnr));
  }
  double b;
  double bp;
  /* data for group wise FEM functions */
  std::vector<double> u;
  std::vector<double> up;
  /* data for PGD specific functions for group wise data for stored enr steps */
  // vector of size of number of energy groups
  // each energy group is a vector w/ entries for each enr step.
  std::vector<std::vector<double>> u_i;
  std::vector<std::vector<double>> up_i;
};

struct EvalFEM2D{
  EvalFEM2D(){};
  EvalFEM2D(const int groups){
    u.resize(groups);
    up_x.resize(groups);
    up_y.resize(groups);
  }
  double b;
  double bp;
  /* data for group wise FEM functions */
  std::vector<double> u;
  std::vector<double> up_x;
  std::vector<double> up_y;
};

class UtilityFunctions{
  /*
  - class for general purpose computing functions
  - in order for a function to be used outside of Solvers (i.e. in PGD algorithm)
    it needs to be public.
  */
  public:
    UtilityFunctions(){}; //empty constructor
    // ----------------------------------------------------------------
    // General purpose numerical functions
    std::vector<double> linspace(const double a, const double b, const int NumPts = 100);
    std::vector<int> linspace(const int a, const int b);
    std::vector<double> logspace(const double a, const double b, const int NumPts = 100);
    std::vector<double> arange(const double a, const double b, const double dx);
    std::vector<double> random(const int NumPts = 100);
    std::vector<double> ones(const int NumPts = 100);
    std::vector<double> zeros(const int NumPts = 100);
    std::vector<double> multiply(const std::vector<double> &v, const double x);
    std::vector<double> multiply(const std::vector<double> &v, const std::vector<double> &u);
    std::vector<double> divide(const std::vector<double> &v, const double x);
    std::vector<double> divide(const std::vector<double> &v, const std::vector<double> &u);
    std::vector<double> subtract(const std::vector<double> &v, const std::vector<double> &u);
    std::vector<double> add(const std::vector<double> &v, const std::vector<double> &u);
    std::vector<double> add(const std::vector<double> &v, const double x);
    std::vector<double> square(const std::vector<double> &v);
    std::vector<double> absolute(const std::vector<double> &v);
    double sum(const std::vector<double> &v);
    int sum(const std::vector<int> &v);
    // ----------------------------------------------------------------
    /* General FEM functions */
    void FEM1D_Eval(const FEMGrid &grid, ShapeFunction1D &shape, const double dx, const int el, const std::vector<double> &v, EvalFEM1D &mynums);
    void FEM1D_Eval(const FEMGrid &grid, ShapeFunction1D &shape, const double dx, const int el, const std::vector<std::vector<double>> &v, EvalFEM1D &mynums);
    void FEM1D_Eval(const FEMGrid &grid, ShapeFunction1D &shape, const double dx, const int el, const std::vector<std::vector<std::vector<double>>> &v, EvalFEM1D &mynums);
    double SqL2Norm(const std::vector<double> &v, const FEMGrid &grid);
    double FD_PNorm(const std::vector<double> &v, const double p=2.0);
    // ----------------------------------------------------------------
    // General PGD operations
    double FEM_Enr_Norm(const std::vector<double> &v, const FEMGrid &xgrid);
    double FD_Enr_Norm(const std::vector<double> &v, const double a, const double b);
    // ----------------------------------------------------------------
    // Quad information for numerical integration
    void get1D_QPs(int maxord, QuadParams1D &qps);
    void get2D_QPs(int numnodes, QuadParams2D &qps);
    void EvalBasis1D(double x, int n, ShapeFunction1D &shape);
    // ----------------------------------------------------------------
    // General PGD i/o functions
    void make_dir(const char *path);
};

#endif
