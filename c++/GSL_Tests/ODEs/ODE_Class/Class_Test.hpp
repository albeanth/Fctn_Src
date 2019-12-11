#include <stdio.h>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_splinalg.h>

struct Solution{
  std::vector<double> t;
  std::vector<double> y;
};

struct Params{
  double mu0;
  double mu1;
  double mu3;
};

class ODEClass_Test{
  public:
    ODEClass_Test(const double atol, const double rtol, const std::vector<double> &tbnds, const int d)
    : TBnds {tbnds}, AbsTol {atol}, RelTol {rtol}, dim {d} {
      mymus.mu0 = 0.35;
      mymus.mu1 = 0.15;
    }
    // public member variables
    Solution soln;             // struct for solution
    std::vector<double> TBnds; // bounds of time domain
    double AbsTol, RelTol;     // absolute and relative tolerances for IVP solver
    Params mymus;              // struct of params for ODE
    int dim;                   // dimension of ODE system
    double analy_soln;         // analytic solution for a given time, t
    // public member functions
    void compute_ODE();
    void compute_Sparse_ODE();
    void compute_error(std::vector<double> &error);
};
