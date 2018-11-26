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
    ODEClass_Test(const double atol, const double rtol, const std::vector<double> &tbnds) 
    : AbsTol {atol}, RelTol {rtol}, tmp {1.0}, TBnds {tbnds} {
      mymus.mu0 = 0.35;
      mymus.mu1 = 0.15;
    }
    std::vector<double> TBnds;
    double AbsTol, RelTol;
    Params mymus;
    double tmp;
    Solution compute_ODE(int dim);
    Solution compute_Sparse_ODE(int dim);
};
