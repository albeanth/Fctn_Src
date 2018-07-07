#include <stdio.h>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_odeiv2.h>

struct Solution{
  std::vector<double> t;
  std::vector<double> y;
};

struct Params{
  double mu0;
  double mu1;
};

class ODEClass_Test{
  public:
    Solution compute_ODE(std::vector<double> TBnds, int dim, double RelTolm, double AbsTol, Params mus);//double mu);
};
