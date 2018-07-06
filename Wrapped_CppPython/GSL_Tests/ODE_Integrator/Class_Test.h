#include <stdio.h>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_odeiv2.h>

struct Solution{
  std::vector<double> t;
  std::vector<double> y;
};

class ODEClass_Test{
  public:
    Solution compute_ODE(std::vector<double> TBnds, int dim, double RelTolm, double AbsTol, double mu);
};
