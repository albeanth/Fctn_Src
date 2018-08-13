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
  double mu3;
};

class ODEClass_Test{
  public:
    ODEClass_Test(){
      mymus.mu0 = 0.35;
      mymus.mu1 = 0.15;
    }
    Params mymus;
    double tmp = 1.0;
    // double mu0 = 0.35;
    // double mu1 = 0.15;
    double mu3;
    Solution compute_ODE(std::vector<double> TBnds, int dim, double RelTol, double AbsTol);
};



// gsl_odeiv2_system sys
// {
//     &my_model::ode_gsl
// ,   nullptr
// ,   chemosc.NEQ
// ,   reinterpret_cast<void *>(::std::addressof(chemosc))
// };
