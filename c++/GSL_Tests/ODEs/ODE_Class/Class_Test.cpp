#include "Class_Test.h"

void adjust_mus(ODEClass_Test *test, double t){
  test->mu3 = (test->mymus.mu0+test->mymus.mu1)/test->tmp;
}

int ODE_fun(double t, const double y[], double f[], void *params){
  ODEClass_Test * mymus = static_cast<ODEClass_Test*>(params);
  adjust_mus(mymus, t);
  f[0] = -(mymus->mu3)*y[0];
  return GSL_SUCCESS;
}

Solution ODEClass_Test::compute_ODE(std::vector<double> TBnds, int dim, double RelTol, double AbsTol){//, Params mu){// double mu){
  printf("Computing ODE...\n");

  const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rkf45;

  gsl_odeiv2_system sys = {ODE_fun, NULL, static_cast<size_t>(dim), this};
  gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, T,1e-6, RelTol, AbsTol);

  double t = TBnds[0];
  double t1 = TBnds[1];
  double y[1] = { 1.0 };

  int status;
  Solution results;
  printf("%.12e %.12e\n", t, 1.0);
  for (int i = 1; i<=15; i++){
    double ti = i * t1 / 15.0;
    int status = gsl_odeiv2_driver_apply (d, &t, ti, y);
    if (status != GSL_SUCCESS)
        break;
    results.t.push_back(t);
    results.y.push_back(y[0]);
    printf("%.12e %.12e\n", t, y[0]);
  }
  gsl_odeiv2_driver_free(d);
  return results;
}
