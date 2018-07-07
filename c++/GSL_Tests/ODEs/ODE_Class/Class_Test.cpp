#include "Class_Test.h"

int ODE_fun(double t, const double y[], double f[], void *params){
  (void)(t); /* avoid unused parameter warning */
  Params * mymus = static_cast<Params*>(params);
  f[0] = -(mymus->mu0+mymus->mu1)*y[0];
  return GSL_SUCCESS;
}

Solution ODEClass_Test::compute_ODE(std::vector<double> TBnds, int dim, double RelTol, double AbsTol, Params mu){// double mu){
  // function to compute function in ODE_fun
  printf("Computing ODE...\n");

  gsl_vector * v = gsl_vector_alloc(2);
  gsl_vector_set(v,0,TBnds[0]);
  gsl_vector_set(v,1,TBnds[1]);
  const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rkf45;
  gsl_odeiv2_step * s = gsl_odeiv2_step_alloc (T, dim);
  gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (RelTol, AbsTol);
  gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc (dim);
  gsl_odeiv2_system sys = {ODE_fun, NULL, dim, &mu};

  double t = TBnds[0];
  double t1 = TBnds[1];
  double h = 1e-6;
  double y[1] = { 1.0 };

  Solution results;
  while (t < t1){
      int status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &t, t1, &h, y);
      if (status != GSL_SUCCESS)
          break;
      results.t.push_back(t);
      results.y.push_back(y[0]);
      // printf ("%.5e %.5e\n", t, y[0]);
  }
  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);
  return results;
}
