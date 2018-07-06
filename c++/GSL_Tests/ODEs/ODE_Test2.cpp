#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

int func (double t, const double y[], double f[], void *params){
  (void)(t); /* avoid unused parameter warning */
  double mu = *(double *)params;
  f[0] = -mu*y[0];
  return GSL_SUCCESS;
}

int main (void){
  const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rkf45;

  gsl_odeiv2_step * s = gsl_odeiv2_step_alloc (T, 1);
  gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (1e-6, 0.0);
  gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc (1);

  double mu = 0.5;
  gsl_odeiv2_system sys = {func, NULL, 1, &mu};

  double t = 0.0, t1 = 1.0;
  double h = 1e-6;
  double y[1] = { 1.0 };

  while (t < t1){
      int status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &t, t1, &h, y);
      if (status != GSL_SUCCESS)
          break;
      printf ("%.5e %.5e\n", t, y[0]);
  }
  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);
  return 0;
}
