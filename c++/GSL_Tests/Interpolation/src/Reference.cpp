#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

int main (void){
  int i;
  double xi, yi, x[10], y[10];

  printf ("#m=0,S=17\n");

  for (i = 0; i < 10; i++){
    x[i] = i + 0.5 * sin (i);
    y[i] = i + cos (i * i);
    printf ("%g %g\n", x[i], y[i]);
  }
  printf ("#m=1,S=0\n");
  {
    const gsl_interp_type * Type = gsl_interp_cspline;
    gsl_interp_accel *acc = gsl_interp_accel_alloc (); // used in the evaluation func, gsl_spline_eval
    gsl_spline *spline = gsl_spline_alloc (Type, 10); // gets cspline method and sets 10 points

    gsl_spline_init (spline, x, y, 10);

    for (xi = x[0]; xi < x[9]; xi += 0.5){
      yi = gsl_spline_eval (spline, xi, acc);
      printf ("%g %g\n", xi, yi);
    }
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
  }
  return 0;
}
