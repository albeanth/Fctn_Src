#include "Class_Test.hpp"

void adjust_mus(ODEClass_Test *test, double t){
  test->mymus.mu3 = test->mymus.mu0+test->mymus.mu1;
}

int ODE_fun(double t, const double y[], double f[], void *params){
  ODEClass_Test * test = static_cast<ODEClass_Test*>(params);
  adjust_mus(test, t);
  f[0] = -(test->mymus.mu3) * y[0];
  return GSL_SUCCESS;
}

void ODEClass_Test::compute_error(std::vector<double> &error){
  // compute relative error from analytic solution
  error.resize(soln.t.size());
  for (int i = 0; i<soln.t.size(); i++) {
    adjust_mus(this, soln.t[i]);
    analy_soln = 1.0 * exp(-this->mymus.mu3 * soln.t[i]);
    error[i] = (soln.y[i] - analy_soln)/analy_soln;
  }
}

void ODEClass_Test::compute_ODE(){
  printf("Computing ODE...\n");
  
  // clear out solutions in solution struct
  soln.t.clear();
  soln.y.clear();

  const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rkf45;
  gsl_odeiv2_system sys = {ODE_fun, NULL, static_cast<size_t>(dim), this};
  gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new (&sys, T, 1.0E-6, AbsTol, RelTol);

  double t {TBnds[0]};
  double t1 {TBnds[1]};
  double y[1] { 1.0 };

  soln.t.push_back(t);
  soln.y.push_back(y[0]);

  int status;
  double ti;
  for (int i = 1; i<=15; i++){
    ti = i * t1 / 15.0;
    status = gsl_odeiv2_driver_apply (d, &t, ti, y);
    if (status != GSL_SUCCESS) {
      printf("\nodeiv2 failed, gsl_errno = %s\n\n", gsl_strerror(status));
    }
    soln.t.push_back(t);
    soln.y.push_back(y[0]);
  }
  gsl_odeiv2_driver_free(d);
}

void ODEClass_Test::compute_Sparse_ODE(){
  /*
  Demonstrates the sparse linear algebra routine (GMRES) on a 1D Poisson eqauation over [0,1]
          \frac{\pd^2 u}{\pd x^2} = -\pi^2\sin(\pi x)
  with BC, u(0)=u(1)=0.0;
  
  The analytic solution is
          u(x) = sin(\pi x)

  The problem is discretized using a O(2) finite difference method (central differencing).
  */

  printf("Computing sparse ODE...\n");

  const size_t N = 100;                       /* number of grid points */
  const size_t n = N - 2;                     /* subtract 2 to exclude boundaries */
  const double h = 1.0 / (N - 1.0);           /* grid spacing */
  gsl_spmatrix *A = gsl_spmatrix_alloc(n, n); /* triplet format - used for creation */
  gsl_spmatrix *C;                            /* compressed format - used for linalg */
  gsl_vector *f = gsl_vector_alloc(n);        /* right hand side vector */
  gsl_vector *u = gsl_vector_alloc(n);        /* solution vector */
  size_t i;

  /* construct the sparse matrix for the finite difference equation */

  /* construct first row */
  gsl_spmatrix_set(A, 0, 0, -2.0);
  gsl_spmatrix_set(A, 0, 1, 1.0);

  /* construct rows [1:n-2] */
  for (i = 1; i < n - 1; ++i) {
    gsl_spmatrix_set(A, i, i + 1, 1.0);
    gsl_spmatrix_set(A, i, i, -2.0);
    gsl_spmatrix_set(A, i, i - 1, 1.0);
  }

  /* construct last row */
  gsl_spmatrix_set(A, n - 1, n - 1, -2.0);
  gsl_spmatrix_set(A, n - 1, n - 2, 1.0);

  /* scale by h^2 */
  gsl_spmatrix_scale(A, 1.0 / (h * h));

  /* construct right hand side vector */
  for (i = 0; i < n; ++i) {
    double xi = (i + 1) * h;
    double fi = -M_PI * M_PI * sin(M_PI * xi);
    gsl_vector_set(f, i, fi);
  }

  /* convert to compressed column format */
  C = gsl_spmatrix_ccs(A);

  /* now solve the system with the GMRES iterative solver */
  {
    const double tol = 1.0e-6;  /* solution relative tolerance */
    const size_t max_iter = 10; /* maximum iterations */
    const gsl_splinalg_itersolve_type *T = gsl_splinalg_itersolve_gmres;
    /* allocates workspace for gmres solver. 0 indicates default value for subspace dim */
    gsl_splinalg_itersolve *work = gsl_splinalg_itersolve_alloc(T, n, 0);
    size_t iter = 0;
    double residual;
    int status;

    /* initial guess u = 0 */
    gsl_vector_set_zero(u);

    /* solve the system A u = f */
    do {
      status = gsl_splinalg_itersolve_iterate(C, f, tol, u, work);

      /* print out residual norm ||A*u - f|| */
      residual = gsl_splinalg_itersolve_normr(work);
      fprintf(stderr, "iter %zu residual = %.12e\n", iter, residual);

      if (status == GSL_SUCCESS)
        fprintf(stderr, "Converged\n");
    } while (status == GSL_CONTINUE && ++iter < max_iter);

    /* output solution */
    for (i = 0; i < n; ++i) {
      double xi = (i + 1) * h;
      double u_exact = sin(M_PI * xi);
      double u_gsl = gsl_vector_get(u, i);

      printf("%f %.12e %.12e\n", xi, u_gsl, u_exact);
    }

    gsl_splinalg_itersolve_free(work);
  }

  gsl_spmatrix_free(A);
  gsl_spmatrix_free(C);
  gsl_vector_free(f);
  gsl_vector_free(u);

}
