#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "GaussianIntegration.h"

#ifndef _OPENMP
fprintf( stderr, "OpenMP is not supported here -- sorry.\n" );
return 1;
#endif
#include <omp.h>

using namespace std;

int main()
{

  int NTHREADS = 4;
  double time0;
  double time1;
  GaussianIntegration mytest;
  double test_integral_serial; double test_integral_parallel;

  // test_integral_serial = mytest.GaussInt_1D(0,M_PI,100);

  //test_integral_serial = mytest.GaussInt_2D_Serial(0,M_PI,100, 0,1,100);;
  //test_integral_parallel = mytest.GaussInt_2D_Parallel(0,M_PI,100, 0,1,100);
  
  time0 = omp_get_wtime();
  test_integral_serial = mytest.GaussInt_3D_Serial(0,M_PI,100, 0,M_PI,100, 0,M_PI,100);
  time1 = omp_get_wtime();
  printf("Serial Calculation took %.4lf \n",time1-time0);
  time0 = omp_get_wtime();
  test_integral_parallel = mytest.GaussInt_3D_Parallel(NTHREADS, 0,M_PI,100, 0,M_PI,100, 0,M_PI,100);
  time1 = omp_get_wtime();
  printf("Parallel Calculation took %.4lf \n",time1-time0);
  
  printf("\n\n");
  printf("Serial   -> %.12f\n",test_integral_serial);
  printf("Parallel -> %.12f\n",test_integral_parallel);
}
