#include "Interp_Test.h"
#include <iostream>
#include <vector>

int main(){
  const int size = 10;
  std::vector<double> Bnds;
  Bnds.push_back(0.0);
  Bnds.push_back(M_PI);
  const gsl_interp_type * Type = gsl_interp_cspline;

  InterpClass_Test test;// get object instance of InterpClass_Test
  Data ScatterData; // for reference data
  Data Approx; // for spline approximated data
  double Integral;

  ScatterData = test.set_xy(size, Bnds);
  gsl_spline * MySpline = test.get_spline(Type, size);
  Approx = test.eval_spline(MySpline, ScatterData, size, 1);
  Integral = test.integrate_spline(MySpline, ScatterData, size);
  printf("# Spline Approx.\n");
  for (int i=0; i<Approx.x.size(); i++){
    printf("%.5e % .5e % .5e\n",Approx.x[i],Approx.y[i],Approx.y_p[i]);
  }
  printf("\n%.8e\n",Integral);
}
