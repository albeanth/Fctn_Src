#include "Interp_Test.hpp"
#include <iostream>

int main(){
  const int size = 10;
  std::vector<double> Bnds;
  Bnds.push_back(0.0);
  Bnds.push_back(M_PI);


  InterpClass_Test test;// get object instance of InterpClass_Test
  Data ScatterData; // for reference data
  Data Approx; // for spline approximated data
  double Integral;
  std::string flag = "cspline";

  ScatterData = test.set_xy(size, Bnds);
  gsl_spline * MySpline = test.get_spline(size, flag);
  Approx = test.eval_spline(MySpline, ScatterData, size, 1);
  Integral = test.integrate_spline(MySpline, ScatterData, size);
  gsl_spline_free(MySpline);
  printf("# Spline Approx.\n");
  for (int i=0; i<Approx.x.size(); i++){
    printf("%.5e % .5e % .5e\n",Approx.x[i],Approx.y[i],Approx.y_p[i]);
  }
  printf("\n%.8e\n",Integral);
}
