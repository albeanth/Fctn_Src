#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <stdexcept>

struct Data{
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> y_p;
};

struct Test{
  int size;
  std::string method;
};

class InterpClass_Test{
  public:
    Data set_xy(const int size, std::vector<double> Bnds);
    gsl_spline* get_spline(Test mydata);//const int size, std::string method);
    Data eval_spline(gsl_spline * MySpline, Data ScatterData, const int size, int flag);
    double integrate_spline(gsl_spline * MySpline, Data ScatterData, const int size);
};
