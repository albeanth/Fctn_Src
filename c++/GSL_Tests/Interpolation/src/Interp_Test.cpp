#include "Interp_Test.h"

Data InterpClass_Test::set_xy(const int size, std::vector<double> Bnds){
  Data mydata; // get struc instance of Data called mydata
  double dx = (Bnds[1]-Bnds[0])/(size-1);
  for (int i=0; i<size; i++){
    mydata.x.push_back( dx*i ); //i + 0.5 * sin (i) );//
    mydata.y.push_back( sin(mydata.x[i]) ); //i + cos (i * i) );//
  }
  return mydata;
}

gsl_spline* InterpClass_Test::get_spline(const gsl_interp_type * Type, const int size){
  return gsl_spline_alloc (Type, size);
}

Data InterpClass_Test::eval_spline(gsl_spline * MySpline, Data ScatterData, const int size, int flag){
  double yi; double ypi;
  Data mydata;
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline_init(MySpline, &ScatterData.x[0], &ScatterData.y[0], size);
  for (double xi = ScatterData.x[0]; xi <= ScatterData.x[size-1]; xi+=0.1){
    yi = gsl_spline_eval(MySpline, xi, acc);
    mydata.x.push_back(xi);
    mydata.y.push_back(yi);
    if (flag==1){
      ypi = gsl_spline_eval_deriv(MySpline, xi, acc);
      mydata.y_p.push_back(ypi);
    }
  }
  return mydata;
}

double InterpClass_Test::integrate_spline(gsl_spline * MySpline, Data ScatterData, const int size){
  double integral;
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  printf("\n%.5e, %.5e\n",ScatterData.x[0], ScatterData.x[size-1]);
  integral = gsl_spline_eval_integ(MySpline, ScatterData.x[0], ScatterData.x[size-1], acc);
  gsl_spline_free (MySpline);
  gsl_interp_accel_free (acc);
  return integral;
}
