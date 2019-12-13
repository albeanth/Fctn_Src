#include "TestFunction.hpp"

void TestFunction::display_selection(){
    /*
    print selection to screen
    */
    try {
      if (selection == 0) {
        strcpy(help, "\nLinear Solution test, Dirichlet BCs left and right, ");
      } 
      else if (selection == 1) {
        strcpy(help, "\nQuadratic Solution test, Dirichlet BCs left and right, ");
      }
      else if (selection == 100){
          strcpy(help, "\nNonlinear Solution test, Dirichlet BCs left and right, ");
      }
      else {
        throw std::invalid_argument("unknown problem selection");
      }
    }
    catch (const std::invalid_argument &e) {
        std::cerr << "Invalid argument: " << e.what() << std::endl;
        std::terminate();
    }

    if (hetero == false) {
        strcat(help, "with homogeneous cross sections.\n\n");
    } 
    else if (hetero == true) {
        strcat(help, "with heterogeneous cross sections.\n\n");
    }
}

double TestFunction::u(const double x){
    /*
    base function definition
    */
    double val{NAN};
    if (selection == 0){
        val = 1.0 - x;
    }
    else if (selection == 1){
        val = pow(x,2) - x + 1.0;
    }
    else if (selection == 100){
        val = 2.0 - x;
    }
    return val;
}

double TestFunction::up(const double x){
    /*
    derivative of base function
    */
   double val{NAN};
   if (selection == 0){
       val = -1.0;
   }
   else if (selection == 1){
       val = 2*x - 1.0;
   }
   else if (selection == 100){
       val = -1.0;
   }
   return val;
}

double TestFunction::upp(const double x){
    /*
    second derivative of base function
    */
   double val{NAN};
   if (selection == 0){
       val = 0.0;
   }
   else if (selection == 1){
       val = 2.0;
   }
   return val;
}

double TestFunction::rho(const double x){
  /*
  function definition for fluid density
  */
  double val{NAN};
  if (selection == 100){
    val = 1.0 + x;
  }
  return val;
}

double TestFunction::rhop(const double x) {
  /*
  first derivative of fluid density
  */
  double val{NAN};
  if (selection == 100) {
    val = 1.0;
  }
  return val;
}

double TestFunction::pressure(const double x){
  /*
  function definition for fluid pressure
  */
  double val{NAN};
  if (selection == 100){
    val = 1.0 + x;
  }
  return val;
}

double TestFunction::pressurep(const double x){
  /*
  first derivatice of fluid pressure
  */
  double val{NAN};
  if (selection == 100){
    val = 1.0;
  }
  return val;
}

double TestFunction::SigA(const double x){
    /*
    absorption cross section
    */
   double val {NAN};
   if (hetero == false){
       val = 1.0;
   }
   else if (hetero == true){
       val = x + 1.0;
   }
   return val;
}

double TestFunction::D(const double x){
    /*
    diffusion coefficient
    */
    return 1.0 / (3.0 * SigA(x));
}

double TestFunction::Dx(const double x){
    /*
    first derivative of diffusion coefficient
    */
    double val {NAN};
    if (hetero == false){
        val = 0.0;
    }
    else{
        val = -1.0 / (3.0*pow(x + 1.0,2.0));
    }
    return val;
}

double TestFunction::MMS_Src(const double x) {
  /*
   * source that defines the MMS problem for 
   * one-group neutron diffusion
   */
  return -(Dx(x) * up(x) + D(x) * upp(x)) + SigA(x) * u(x);
}

double TestFunction::MMS_Src_Mass(const double x){
    /*
     * source to define the MMS problem for 
     * the conservation of mass
     */
    return rho(x) * up(x) + rhop(x) * u(x);
}

double TestFunction::MMS_Src_Momentum(const double x){
    /*
     * source to define the MMS problem for 
     * the conservation of momentum
     */
    return rho(x) * 2.0 * u(x) * up(x) + rhop(x) * pow(u(x), 2.0);//- pressurep(x);
}