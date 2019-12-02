#include "TestFunction.hpp"

void TestFunction::display_selection(){
    /*
    print selection to screen
    */
    try {
      if (selection == 0) {
        printf("\nLinear Solution test, Dirichlet BCs left and right, ");
      } 
      else if (selection == 1) {
        printf("\nNonlinear Solution test, Dirichlet BCs left and right, ");
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
        printf("with homogeneous cross sections.\n\n");
    } 
    else if (hetero == true) {
        printf("with heterogeneous cross sections.\n\n");
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

double TestFunction::f(const double x){
    /*
    source that defines the MMS problem
    */
    return -(Dx(x)*up(x) + D(x)*upp(x)) + SigA(x)*u(x);
}