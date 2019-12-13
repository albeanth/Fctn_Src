// INCLUDE GUARD
#ifndef __TestFunction_H_INCLUDED__
#define __TestFunction_H_INCLUDED__

// headers for dependents (what this class needs to function)
#include "UtilityFunctions.hpp" // <- parent class
#include <iostream>
#include <stdexcept>

class TestFunction{
    /*
    Defines the function, derivative, and source required for
    a method of manufactured solutions (MMS) problem set up.
    Problem #'s < 100:
        a simple one-group neutron diffusion equation type problem
        -(d/dx * D(x) du/dx) + Sigma_a(x) u(x)= f(x)
    Problem #'s \geq 100
        a compressible mass/momentum equations
        \frac{d}{dx}(\rho u) = f_m(x)
        \frac{d}{dx}(\rho u^2) = f_{\rho} - \frac{dP_m}{dx}
    */
    public:
        TestFunction(const int problem, const bool material) :
        /*
            problem -> problem selection flag
            material -> homogeneous/heterogeneous/ringing
        */
        selection{problem}, hetero{material}{ // initialize member variables
            display_selection();              // call display selection function
        }
        public:
            /* Public member variables */
            const int selection;
            const bool hetero;
            char help[120];
            /* Public member functions */
            void display_selection();
            // neutron diffusion or fluid velocity
            double u(const double x);
            double up(const double x);
            double upp(const double x);
            // fluid density and pressure
            double rho(const double x);
            double rhop(const double x);
            double pressure(const double x);
            double pressurep(const double x);
            // cross sections for neutron diffusion
            double SigA(const double x);
            double D(const double x);
            double Dx(const double x);
            // MMS sources
            double MMS_Src(const double x);          // neutron diffusion
            double MMS_Src_Mass(const double x);     // cons. of mass
            double MMS_Src_Momentum(const double x); // cons. of momentum
};
#endif