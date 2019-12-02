// INCLUDE GUARD
#ifndef __TestFunction_H_INCLUDED__
#define __TestFunction_H_INCLUDED__

// headers for dependents (what this class needs to function)
#include "UtilityFunctions.hpp" // <- parent class
#include <iostream>
#include <stdexcept>

class TestFunction : public UtilityFunctions{
    /*
    Defines the function, derivative, and source required for
    a method of manufactured solutions (MMS) problem set up to
    a simple one-group neutron diffusion equation type problem
    - (d/dx * D(x) du/dx) + Sigma_a(x) u(x)= f(x)
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
        // Public member variables and functions
        const int selection;
        const bool hetero;
        void display_selection();
        double u(const double x);
        double up(const double x);
        double upp(const double x);
        double SigA(const double x);
        double D(const double x);
        double Dx(const double x);
        double f(const double x);
};
#endif