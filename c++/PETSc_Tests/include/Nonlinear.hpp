// INCLUDE GUARD
#ifndef __NL_H_INCLUDED__
#define __NL_H_INCLUDED__

// headers for dependents (what this class needs to function)
#include "SetUpGrid.hpp"    // <- parent class
#include <petscsnes.h>
#include <vector>

// ANSI Color Codes for color output
#define GREEN "\x1b[32m"
#define RED "\x1b[31m"
#define YELLOW "\x1b[33m"
#define GRAY "\x1b[37m"
#define RESET "\x1b[0m"

class NonLinear : public SetUpGrid{
    /*
    uses petsc snes library to solve nonlinear differential equations
    */
    public:
        NonLinear(){
            strcpy(help, "\nNonlinear Solution test, Dirichlet BCs left and right, ");
        };
        // Public Member Variables
        // Public Member Functions
        PetscErrorCode NL_1D(int argc, char **args);
    private:
        // Private Member Variables
        PetscErrorCode ierr;
        char help[120];
        // Private Member Functions
};

#endif