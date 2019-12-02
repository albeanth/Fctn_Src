// INCLUDE GUARD
#ifndef __BVP_H_INCLUDED__
#define __BVP_H_INCLUDED__

// headers for dependents (what this class needs to function)
#include "SetUpGrid.hpp"    // <- parent class
#include "TestFunction.hpp" // <- parent class
#include <petscksp.h>
#include <vector>

class BVP : public TestFunction, public SetUpGrid{
    /*
    solves a two-point boundary value problem (BVP) using continuous
    or discontinuous Galerkin finite elements
    */
    public:
    BVP(const int test_num, const bool hetero): TestFunction(test_num, hetero), SetUpGrid(){};
    // Public member variables
    std::vector<double> soln;
    // Public member functions
    PetscErrorCode CFEM_1D(int argc, char **args);
};

#endif