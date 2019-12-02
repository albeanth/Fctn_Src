// INCLUDE GUARD
#ifndef __BVP_H_INCLUDED__
#define __BVP_H_INCLUDED__

// headers for dependents (what this class needs to function)
#include "TestFunction.hpp" // <- parent class
#include "SetUpGrid.hpp" // <- parent class
#include <vector>

class BVP : public TestFunction, public SetUpGrid{
    /*
    solves a two-point boundary value problem (BVP) using continuous
    or discontinuous Galerkin finite elements
    */
    public:
    BVP(const int a, const bool b): TestFunction(a, b), SetUpGrid(){};
    // Public member variables
    std::vector<double> soln;
    // Public member functions
    void CFEM_1D();
};

#endif