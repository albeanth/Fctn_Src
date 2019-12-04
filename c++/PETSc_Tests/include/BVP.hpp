// INCLUDE GUARD
#ifndef __BVP_H_INCLUDED__
#define __BVP_H_INCLUDED__

// headers for dependents (what this class needs to function)
#include "SetUpGrid.hpp"    // <- parent class
#include "TestFunction.hpp" // <- parent class
#include <petscksp.h>
#include <vector>

struct PetscShapeFunction1D {
  Vec psi, dpsi;
};

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
    private:
        // PetscError code
        PetscErrorCode ierr;
        // Global PETSc matrices (n x n)
        Mat stiff, mass;
        Vec rhsf;
        // Local PETSc matrices (ord x ord)
        Mat m,k;
        Vec f;
        // matrices used for generation of local petsc matrices (m & k)
        Mat psi, dpsi, dpsiMat, psiMat;

        // Private member functions
        PetscErrorCode AssignEvaldBasis(const int ord, const double dx, const PetscShapeFunction1D &shape1d);
        PetscErrorCode PetscEvalBasis1D(const double x, const int n, PetscShapeFunction1D &shape);
};

#endif