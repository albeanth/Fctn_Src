// INCLUDE GUARD
#ifndef __BVP_H_INCLUDED__
#define __BVP_H_INCLUDED__

// headers for dependents (what this class needs to function)
#include "PetscFEFuncs.hpp" // <- parent class
#include "SetUpGrid.hpp"    // <- parent class
#include "TestFunction.hpp" // <- parent class
#include <petscksp.h>
#include <vector>

// ANSI Color Codes for color output
#define GREEN "\x1b[32m"
#define RED "\x1b[31m"
#define YELLOW "\x1b[33m"
#define GRAY "\x1b[37m"
#define RESET "\x1b[0m"

class BVP : public TestFunction, public SetUpGrid, public PetscFEFuncs{
    /*
    solves a two-point boundary value problem (BVP) using continuous
    or discontinuous Galerkin finite elements
    */
    public:
        BVP(const int test_num, const bool hetero, const std::vector<double> &Bnds)
        : TestFunction(test_num, hetero), SetUpGrid(Bnds), PetscFEFuncs()
        {
          petsc_one = 1.0;
        };
        // Public member variables
        std::vector<double> solution;
        // Public member functions
        PetscErrorCode CFEM_1D(int argc, char **args);
    private:
        // tmp dummy petsc variables.
        PetscScalar petsc_one;
        //
        PetscInt N, ord;
        PetscScalar *zero;
        PetscInt *i, *j, *ii, *jj;
        // PetscError code
        PetscErrorCode ierr;
        // Global PETSc matrices (n x n)
        Mat stiff, mass;
        Vec rhsf, rhsf_tmp;
        Vec soln;
        // Local PETSc matrices (ord x ord)
        Mat m,k;
        Vec f;
        // matrices used for generation of local petsc matrices (m & k)
        Mat psi, dpsi, dpsiMat, psiMat;
        // Linear system and preconditioner declarations
        Mat A;
        Vec b, xVec;
        KSP ksp;
        PC pc;
        // exact solution for error checking
        Vec ExactSoln;
        PetscScalar l2Err;
        PetscScalar h1Err;

        // Private member functions
        PetscErrorCode L2Error();
        PetscErrorCode DoLinearAlgebra();
        PetscErrorCode AssignLocalToGlobal(const std::vector<int> &tmp);
        PetscErrorCode InitializeLocalMatrices();
        PetscErrorCode AssignEvaldBasis(const double dx, const PetscShapeFunction1D &shape1d);
};

#endif