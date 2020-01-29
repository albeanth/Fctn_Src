// INCLUDE GUARD
#ifndef __NL_H_INCLUDED__
#define __NL_H_INCLUDED__

// headers for dependents (what this class needs to function)
#include "PetscFEFuncs.hpp" // <- parent class
#include "SetUpGrid.hpp"    // <- parent class
#include "TestFunction.hpp" // <- parent class
#include <petscsnes.h>
#include <vector>

// ANSI Color Codes for color output
#define GREEN "\x1b[32m"
#define RED "\x1b[31m"
#define YELLOW "\x1b[33m"
#define GRAY "\x1b[37m"
#define RESET "\x1b[0m"

// strcut definitions
struct ApplicationCTX{
    PetscScalar gamma_s = 5.0/3.0 - 1.0;
    Vec mass_src;      /* mms source for cons. of mass */
    Vec momen_src;     /* mms source for cons. of momentum */
    Vec energy_src;    /* mms source for cons. of momentum */
    PetscScalar mass_upwind;   /* upwind source for cons. of mass */
    PetscScalar momen_upwind;  /* upwind sources for cons. of momentum */
    PetscScalar energy_upwind; /* upwind source for cons. of momentum */
};

// function declarations for nonlinear function and jacobian forms
extern PetscErrorCode FormFunction(SNES snes, Vec x, Vec f, void *ctx);
extern PetscErrorCode FormJacobian(SNES snes, Vec x, Mat jac, Mat B, void *ctx);

class NonLinear : public TestFunction, public SetUpGrid, public PetscFEFuncs{
    /*
    uses petsc snes library to solve nonlinear differential equations
    */
    public:
        NonLinear(const int test_num, const bool hetero, const std::vector<double> &Bnds)
        : TestFunction(test_num, hetero), SetUpGrid(Bnds), PetscFEFuncs(){};
        // Public Member Variables
        Vec velocity; /* global solution */
        Vec density;  /* global solution */
        Vec energy;  /* global solution */
        PetscScalar l2Err_Vel, l2Err_Rho, l2Err_Em; /* l2 error */
        PetscScalar h1Err_Vel, h1Err_Rho, h1Err_Em; /* h1 error */
        // Public Member Functions
        PetscErrorCode NL_1D(int argc, char **args);
    private:
        // Private Member Variables
        PetscErrorCode ierr;    /* petsc error code */
        SNES snes;              /* nonlinear solver context */
        Vec x, r;               /* local solution and residual vectors */
        Mat J;                  /* Jacobian matrix */
        ApplicationCTX ctx;     /* Instance of ApplicationCTX struct */
        PetscInt index;         /* Index number for where to pull upwind values from */
        KSP ksp;                /* linear solver context */
        PC pc;                  /* preconditioner context */
        // Private Member Functions
        PetscErrorCode NLSolve(const int elem, PetscScalar vL, PetscScalar vR, PetscScalar rL, PetscScalar rR, PetscScalar eL, PetscScalar eR);
        PetscErrorCode InitializeLocalRHSF();
        PetscErrorCode VelRho_L2Error();
};

#endif