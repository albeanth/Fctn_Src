// INCLUDE GUARD
#ifndef __NL_H_INCLUDED__
#define __NL_H_INCLUDED__

// headers for dependents (what this class needs to function)
#include "PetscFEFuncs.hpp" // <- parent class
#include "SetUpGrid.hpp"    // <- parent class
#include "TestFunction.hpp" // <- parent class
#include <petscsnes.h>
#include <petscdraw.h>
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
    Vec loc_mass_src, glo_mass_src;      // mms source for cons. of mass
    Vec loc_momen_src, glo_momen_src;    // mms source for cons. of momentum
    Vec loc_efluid_src, glo_efluid_src;  // mms source for cons. of fluid energy
};

struct MonitorCTX{
  PetscViewer viewer1, viewer2, viewer3;
  PetscDraw draw;
};

// function declarations for nonlinear function and jacobian forms
extern PetscErrorCode FormFunction(SNES snes, Vec x, Vec f, void *ctx);
extern PetscErrorCode FormJacobian(SNES snes, Vec x, Mat jac, Mat B, void *ctx);
extern PetscErrorCode Monitor(SNES snes, PetscInt a, PetscReal b, void* ctx);

class NonLinear : public TestFunction, public SetUpGrid, public PetscFEFuncs{
    /*
    uses petsc snes library to solve nonlinear differential equations
    */
    public:
        NonLinear(int argc, char **args, const int test_num, const bool hetero, const std::vector<double> &Bnds)
        : TestFunction(test_num, hetero), SetUpGrid(Bnds), PetscFEFuncs()
        {
          /* ------ PETSc initialization information ------ */
          PetscInitialize(&argc, &args, (char *)0, help);
          MPI_Comm_size(PETSC_COMM_WORLD, &size);
        };
        // Public Member Variables
        Vec velocity; /* global solution */
        Vec density;  /* global solution */
        Vec energy;  /* global solution */
        PetscScalar l2Err_Vel, l2Err_Rho, l2Err_Em; /* l2 error */
        PetscScalar h1Err_Vel, h1Err_Rho, h1Err_Em; /* h1 error */
        ApplicationCTX ctx; /* Instance of ApplicationCTX struct */
        PetscInt N; /* number of nodes in stencil * number of governing equations */
        MonitorCTX monCTX;
        // Public Member Functions
        PetscErrorCode NL_1D();
    private:
        // Private Member Variables
        double xL, xR, dx, x; /* cell specific information */
        PetscScalar src_mass, src_momen, src_energy;
        PetscMPIInt size;
        PetscErrorCode ierr;    /* petsc error code */
        SNES snes;              /* nonlinear solver context */
        Vec soln, residual;     /* local solution and residual vectors */
        Mat J;                  /* Jacobian matrix */
        KSP ksp;                /* linear solver context */
        PC pc;                  /* preconditioner context */
        /* Petsc Vectors for evaluated Shape Functions */
        Vec mass_basis_src;
        Vec momen_basis_src;
        Vec efluid_basis_src;
        /* Initialize Quadrature Parameters */
        QuadParams1D qps1d;
        // Private Member Functions
        PetscErrorCode NLSolve();
        PetscErrorCode InitializeLocalRHSF();
        PetscErrorCode VelRho_L2Error();
        PetscErrorCode EvalBasis(const double x, const int ord);
        PetscErrorCode Initialize_NL_1D();
        PetscErrorCode Local2Global(const int el);
};

#endif