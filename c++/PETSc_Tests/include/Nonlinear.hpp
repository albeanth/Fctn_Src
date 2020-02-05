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
    /* conservation of mass */
    Vec loc_mass_i, glo_mass_i;
    Vec loc_mass_ii, glo_mass_ii;
    Vec loc_mass_iii, glo_mass_iii;
    Vec loc_mass_iv, glo_mass_iv;
    Vec loc_mass_src, glo_mass_src;      // mms source for cons. of mass
    Vec mass_upwind;   // GLOBAL upwind source for cons. of mass
    /* conservation of momentum */
    Vec loc_momen_i, glo_momen_i;
    Vec loc_momen_ii, glo_momen_ii;
    Vec loc_momen_iii, glo_momen_iii;
    Vec loc_momen_iv, glo_momen_iv;
    Vec loc_momen_v, glo_momen_v;
    Vec loc_momen_vi, glo_momen_vi;
    Vec loc_momen_vii, glo_momen_vii;
    Vec loc_momen_viii, glo_momen_viii;
    Vec loc_momen_src, glo_momen_src;     // mms source for cons. of momentum
    Vec momen_upwind;  // GLOBAL upwind sources for cons. of momentum
    /* conservation of fluid energy */
    Vec loc_efluid_i, glo_efluid_i;
    Vec loc_efluid_ii, glo_efluid_ii;
    Vec loc_efluid_iii, glo_efluid_iii;
    Vec loc_efluid_iv, glo_efluid_iv;
    Vec loc_efluid_v, glo_efluid_v;
    Vec loc_efluid_vi, glo_efluid_vi;
    Vec loc_efluid_vii, glo_efluid_vii;
    Vec loc_efluid_viii, glo_efluid_viii;
    Vec loc_efluid_ix, glo_efluid_ix;
    Vec loc_efluid_x, glo_efluid_x;
    Vec loc_efluid_xi, glo_efluid_xi;
    Vec loc_efluid_xii, glo_efluid_xii;
    Vec loc_efluid_src, glo_efluid_src;    // mms source for cons. of momentum
    Vec efluid_upwind; // GLOBAL upwind source for cons. of momentum
};

struct MassBasis{
  Vec i, ii, iii, iv;
  Vec src;
};

struct MomentumBasis{
  Vec i, ii, iii, iv, v, vi ,vii, viii;
  Vec src;
};

struct EFluidBasis {
  Vec i, ii, iii, iv, v, vi, vii, viii, ix, x, xi, xii;
  Vec src;
};

// function declarations for nonlinear function and jacobian forms
extern PetscErrorCode FormFunction(SNES snes, Vec x, Vec f, void *ctx);
extern PetscErrorCode FormJacobian(SNES snes, Vec x, Mat jac, Mat B, void *ctx);

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
        // Public Member Functions
        PetscErrorCode NL_1D();
    private:
        // Private Member Variables
        double xL, xR, dx, x; /* cell specific information */
        PetscScalar src_mass, src_momen, src_energy;
        PetscMPIInt size;
        PetscInt N;
        PetscErrorCode ierr;    /* petsc error code */
        SNES snes;              /* nonlinear solver context */
        Vec soln, residual;     /* local solution and residual vectors */
        Mat J;                  /* Jacobian matrix */
        PetscInt index;         /* Index number for where to pull upwind values from */
        KSP ksp;                /* linear solver context */
        PC pc;                  /* preconditioner context */
        /* Petsc Vectors for evaluated Shape Functions */
        MassBasis mass_basis;
        MomentumBasis momen_basis;
        EFluidBasis efluid_basis;
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