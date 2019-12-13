// INCLUDE GUARD
#ifndef __PetscFEFuncs_H_INCLUDED__
#define __PetscFEFuncs_H_INCLUDED__

// headers for dependents (what this class needs to function)
#include <petscksp.h>

// Structure declarations
struct PetscShapeFunction1D {
  Vec psi, dpsi;
};

class PetscFEFuncs{
    /*
     *  contains common FE functions for use with Petsc
     */
    public:
        PetscFEFuncs(){};
        // public member functions
        PetscErrorCode PetscEvalBasis1D(const double x, const int n, PetscShapeFunction1D &shape);
    private:
        // private member variables
        PetscErrorCode ierr;
};

#endif