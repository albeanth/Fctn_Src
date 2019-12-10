#include "Nonlinear.hpp"

PetscErrorCode NonLinear::NL_1D(int argc, char **args){
    /*
    1D nonlinear solver using PETSc snes functionality
    */
    // PETSc initialization information
    PetscMPIInt size;
    ierr = PetscInitialize(&argc, &args, (char *)0, help); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
    if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only!");

    printf("in NonLinear.cpp\n");
    return ierr;
}