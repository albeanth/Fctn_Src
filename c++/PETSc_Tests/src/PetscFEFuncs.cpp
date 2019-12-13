#include "PetscFEFuncs.hpp"

PetscErrorCode PetscFEFuncs::PetscEvalBasis1D(const double x, const int n, PetscShapeFunction1D &shape){
  // shape function on reference element (-1,1)
  // FEM order (n), has (n+1) quadrature points, and integrates (2n-1) functions EXACTLY
  ierr = VecCreate(PETSC_COMM_WORLD, &shape.psi); CHKERRQ(ierr);
  ierr = VecCreate(PETSC_COMM_WORLD, &shape.dpsi); CHKERRQ(ierr);
  ierr = VecSetSizes(shape.psi, PETSC_DECIDE, n); CHKERRQ(ierr);
  ierr = VecSetSizes(shape.dpsi, PETSC_DECIDE, n); CHKERRQ(ierr);
  ierr = VecSetFromOptions(shape.psi); CHKERRQ(ierr);
  ierr = VecSetFromOptions(shape.dpsi); CHKERRQ(ierr);
  if (n==2){ //for linear FE's (n=1). Integrates up to linear FE's exactly.
    PetscInt i0 = 0, i1 = 1;
    PetscScalar v1,v2,v3,v4;
    v1 = 0.5 * (1.0 - x);
    v2 = 0.5 * (1.0 + x);
    v3 = -0.5;
    v4 = 0.5;
    ierr = VecSetValues(shape.psi, 1, &i0, &v1, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecSetValues(shape.psi, 1, &i1, &v2, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecSetValues(shape.dpsi, 1, &i0, &v3, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecSetValues(shape.dpsi, 1, &i1, &v4, INSERT_VALUES); CHKERRQ(ierr);
  }
  else if (n==3){ //for quadratic FE's (n=2). Integrates up to cubic FE's exactly.
    PetscInt i0 = 0, i1 = 1, i2 = 2;
    PetscScalar v1, v2, v3, v4, v5, v6;
    v1 = (pow(x,2)-x)/2.0;
    v2 = 1-pow(x,2);
    v3 = (pow(x,2)+x)/2.0;
    v4 = x-1./2.;
    v5 = -2.0*x;
    v6 = x+1./2.;
    ierr = VecSetValues(shape.psi, 1, &i0, &v1, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecSetValues(shape.psi, 1, &i1, &v2, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecSetValues(shape.psi, 1, &i2, &v3, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecSetValues(shape.dpsi, 1, &i0, &v4, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecSetValues(shape.dpsi, 1, &i1, &v5, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecSetValues(shape.dpsi, 1, &i2, &v6, INSERT_VALUES); CHKERRQ(ierr);
  }
  else if (n==4){ //for cubic FE's (n=3). Integrates up to 5th order FE's exactly.
    PetscInt i0 = 0, i1 = 1, i2 = 2, i3 = 2;
    PetscScalar v1, v2, v3, v4, v5, v6, v7, v8;
    v1 = 1./16.*(-9.*pow(x,3) + 9.*pow(x,2) + x - 1.);
    v2 = 1./16.*(27.*pow(x,3) - 9.*pow(x,2) - 27.*x + 9.);
    v3 = 1./16.*(-27.*pow(x,3) - 9.*pow(x,2) + 27.*x + 9.);
    v4 = 1./16.*(9.*pow(x,3) + 9.*pow(x,2) - x - 1.);
    v5 =  1./16.*(-27.*pow(x,2) + 18.*x + 1.);
    v6 =  1./16.*(81.*pow(x,2) - 18.*x - 27.);
    v7 =  1./16.*(-81.*pow(x,2) - 18.*x + 27.);
    v8 =  1./16.*(27.*pow(x,2) + 18.*x - 1.);
    ierr = VecSetValues(shape.psi, 1, &i0, &v1, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecSetValues(shape.psi, 1, &i1, &v2, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecSetValues(shape.psi, 1, &i2, &v3, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecSetValues(shape.psi, 1, &i3, &v4, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecSetValues(shape.dpsi, 1, &i0, &v5, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecSetValues(shape.dpsi, 1, &i1, &v6, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecSetValues(shape.dpsi, 1, &i2, &v7, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecSetValues(shape.dpsi, 1, &i3, &v8, INSERT_VALUES); CHKERRQ(ierr);
  }
  else{
    printf("Order = %d shape function does not exist.",n);
    exit(1);
  }
  // Assemble the vectors
  ierr = VecAssemblyBegin(shape.psi); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(shape.dpsi); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(shape.psi); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(shape.dpsi); CHKERRQ(ierr);
  
  // view the vectors
  // printf(YELLOW "\nshape.psi" RESET "\n");
  // ierr = VecView(shape.psi, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr); // petsc_viewer_stdout_world allows for
  // printf(YELLOW "\nshape.dpsi" RESET "\n");
  // ierr = VecView(shape.dpsi, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr); // sequential and parallel vectors

  return ierr;
}