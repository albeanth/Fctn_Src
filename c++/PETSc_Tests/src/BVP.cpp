#include "BVP.hpp"

PetscErrorCode BVP::CFEM_1D(int argc, char **args){
    /*
    1D CFEM solver uing PETSC ksp functionality
    */
    // PETSc initialization information
    PetscInt n = 10;
    PetscMPIInt size;
    ierr = PetscInitialize(&argc, &args, (char *)0, help); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
    if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only!");
    ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);

    /* PETSc matrices and vector declarations */
    PetscInt N = info.nnodes;
    PetscInt ord = info.order[0];
    // Initialize global stiffness and mass matrices and global rhsf vector
    ierr = MatCreate(PETSC_COMM_WORLD, &stiff); CHKERRQ(ierr);
    ierr = MatSetSizes(stiff, PETSC_DECIDE, PETSC_DECIDE, N, N); CHKERRQ(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD, &mass); CHKERRQ(ierr);
    ierr = MatSetSizes(mass, PETSC_DECIDE, PETSC_DECIDE, N, N); CHKERRQ(ierr);
    ierr = VecCreate(PETSC_COMM_WORLD, &rhsf); CHKERRQ(ierr);
    ierr = VecSetSizes(rhsf, PETSC_DECIDE, N); CHKERRQ(ierr);
    // Declare m, k, f local matrices/vectors
    ierr = MatCreate(PETSC_COMM_WORLD, &m); CHKERRQ(ierr);
    ierr = MatSetSizes(m, PETSC_DECIDE, PETSC_DECIDE, ord, ord); CHKERRQ(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD, &k); CHKERRQ(ierr);
    ierr = MatSetSizes(k, PETSC_DECIDE, PETSC_DECIDE, ord, ord); CHKERRQ(ierr);
    ierr = VecCreate(PETSC_COMM_WORLD, &f); CHKERRQ(ierr);
    ierr = VecSetSizes(f, PETSC_DECIDE, ord); CHKERRQ(ierr);
    // Declare shape function matrices
    ierr = MatCreate(PETSC_COMM_WORLD, &dpsi); CHKERRQ(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD, &psi); CHKERRQ(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD, &dpsiMat); CHKERRQ(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD, &psiMat); CHKERRQ(ierr);
    ierr = MatSetSizes(dpsi, PETSC_DECIDE, PETSC_DECIDE, ord, ord); CHKERRQ(ierr);
    ierr = MatSetSizes(psi, PETSC_DECIDE, PETSC_DECIDE, ord, ord); CHKERRQ(ierr);
    ierr = MatSetSizes(dpsiMat, PETSC_DECIDE, PETSC_DECIDE, ord, ord); CHKERRQ(ierr);
    ierr = MatSetSizes(psiMat, PETSC_DECIDE, PETSC_DECIDE, ord, ord); CHKERRQ(ierr);

    QuadParams1D qps1d;
    get1D_QPs(info.maxord, qps1d);
    
    PetscShapeFunction1D shape1d;
    
    double xL, xR, dx, x;
    double Dval, SigAbs, fval;
    for (int elem = 0; elem < info.nels; elem ++){
        xL = info.xnod[ info.nod[elem][0] ];
        xR = info.xnod[ info.nod[elem][info.order[elem]-1] ];
        dx = (xR-xL)/2.0;
        
        for (int l1 = 0; l1 < qps1d.nw; l1++){
            /* map from ref elem to real elem */
            x = xL + (1.0 + qps1d.xw[l1]) * dx;
            /* evaluate basis functions */
            ierr = PetscEvalBasis1D(qps1d.xw[l1], info.order[elem], shape1d); CHKERRQ(ierr);
            /* assign eval'd shape funcs to petsc matrices */

            /* evaluate known functions */
            Dval = D(x);
            SigAbs = SigA(x);
            fval = MMS_Src(x);

    // all petsc based functions need to end with PetscFinalize()
    ierr = PetscFinalize();
    return ierr;
}

PetscErrorCode BVP::PetscEvalBasis1D(const double x, const int n, PetscShapeFunction1D &shape){
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
  VecAssemblyBegin(shape.psi);
  VecAssemblyBegin(shape.dpsi);
  VecAssemblyEnd(shape.psi);
  VecAssemblyEnd(shape.dpsi);
  
  // view the vectors
  VecView(shape.psi, PETSC_VIEWER_STDOUT_WORLD);  // petsc_viewer_stdout_world allows for
  VecView(shape.dpsi, PETSC_VIEWER_STDOUT_WORLD); // sequential and parallel vectors

  return ierr;
}