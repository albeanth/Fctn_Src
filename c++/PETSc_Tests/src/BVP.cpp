#include "BVP.hpp"

PetscErrorCode BVP::CFEM_1D(int argc, char **args){
    /*
    1D CFEM solver uing PETSC ksp functionality
    */
    // PETSc initialization information
    PetscErrorCode ierr;
    PetscInt i, n = 10, col[3], its;
    PetscMPIInt size;
    PetscScalar value[3];
    ierr = PetscInitialize(&argc, &args, (char *)0, help); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
    if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only!");
    ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);

    /* PETSc matrices and vector declarations */
    PetscInt N = info.nnodes;
    PetscInt ord = info.order[0];
    // Initialize global stiffness and mass matrices and global rhsf vector
    Mat stiff, mass;
    ierr = MatCreate(PETSC_COMM_WORLD, &stiff); CHKERRQ(ierr);
    ierr = MatSetSizes(stiff, PETSC_DECIDE, PETSC_DECIDE, N, N); CHKERRQ(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD, &mass); CHKERRQ(ierr);
    ierr = MatSetSizes(mass, PETSC_DECIDE, PETSC_DECIDE, N, N); CHKERRQ(ierr);
    Vec rhsf;
    ierr = VecCreate(PETSC_COMM_WORLD, &rhsf); CHKERRQ(ierr);
    ierr = VecSetSizes(rhsf, PETSC_DECIDE, N); CHKERRQ(ierr);
    // Declare m, k, f local matrices/vectors
    Mat m,k;
    ierr = MatCreate(PETSC_COMM_WORLD, &m); CHKERRQ(ierr);
    ierr = MatSetSizes(m, PETSC_DECIDE, PETSC_DECIDE, ord, ord); CHKERRQ(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD, &k); CHKERRQ(ierr);
    ierr = MatSetSizes(k, PETSC_DECIDE, PETSC_DECIDE, ord, ord); CHKERRQ(ierr);
    Vec f;
    ierr = VecCreate(PETSC_COMM_WORLD, &f); CHKERRQ(ierr);
    ierr = VecSetSizes(f, PETSC_DECIDE, ord); CHKERRQ(ierr);
    // Declare shape function matrices
    Mat dpsiMat, psiMat;
    ierr = MatCreate(PETSC_COMM_WORLD, &dpsiMat); CHKERRQ(ierr);
    ierr = MatSetSizes(dpsiMat, PETSC_DECIDE, PETSC_DECIDE, ord, ord); CHKERRQ(ierr);

    QuadParams1D qps1d;
    get1D_QPs(info.maxord, qps1d);
    
    ShapeFunction1D shape1d;
    
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
            EvalBasis1D(qps1d.xw[l1], info.order[elem], shape1d);
            /* assign eval'd shape funcs to petsc matrices */

            /* evaluate known functions */
            Dval = D(x);
            SigAbs = SigA(x);
            fval = MMS_Src(x);

    // all petsc based functions need to end with PetscFinalize()
    ierr = PetscFinalize();
    return ierr;
}