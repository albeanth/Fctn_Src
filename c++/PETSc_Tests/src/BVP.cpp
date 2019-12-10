#include "BVP.hpp"

PetscErrorCode BVP::CFEM_1D(int argc, char **args){
    /*
    1D CFEM solver uing PETSC ksp functionality
    */
    // PETSc initialization information
    PetscMPIInt size;
    ierr = PetscInitialize(&argc, &args, (char *)0, help); CHKERRQ(ierr);
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
    if (size != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only!");
    /* the following two lines are examples of how to set a PetscInt
     * and then overwrite it by passing something in at runtime
     */
    // PetscInt n = 10;
    // ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);

    /* PETSc matrices and vector declarations */
    N = info.nnodes;
    ord = info.order[0];
    zero = new PetscScalar[ord];
    i = new PetscInt[ord];
    j = new PetscInt[ord];
    for (int idx = 0; idx < ord; idx++) {
      zero[idx] = 0.0;
      i[idx] = idx;
      j[idx] = idx;
    }
    // Initialize global stiffness and mass matrices and global rhsf/_tmp and solution vectors
    ierr = MatCreate(PETSC_COMM_WORLD, &stiff); CHKERRQ(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD, &mass); CHKERRQ(ierr);
    ierr = MatSetSizes(mass, PETSC_DECIDE, PETSC_DECIDE, N, N); CHKERRQ(ierr);
    ierr = MatSetSizes(stiff, PETSC_DECIDE, PETSC_DECIDE, N, N); CHKERRQ(ierr);
    ierr = MatSetUp(stiff); CHKERRQ(ierr);
    ierr = MatSetUp(mass); CHKERRQ(ierr);
    ierr = VecCreate(PETSC_COMM_WORLD, &ExactSoln); CHKERRQ(ierr);
    ierr = VecCreate(PETSC_COMM_WORLD, &soln); CHKERRQ(ierr);
    ierr = VecCreate(PETSC_COMM_WORLD, &rhsf); CHKERRQ(ierr);
    ierr = VecCreate(PETSC_COMM_WORLD, &rhsf_tmp); CHKERRQ(ierr);
    ierr = VecSetSizes(ExactSoln, PETSC_DECIDE, N); CHKERRQ(ierr);
    ierr = VecSetSizes(soln, PETSC_DECIDE, N); CHKERRQ(ierr);
    ierr = VecSetSizes(rhsf, PETSC_DECIDE, N); CHKERRQ(ierr);
    ierr = VecSetSizes(rhsf_tmp, PETSC_DECIDE, N); CHKERRQ(ierr);
    ierr = VecSetFromOptions(ExactSoln); CHKERRQ(ierr);
    ierr = VecSetFromOptions(soln); CHKERRQ(ierr);
    ierr = VecSetFromOptions(rhsf); CHKERRQ(ierr);
    ierr = VecSetFromOptions(rhsf_tmp); CHKERRQ(ierr);
    ierr = VecZeroEntries(soln); CHKERRQ(ierr);
    // Declare m, k, f local matrices/vectors
    ierr = MatCreate(PETSC_COMM_WORLD, &m); CHKERRQ(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD, &k); CHKERRQ(ierr);
    ierr = MatSetSizes(m, PETSC_DECIDE, PETSC_DECIDE, ord, ord); CHKERRQ(ierr);
    ierr = MatSetSizes(k, PETSC_DECIDE, PETSC_DECIDE, ord, ord); CHKERRQ(ierr);
    ierr = MatSetUp(m); CHKERRQ(ierr);
    ierr = MatSetUp(k); CHKERRQ(ierr);
    ierr = VecCreate(PETSC_COMM_WORLD, &f); CHKERRQ(ierr);
    ierr = VecSetSizes(f, PETSC_DECIDE, ord); CHKERRQ(ierr);
    ierr = VecSetFromOptions(f); CHKERRQ(ierr);
    ierr = InitializeLocalMatrices(); CHKERRQ(ierr);
    // Declare shape function matrices
    ierr = MatCreate(PETSC_COMM_WORLD, &dpsi); CHKERRQ(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD, &psi); CHKERRQ(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD, &dpsiMat); CHKERRQ(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD, &psiMat); CHKERRQ(ierr);
    ierr = MatSetSizes(dpsi, PETSC_DECIDE, PETSC_DECIDE, ord, ord); CHKERRQ(ierr);
    ierr = MatSetSizes(psi, PETSC_DECIDE, PETSC_DECIDE, ord, ord); CHKERRQ(ierr);
    ierr = MatSetSizes(dpsiMat, PETSC_DECIDE, PETSC_DECIDE, ord, ord); CHKERRQ(ierr);
    ierr = MatSetSizes(psiMat, PETSC_DECIDE, PETSC_DECIDE, ord, ord); CHKERRQ(ierr);
    ierr = MatSetUp(dpsi); CHKERRQ(ierr);
    ierr = MatSetUp(psi); CHKERRQ(ierr);
    ierr = MatSetUp(dpsiMat); CHKERRQ(ierr);
    ierr = MatSetUp(psiMat); CHKERRQ(ierr);

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
            ierr = AssignEvaldBasis(dx, shape1d); CHKERRQ(ierr);
            /* evaluate known functions */
            Dval = D(x);
            SigAbs = SigA(x);
            fval = MMS_Src(x);
            /* scale and add local FE matrices */
            ierr = VecScale(shape1d.psi, fval * qps1d.w[l1] * dx); CHKERRQ(ierr);
            ierr = VecAXPY(f, petsc_one, shape1d.psi); CHKERRQ(ierr);

            ierr = MatScale(dpsiMat, Dval * qps1d.w[l1] * dx); CHKERRQ(ierr);
            ierr = MatAXPY(k, petsc_one, dpsiMat, SAME_NONZERO_PATTERN); CHKERRQ(ierr);

            ierr = MatScale(psiMat, SigAbs * qps1d.w[l1] * dx); CHKERRQ(ierr);
            ierr = MatAXPY(m, petsc_one, psiMat, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
        }
        /* print local FE matrices and RHS vector */
        // printf(YELLOW "\nf" RESET "\n");
        // ierr = VecView(f, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
        // printf(YELLOW "\nm" RESET "\n");
        // ierr = MatView(m, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
        // printf(YELLOW "\nk" RESET "\n");
        // ierr = MatView(k, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
        // exit(-1);
        
        /* map local matrices to global matrices */
        ierr = AssignLocalToGlobal(info.nod[elem]); CHKERRQ(ierr);

        // set all matrix/vector entries equal to zero (while maintaining structure)
        ierr = InitializeLocalMatrices(); CHKERRQ(ierr);
    }
    /* assemble global mass and stiffness matrices */
    ierr = MatAssemblyBegin(mass, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(mass, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyBegin(stiff, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(stiff, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    
    /* print out global mass matrix and stiffness matrices */
    // printf(YELLOW "\nGlobal Mass" RESET "\n");
    // ierr = MatView(mass, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    // printf(YELLOW "\nGlobal Stiff" RESET "\n");
    // ierr = MatView(stiff, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

    /* add global mass and stiffness matrices */
    ierr = MatAXPY(mass, petsc_one, stiff, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
    // printf(YELLOW "\nGlobal Total" RESET "\n");
    // ierr = MatView(mass, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    // exit(-1);

    /* assign the boundary conditions */
    if ((selection == 0) || (selection == 1)){ /* Dirichlet BCs */
      const PetscInt rows[] = {0, N - 1};
      const PetscScalar vals[] = { u(info.bounds[0]), u(info.bounds[1]) };
      ierr = VecSetValues(soln, ord, rows, vals, INSERT_VALUES);
      ierr = MatMult(mass, soln, rhsf_tmp); CHKERRQ(ierr);
      ierr = VecAXPY(rhsf, -petsc_one, rhsf_tmp); CHKERRQ(ierr);
      ierr = VecAssemblyBegin(rhsf); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(rhsf); CHKERRQ(ierr);

      /* initialize matrix A and rhsf b to pass to linear algebra solver */
      ierr = MatCreate(PETSC_COMM_WORLD, &A); CHKERRQ(ierr);
      ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, N-2, N-2); CHKERRQ(ierr);
      ierr = MatSetUp(A); CHKERRQ(ierr);
      ierr = VecCreate(PETSC_COMM_WORLD, &b); CHKERRQ(ierr);
      ierr = VecSetSizes(b, PETSC_DECIDE, N-2); CHKERRQ(ierr);
      ierr = VecSetFromOptions(b); CHKERRQ(ierr);
      
      /* extract relevant information from global matrix and insert into matrix A */
      PetscInt row_BCMat;                        // row number for matrix without BC row/col (i.e. BCMat)
      PetscInt col_GlobMat[ord], col_BCMat[ord]; // col numbers to extract from global/orig matrix and BC mat
      PetscScalar val[ord+1];                    // specific vals from global/orig matrix corresponding to col_GlobMat
      for (PetscInt idx = 1; idx<N-1; idx++){
        row_BCMat = idx-1;
        if (idx == 1){
          for (PetscInt k = 0; k < 2; k++){
            col_GlobMat[k] = idx+k;
            ierr = MatGetValues(mass, 1, &idx, 1, &col_GlobMat[k], &val[k]); CHKERRQ(ierr);
            col_BCMat[k] = col_GlobMat[k] - 1;
          }
          ierr = MatSetValues(A, 1, &row_BCMat, 2, col_BCMat, val, INSERT_VALUES); CHKERRQ(ierr);
        }
        else if (idx == N-2){
          for (PetscInt k = 0; k < 2; k++){
            col_GlobMat[k] = (idx-1)+k;
            ierr = MatGetValues(mass, 1, &idx, 1, &col_GlobMat[k], &val[k]); CHKERRQ(ierr);
            col_BCMat[k] = col_GlobMat[k] - 1;
          }
          ierr = MatSetValues(A, 1, &row_BCMat, 2, col_BCMat, val, INSERT_VALUES); CHKERRQ(ierr);
        }
        else{
          for (PetscInt k = 0; k < 3; k++) {
            col_GlobMat[k] = (idx - 1) + k;
            ierr = MatGetValues(mass, 1, &idx, 1, &col_GlobMat[k], &val[k]); CHKERRQ(ierr);
            col_BCMat[k] = col_GlobMat[k] - 1;
          }
          ierr = MatSetValues(A, 1, &row_BCMat, 3, col_BCMat, val, INSERT_VALUES);
          CHKERRQ(ierr);
        }
      }
      /* assemble matrix A */
      ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
      ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr); 
      
      /* extract relevant information from global rhsf and insert into vector b */
      ii = new PetscInt[N-2];
      jj = new PetscInt[N-2];
      for (int idx = 0; idx < N-2; idx++) {
        ii[idx] = idx + 1;
        jj[idx] = idx;
      }
      PetscScalar VecVals[N - 2];
      ierr = VecGetValues(rhsf, N-2, ii, VecVals); CHKERRQ(ierr);
      ierr = VecSetValues(b, N-2, jj, VecVals, ADD_VALUES); CHKERRQ(ierr);
      ierr = VecDuplicate(b, &xVec); CHKERRQ(ierr);

      /* assemble linear system and RHS */
      ierr = VecAssemblyBegin(b); CHKERRQ(ierr);
      ierr = VecAssemblyEnd(b); CHKERRQ(ierr);
      
      /* Solve the linear algebra */
      ierr = DoLinearAlgebra(); CHKERRQ(ierr);

      /* combine numerical solution and BCs */
      ierr = VecGetValues(xVec, N-2, jj, VecVals); CHKERRQ(ierr);
      ierr = VecSetValues(soln, N-2, ii, VecVals, INSERT_VALUES); CHKERRQ(ierr);
      // printf(YELLOW "\nNumerical Solution w/ BCs" RESET "\n");
      // ierr = VecView(soln, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    }
    /* Check numerical solution */
    ierr = L2Error(); CHKERRQ(ierr);

    // all petsc based functions need to end with PetscFinalize()
    ierr = PetscFinalize();
    return ierr;
}

PetscErrorCode BVP::L2Error(){
  /* 
  Check numerical solution 
  */
  // printf(YELLOW "Exact and Numerical Solution" RESET "\n");
  // PetscScalar uh, val;
  // for (PetscInt idx = 0; idx < N; idx++){
  //   val = u(info.xnod[idx]);
  //   ierr = VecSetValue(ExactSoln, idx, val, INSERT_VALUES); CHKERRQ(ierr);
  //   ierr = VecGetValues(soln, 1, &idx, &uh);
  //   ierr = PetscPrintf(PETSC_COMM_WORLD, "%.4e\t%.4e\n", double(val), double(uh));
  // }
  PetscScalar uval, duval;
  PetscScalar uhval, duhval;
  PetscInt mynum;
  PetscScalar tmp_soln, tmp_b, tmp_bp;
  l2Err = 0.0;
  h1Err = 0.0;
  QuadParams1D qps1d;
  get1D_QPs(info.maxord, qps1d);
  PetscShapeFunction1D shape1d;
  double xL, xR, dx, x;
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
      ierr = AssignEvaldBasis(dx, shape1d); CHKERRQ(ierr);
      /* evaluate known functions */
      uval = u(x);
      duval = up(x);
      uhval = 0.0;
      duhval = 0.0;
      for (int k = 0; k < ord; k++){
        mynum = info.nod[elem][k];
        ierr = VecGetValues(soln, 1, &mynum, &tmp_soln); CHKERRQ(ierr);
        ierr = VecGetValues(shape1d.psi, 1, &k, &tmp_b); CHKERRQ(ierr);
        ierr = VecGetValues(shape1d.dpsi, 1, &k, &tmp_bp); CHKERRQ(ierr);
        uhval += tmp_soln * tmp_b;
        duhval += tmp_soln * tmp_bp / dx;
      }
      l2Err += pow(uval - uhval, 2.0) * qps1d.w[l1] * dx;
      h1Err += pow(duval - duhval, 2.0) * qps1d.w[l1] * dx;
    }
  }
  l2Err = sqrt(l2Err);
  h1Err = sqrt(l2Err + h1Err);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "L2 Norm of Error %g\n", (double)l2Err); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "H1 Norm of Error %g\n", (double)h1Err); CHKERRQ(ierr);
  return ierr;
}

PetscErrorCode BVP::DoLinearAlgebra(){
  /* print out combined mass matrix with BCs */
  // printf(YELLOW "\nGlobal Mass" RESET "\n");
  // ierr = MatView(A, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  // printf(YELLOW "\nTotal RHSF" RESET "\n");
  // ierr = VecView(b, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  // exit(-1);

  /* Create KSP object */
  ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
  
  /* Set up solvers */
  ierr = KSPSetOperators(ksp, A, A); // setting Amat = Pmat is typical - see 73/278 of petsc user manual
  // set linear solver defaults for the problem
  ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
  ierr = PCSetType(pc, PCJACOBI); CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp, 1e-12, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT); CHKERRQ(ierr);
  // set runtime options THESE WILL OVERRIDE THE ABOVE THREE COMMANDS
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
  
  /* SOLVE! */
  ierr = KSPSolve(ksp, b, xVec); CHKERRQ(ierr);
  ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr); // view solver info
  
  /* clear up ksp */
  ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
  return ierr;
}

PetscErrorCode BVP::AssignLocalToGlobal(const std::vector<int> &tmp) {
  /*
  function to assign local FE matrices and rhsf vector to global 
  linear system and RHS
  */
  // get indices for global matrix
  ii = new PetscInt [ord];
  jj = new PetscInt [ord];
  for (int idx = 0; idx < ord; idx++) {
    ii[idx] = tmp[idx];
    jj[idx] = tmp[idx];
  }
  PetscScalar MatVals[ord*ord], VecVals[ord]; // to hold matrix/vector values in c arrays
  // get values and assign for mass...
  ierr = MatGetValues(m, ord, i, ord, j, MatVals); CHKERRQ(ierr);
  ierr = MatSetValues(mass, ord, ii, ord, jj, MatVals, ADD_VALUES); CHKERRQ(ierr);
  // ... stiffness ...
  ierr = MatGetValues(k, ord, i, ord, j, MatVals); CHKERRQ(ierr);
  ierr = MatSetValues(stiff, ord, ii, ord, jj, MatVals, ADD_VALUES); CHKERRQ(ierr);
  // ... and rhsf ...
  ierr = VecGetValues(f, ord, i, VecVals); CHKERRQ(ierr);
  ierr = VecSetValues(rhsf, ord, ii, VecVals, ADD_VALUES); CHKERRQ(ierr);

  return ierr;
}

PetscErrorCode BVP::InitializeLocalMatrices(){
  /*
  used to initialize local FE matrices. Sets everything to 0.0. 
  */
  /* set local FE matrices m and k, row by row */
  for (PetscInt idx = 0; idx < ord; idx++){
    ierr = MatSetValues(m, 1, &idx, ord, j, zero, INSERT_VALUES); CHKERRQ(ierr);
    ierr = MatSetValues(k, 1, &idx, ord, j, zero, INSERT_VALUES); CHKERRQ(ierr);
  }
  ierr = VecSet(f, 0.0); CHKERRQ(ierr);
  /* assemble matrices and vector */
  ierr = MatAssemblyBegin(m, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(m, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyBegin(k, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(k, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = VecAssemblyBegin(f); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(f); CHKERRQ(ierr);
  return ierr;
}

PetscErrorCode BVP::AssignEvaldBasis(const double dx, const PetscShapeFunction1D &shape1d){
    /*
    function used to take outer products of basis functions
    for the local matrices over each FE
    */
    PetscScalar x = 1.0/pow(dx, 2.0);
    PetscScalar v[ord];
    /* get all vector entries (vector of length ord) and store them in v*/ 
    ierr = VecGetValues(shape1d.psi, ord, i, v); CHKERRQ(ierr);
    /* insert vector, v, into -> (2 rows, rows 0 and 1, 1 col, cols 0, ...)*/
    ierr = MatSetValues(psi, ord, i, 1, &j[0], v, INSERT_VALUES); CHKERRQ(ierr);
    ierr = MatSetValues(psi, ord, i, 1, &j[1], zero, INSERT_VALUES); CHKERRQ(ierr);
    /* repeat process for dpsi matrix  */
    ierr = VecGetValues(shape1d.dpsi, ord, i, v); CHKERRQ(ierr);
    ierr = MatSetValues(dpsi, ord, i, 1, &j[0], v, INSERT_VALUES); CHKERRQ(ierr);
    ierr = MatSetValues(dpsi, ord, i, 1, &j[1], zero, INSERT_VALUES); CHKERRQ(ierr);

    /* assemle matrices */
    ierr = MatAssemblyBegin(psi, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(psi, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyBegin(dpsi, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(dpsi, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    /* generate local matrices, psiMat and dpsiMat */
    ierr = MatMatTransposeMult(psi, psi, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &psiMat); CHKERRQ(ierr);
    ierr = MatMatTransposeMult(dpsi, dpsi, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &dpsiMat); CHKERRQ(ierr);
    ierr = MatScale(dpsiMat, x); CHKERRQ(ierr);

    ierr = MatAssemblyBegin(psiMat, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(psiMat, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyBegin(dpsiMat, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(dpsiMat, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

    /* View/print to screen matrices
       format:
       row 0: (col #, <value>) (col #, <value>)
       row 1: (col #, <value>) (col #, <value>)
    */
    // printf(YELLOW "\npsi" RESET "\n");
    // ierr = MatView(psi, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    // printf(YELLOW "\ndpsi" RESET "\n");
    // ierr = MatView(dpsi, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    // printf(YELLOW "\npsiMat" RESET "\n");
    // ierr = MatView(psiMat, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    // printf(YELLOW "\ndpsiMat" RESET "\n");
    // ierr = MatView(dpsiMat, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
    // exit(-1);

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