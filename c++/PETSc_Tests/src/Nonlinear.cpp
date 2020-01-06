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

    /* Initialize snes */
    // create context
    ierr = SNESCreate(PETSC_COMM_WORLD, &snes); CHKERRQ(ierr);
    // create vectors for solution and residual
    ierr = VecCreate(PETSC_COMM_WORLD, &x); CHKERRQ(ierr);
    ierr = VecSetSizes(x, PETSC_DECIDE, 4); CHKERRQ(ierr);
    ierr = VecSetFromOptions(x); CHKERRQ(ierr);
    ierr = VecDuplicate(x, &r); CHKERRQ(ierr);
    /*
     *   Set SNES/KSP/KSP/PC runtime options, e.g.,
     *       -snes_view -snes_monitor -ksp_type <ksp> -pc_type <pc>
     *   These options will override those specified above as long as
     *   SNESSetFromOptions() is called _after_ any other customization
     *   routines.
     */
    // create Jacobian matrix data structure
    ierr = MatCreate(PETSC_COMM_WORLD,&J);CHKERRQ(ierr);
    ierr = MatSetSizes(J,PETSC_DECIDE,PETSC_DECIDE,4,4);CHKERRQ(ierr);
    ierr = MatSetFromOptions(J);CHKERRQ(ierr);
    ierr = MatSetUp(J);CHKERRQ(ierr);
    // set function evaluation routing and vector
    ierr = SNESSetFunction(snes, r, FormFunction, &ctx); CHKERRQ(ierr);
    // set Jacobian matrix data strcuture and evaluation routine
    ierr = SNESSetJacobian(snes, J, J, FormJacobian, &ctx);CHKERRQ(ierr);
    /*
     * Customize nonlinear solver; set runtime options
     *   Set linear solver defaults for this problem. By extracting the
     *   KSP and PC contexts from the SNES context, we can then
     *   directly call any KSP and PC routines to set various options.
     */
    // set linear solver defaults for the problem
    ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
    ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
    ierr = PCSetType(pc, PCNONE); CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp, 1e-12, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT); CHKERRQ(ierr);
    // set runtime options THESE WILL OVERRIDE THE ABOVE THREE COMMANDS
    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
    ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);

    /* Initialize local rhs vectors */
    ierr = VecCreate(PETSC_COMM_WORLD, &ctx.mass_src); CHKERRQ(ierr);
    ierr = VecSetSizes(ctx.mass_src, PETSC_DECIDE, info.order[0]); CHKERRQ(ierr);
    ierr = VecSetFromOptions(ctx.mass_src); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.mass_src, &ctx.momen_src); CHKERRQ(ierr);
    ierr = VecCopy(ctx.mass_src, ctx.momen_src); CHKERRQ(ierr);
    ierr = InitializeLocalRHSF(); CHKERRQ(ierr);
    /* Initialize global solution vectors */
    ierr = VecCreate(PETSC_COMM_WORLD, &velocity); CHKERRQ(ierr);
    ierr = VecSetSizes(velocity, PETSC_DECIDE, info.nnodes); CHKERRQ(ierr);
    ierr = VecSetFromOptions(velocity); CHKERRQ(ierr);
    ierr = VecDuplicate(velocity, &density); CHKERRQ(ierr);
    /* Assign left side BC */
    ierr = VecSetValue(velocity, 0, u(info.bounds[0]), INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecSetValue(density, 0, rho(info.bounds[0]), INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(velocity); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(velocity); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(density); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(density); CHKERRQ(ierr);

    QuadParams1D qps1d;
    get1D_QPs(info.maxord, qps1d);

    PetscShapeFunction1D shape1d;
    Vec TmpEvaldShape;
    ierr = VecCreate(PETSC_COMM_WORLD, &shape1d.psi); CHKERRQ(ierr);
    ierr = VecCreate(PETSC_COMM_WORLD, &shape1d.dpsi); CHKERRQ(ierr);
    ierr = VecSetSizes(shape1d.psi, PETSC_DECIDE, info.order[0]); CHKERRQ(ierr);
    ierr = VecSetSizes(shape1d.dpsi, PETSC_DECIDE, info.order[0]); CHKERRQ(ierr);
    ierr = VecSetFromOptions(shape1d.psi); CHKERRQ(ierr);
    ierr = VecSetFromOptions(shape1d.dpsi); CHKERRQ(ierr);
    ierr = VecDuplicate(shape1d.psi, &TmpEvaldShape);CHKERRQ(ierr);

    PetscScalar src_mass, src_momen;
    double xL, xR, dx, x;
    PetscScalar tmp_vel, tmp_rho;
    /* Sweep over elements and solve */
    for (int elem = 0; elem < info.nels; elem ++){
        xL = info.xnod[ info.nod[elem][0] ];
        xR = info.xnod[ info.nod[elem][info.order[elem]-1] ];
        dx = (xR-xL)/2.0;
        
        for (int l1 = 0; l1 < qps1d.nw; l1++){
            /* map from ref elem to real elem */
            x = xL + (1.0 + qps1d.xw[l1]) * dx;
            /* evaluate basis functions */
            ierr = PetscEvalBasis1D(qps1d.xw[l1], info.order[elem], shape1d); CHKERRQ(ierr);CHKERRQ(ierr);
            ierr = VecCopy(shape1d.psi, TmpEvaldShape);CHKERRQ(ierr);
            /* evaluate known functions and integrate */
            src_mass = MMS_Src_Mass(x);
            src_momen = MMS_Src_Momentum(x);
            ierr = VecScale(shape1d.psi, src_mass * qps1d.w[l1] * dx); CHKERRQ(ierr);
            ierr = VecScale(TmpEvaldShape, src_momen * qps1d.w[l1] * dx); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.mass_src, 1.0, shape1d.psi); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.momen_src, 1.0, TmpEvaldShape); CHKERRQ(ierr);
        }
        /* get upwind info */
        if (elem == 0) { // use BC info
          tmp_vel = u(info.bounds[0]);
          tmp_rho = rho(info.bounds[0]);
        } 
        else { // use previous element info
          index = elem * info.order[elem] - 1;
          VecGetValues(velocity, 1, &index, &tmp_vel);
          VecGetValues(density, 1, &index, &tmp_rho);
        }
        // compute upwind values for mass and momentum
        ctx.mass_upwind = tmp_vel * tmp_rho;
        ctx.momen_upwind = tmp_rho * pow(tmp_vel, 2.0);
        
        /* Solve nonlinear system of equations over elem */
        ierr = NLSolve(elem); CHKERRQ(ierr);
        
        /* set vector entries equal to zero (while maintaining structure) */
        ierr = InitializeLocalRHSF(); CHKERRQ(ierr);
    }
    /* print global solution */
    PetscPrintf(PETSC_COMM_WORLD, "cm\tvelocity\tdensity\n");
    for (int index = 0; index < info.nnodes; index++){
        VecGetValues(velocity, 1, &index, &tmp_vel);
        VecGetValues(density, 1, &index, &tmp_rho);
        PetscPrintf(PETSC_COMM_WORLD, "%.4e\t% .8e\t% .8e\t% .8e\t% .8e\n", info.xnod[index], tmp_vel, tmp_rho, u(info.xnod[index]), rho(info.xnod[index]) );
    }
    ierr = PetscFinalize();
    return ierr;
}

PetscErrorCode NonLinear::InitializeLocalRHSF(){
    /*
    used to initialize local FE vecters. Sets everything to 0.0. 
    */
    ierr = VecSet(ctx.mass_src, 0.0); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.mass_src); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.mass_src); CHKERRQ(ierr);
    ierr = VecSet(ctx.momen_src, 0.0); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.momen_src); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.momen_src); CHKERRQ(ierr);
    return ierr;
}

PetscErrorCode NonLinear::NLSolve(const int elem){
    /*
     *  create nonlinear solver context and solve equations
     */
    /* Set initial guess */
    ierr = VecSet(x, 1.0);CHKERRQ(ierr);

    /* Solve nonlinear system */
    ierr = SNESSolve(snes, NULL, x); CHKERRQ(ierr);
    // ierr = SNESView(snes, PETSC_VIEWER_STDOUT_WORLD);

    /* Map solution to global solution */
    PetscScalar value[4];         // array of values computed from NL solve
    PetscInt idx[4] = {0,1,2,3};  // indices of values to pull from NL solve
    PetscInt el = elem * 2;
    ierr = VecGetValues(x, 4, idx, value); CHKERRQ(ierr);
    ierr = VecSetValue(velocity, el, value[0], INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecSetValue(velocity, el+1, value[1], INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecSetValue(density, el, value[2], INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecSetValue(density, el+1, value[3], INSERT_VALUES); CHKERRQ(ierr);
    // PetscPrintf(PETSC_COMM_WORLD, "%.4e\t% .8e\t% .8e\n", info.xnod[el], value[0], value[2]);
    // PetscPrintf(PETSC_COMM_WORLD, "%.4e\t% .8e\t% .8e\n", info.xnod[el+1], value[1], value[3]);
    // exit(-1);
    return ierr;
}

PetscErrorCode FormFunction(SNES snes, Vec x, Vec f, void *ctx) {
    /*
     * evaluates nonlinear function, F(x)
     * INPUTS:
     *    snes - the PETSc SNES context
     *    x    - input vector
     *    ctx  - optional user-defined content
     * OUTPUTS:
     *    f    - function vector
     */
    // printf("  In FormFunction()\n");
    ApplicationCTX *user = (ApplicationCTX *)ctx;
    PetscErrorCode ierr;
    const PetscScalar *xx;
    PetscScalar       *ff;
    PetscScalar mass_src[2];
    PetscScalar momen_src[2];
    PetscInt idx[2] = {0, 1};
    /* Get values from local rhsf vectors and store into petsc scalars */
    ierr = VecGetValues(user->mass_src, 2, idx, mass_src);
    ierr = VecGetValues(user->momen_src, 2, idx, momen_src);

    // printf("  user->mass_upwind = %.3e\n", user->mass_upwind);
    // printf("  user->mass_src[0] = %.3e\n", mass_src[0]);
    // printf("  user->mass_src[1] = %.3e\n", mass_src[1]);
    // printf("  user->momen_upwind = %.3e\n", user->momen_upwind);
    // printf("  user->momen_src[0] = %.3e\n", momen_src[0]);
    // printf("  user->momen_src[1] = %.3e\n", momen_src[1]);
    // exit(-1);

    /*
    Get pointers to vector data.
        - For default PETSc vectors, VecGetArray() returns a pointer to
            the data array.  Otherwise, the routine is implementation dependent.
        - You MUST call VecRestoreArray() when you no longer need access to
            the array.
    */
    ierr = VecGetArrayRead(x,&xx);CHKERRQ(ierr);
    ierr = VecGetArray(f,&ff);CHKERRQ(ierr);

    /* Compute function */
    ff[0] = (xx[0]*xx[2])/3.0 + (xx[0]*xx[3])/6.0 + (xx[1]*xx[2])/6.0 + (xx[1]*xx[3])/3.0 - user->mass_upwind - mass_src[0];
    ff[1] =-(xx[0]*xx[2])/3.0 - (xx[0]*xx[3])/6.0 - (xx[1]*xx[2])/6.0 - (xx[1]*xx[3])/3.0 + xx[1]*xx[3] - mass_src[1];
    ff[2] = (xx[2]*pow(xx[0],2.0))/4.0 + (xx[3]*pow(xx[0],2.0))/12.0 + (xx[2]*xx[0]*xx[1])/6.0 + (xx[3]*xx[0]*xx[1])/6.0 + (xx[2]*pow(xx[1],2.0))/12.0 + (xx[3]*pow(xx[1],2.0))/4.0 - user->momen_upwind - momen_src[0];
    ff[3] =-(xx[2]*pow(xx[0],2.0))/4.0 - (xx[3]*pow(xx[0],2.0))/12.0 - (xx[2]*xx[0]*xx[1])/6.0 - (xx[3]*xx[0]*xx[1])/6.0 - (xx[2]*pow(xx[1],2.0))/12.0 - (xx[3]*pow(xx[1],2.0))/4.0 + xx[3]*pow(xx[1],2.0) - momen_src[1];

    /* Restore vectors */
    ierr = VecRestoreArrayRead(x,&xx);CHKERRQ(ierr);
    ierr = VecRestoreArray(f,&ff);CHKERRQ(ierr);

    return 0;
}

PetscErrorCode FormJacobian(SNES snes, Vec x, Mat jac, Mat B, void *dummy) {
    /*
     *  compute Jacobia entries and insert into matrix
     */
    const PetscScalar *xx;
    PetscScalar A[16];
    PetscErrorCode ierr;
    PetscInt idx[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

    /* Get pointer to vector data */
    ierr = VecGetArrayRead(x,&xx);CHKERRQ(ierr);

    /* Compute Jacobian entries */
    // over f0
    A[0] = xx[2]/3.0 + xx[3]/6.0;
    A[1] = xx[2]/6.0 + xx[3]/3.0;
    A[2] = xx[0]/3.0 + xx[1]/6.0;
    A[3] = xx[0]/6.0 + xx[1]/3.0;
    // over f1
    A[4] = -A[0];         //-xx[2]/3.0 - xx[3]/6.0;
    A[5] = -A[1] + xx[3]; //-xx[2]/6.0 - xx[3]/3.0 + xx[3];
    A[6] = -A[2];         //-xx[0]/3.0 - xx[1]/6.0;
    A[7] = -A[3] + xx[1]; //-xx[0]/6.0 - xx[1]/3.0 + xx[1];
    // over f2
    A[8] = (2.0*xx[0]*xx[2])/4.0 + (2.0*xx[0]*xx[3])/12.0 + (xx[2]*xx[1])/6.0 + (xx[3]*xx[1])/6.0;
    A[9] = (2.0*xx[1]*xx[3])/4.0 + (2.0*xx[1]*xx[2])/12.0 + (xx[3]*xx[0])/6.0 + (xx[2]*xx[0])/6.0;
    A[10] = pow(xx[0], 2.0)/4.0 + (xx[0]*xx[1])/6.0 + pow(xx[1],2.0)/12.0;
    A[11] = pow(xx[0], 2.0)/12.0 + (xx[0]*xx[1])/6.0 + pow(xx[1], 2.0)/4.0;
    // over f3
    A[12] = -A[8];                     //-(2.0*xx[0]*xx[2])/4.0 - (2.0*xx[0]*xx[3])/12.0 - (xx[2]*xx[1])/6.0 - (xx[3]*xx[1])/6.0;
    A[13] = -A[9] + 2.0*xx[1]*xx[3] ;  //-(2.0*xx[1]*xx[3])/4.0 - (2.0*xx[1]*xx[2])/12.0 - (xx[3]*xx[0])/6.0 - (xx[2]*xx[0])/6.0 + 2.0*xx[1]*xx[3];
    A[14] = -A[10];                    //-pow(xx[0], 2.0)/4.0 - (xx[0]*xx[1])/6.0 - pow(xx[0],2.0)/12.0;
    A[15] = -A[11] + pow(xx[1], 2.0);  //-pow(xx[0], 2.0)/12.0 - (xx[0]*xx[1])/6.0 - pow(xx[1], 2.0)/4.0 + pow(xx[1], 2.0);

    /* and insert into matrix B */
    ierr = MatSetValues(B, 4, idx, 4, idx, A, INSERT_VALUES); CHKERRQ(ierr);

    /* Restor vector */
    ierr = VecRestoreArrayRead(x, &xx);CHKERRQ(ierr);
    
    /* Assemble matrix */
    ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    if (jac != B) {
        ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    }
    return 0;
}