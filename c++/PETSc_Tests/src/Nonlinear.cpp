#include "Nonlinear.hpp"

PetscErrorCode NonLinear::Initialize_NL_1D(){
    /*
     *  initialize all PETSC data structures for NonLinear::NL_1D
     */
    if (info.nnodes == NAN){
        SETERRQ(PETSC_COMM_WORLD,1,"Please build mesh via SetUpGrid::add_DFEMGrid before executing Initialize_NL_1D()!");
    }
    N = info.nnodes * 3;
    /* ------ Initialize SNES ------ */
    ierr = SNESCreate(PETSC_COMM_WORLD, &snes); CHKERRQ(ierr); // create context 
    ierr = SNESSetTolerances(snes, 1e-08, 1e-50, 1e-50, PETSC_DEFAULT, PETSC_DEFAULT); CHKERRQ(ierr);
    // create solution vector
    ierr = VecCreate(PETSC_COMM_WORLD, &soln); CHKERRQ(ierr);
    ierr = VecSetSizes(soln, PETSC_DECIDE, N); CHKERRQ(ierr);
    ierr = VecSetFromOptions(soln); CHKERRQ(ierr);
    ierr = VecDuplicate(soln, &residual); CHKERRQ(ierr);
    // create Jacobian matrix data structure
    ierr = MatCreate(PETSC_COMM_WORLD, &J); CHKERRQ(ierr);
    ierr = MatSetSizes(J, PETSC_DECIDE, PETSC_DECIDE, N, N); CHKERRQ(ierr);
    ierr = MatSetFromOptions(J); CHKERRQ(ierr);
    ierr = MatSetUp(J); CHKERRQ(ierr);
    // set function and jacobian functions
    ierr = SNESSetFunction(snes, residual, FormFunction, this); CHKERRQ(ierr); // set function evaluation routing and vector 
    ierr = SNESSetJacobian(snes, J, J, FormJacobian, this); CHKERRQ(ierr); // set Jacobian matrix data strcuture and evaluation routine 
    /*
    *  Customize nonlinear solver; set runtime options
    *    Set linear solver defaults for this problem. By extracting the
    *    KSP and PC contexts from the SNES context, we can then
    *    directly call any KSP and PC routines to set various options.
    *   Set SNES/KSP/KSP/PC runtime options, e.g.,
    *       -snes_view -snes_monitor -ksp_type <ksp> -pc_type <pc>
    *   These options will override those specified above as long as
    *   SNESSetFromOptions() is called _after_ any other customization
    *   routines.
    */
    // set linear solver defaults for the problem
    ierr = SNESGetKSP(snes,&ksp); CHKERRQ(ierr);
    ierr = KSPSetType(ksp, KSPGMRES); CHKERRQ(ierr);
    ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
    ierr = PCSetType(pc, PCILU); CHKERRQ(ierr);
    ierr = KSPSetTolerances(ksp, 1e-12, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT); CHKERRQ(ierr);
    // set runtime options THESE WILL OVERRIDE THE ABOVE THREE COMMANDS
    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);
    ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);

    // ------ local FE vec initialization ------
    ierr = VecCreate(PETSC_COMM_WORLD, &mass_basis_src); CHKERRQ(ierr);
    ierr = VecSetSizes(mass_basis_src, PETSC_DECIDE, info.order[0]); CHKERRQ(ierr);
    ierr = VecSetFromOptions(mass_basis_src); CHKERRQ(ierr);
    ierr = VecDuplicate(mass_basis_src, &momen_basis_src); CHKERRQ(ierr);
    ierr = VecDuplicate(mass_basis_src, &efluid_basis_src); CHKERRQ(ierr);
    
    ierr = VecCreate(PETSC_COMM_WORLD, &ctx.loc_mass_src); CHKERRQ(ierr);
    ierr = VecSetSizes(ctx.loc_mass_src, PETSC_DECIDE, info.order[0]); CHKERRQ(ierr);
    ierr = VecSetFromOptions(ctx.loc_mass_src); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_src, &ctx.loc_momen_src); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_src, &ctx.loc_efluid_src); CHKERRQ(ierr);

    // ------ Global FE vec initialization ------
    // source information
    ierr = VecCreate(PETSC_COMM_WORLD, &ctx.glo_mass_src); CHKERRQ(ierr);
    ierr = VecSetSizes(ctx.glo_mass_src, PETSC_DECIDE, info.nnodes); CHKERRQ(ierr);
    ierr = VecSetFromOptions(ctx.glo_mass_src); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_src, &ctx.glo_momen_src); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_src, &ctx.glo_efluid_src); CHKERRQ(ierr);
    // solution vectors
    ierr = VecDuplicate(ctx.glo_mass_src, &velocity); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_src, &density); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_src, &energy); CHKERRQ(ierr);

    PetscViewerDrawOpen(PETSC_COMM_WORLD, 0,0,0,0,400,400, &monCTX.viewer1);
    PetscViewerDrawOpen(PETSC_COMM_WORLD, 0,0,PETSC_DECIDE,PETSC_DECIDE,400,400, &monCTX.viewer2);
    PetscViewerDrawOpen(PETSC_COMM_WORLD, 0,0,PETSC_DECIDE,PETSC_DECIDE,400,400, &monCTX.viewer3);
    SNESMonitorSet(snes, Monitor, this,0);

    return ierr;
}

PetscErrorCode NonLinear::NL_1D(){
    /*
    1D nonlinear solver using PETSc snes functionality
    */
    if (size != 1){
        SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only!");
    }
    ierr = Initialize_NL_1D(); CHKERRQ(ierr);

    /* --------- Assign and Assembly Boundary Conditions --------- */
    // assign left side BC (hydro equations)
    ierr = VecSetValue(velocity, 0, u(info.bounds[0]), INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecSetValue(density, 0, rho(info.bounds[0]), INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecSetValue(energy, 0, efluid(info.bounds[0]), INSERT_VALUES); CHKERRQ(ierr);
    // assemble vectors
    ierr = VecAssemblyBegin(velocity); CHKERRQ(ierr); 
    ierr = VecAssemblyEnd(velocity); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(density); CHKERRQ(ierr); 
    ierr = VecAssemblyEnd(density); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(energy); CHKERRQ(ierr); 
    ierr = VecAssemblyEnd(energy); CHKERRQ(ierr);   

    ierr = InitializeLocalRHSF(); CHKERRQ(ierr);
    get1D_QPs(info.maxord, qps1d); // set qps

    /* Sweep over elements and solve */
    for (int elem = 0; elem < info.nels; elem ++){
        xL = info.xnod[ info.nod[elem][0] ];
        xR = info.xnod[ info.nod[elem][info.order[elem]-1] ];
        dx = (xR-xL)/2.0;

        for (int l1 = 0; l1 < qps1d.nw; l1++){
            /* map from ref elem to real elem */
            x = xL + (1.0 + qps1d.xw[l1]) * dx;
            /* evaluate basis functions */
            ierr = EvalBasis(qps1d.xw[l1], info.order[elem]); CHKERRQ(ierr);
            /* evaluate known functions and integrate */
            src_mass = MMS_Src_Mass(x);
            src_momen = MMS_Src_Momentum(x);
            src_energy = MMS_Src_Energy(x);

            ierr = VecAXPY(ctx.loc_mass_src, src_mass * qps1d.w[l1] * dx, mass_basis_src); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_momen_src, src_momen * qps1d.w[l1] * dx, momen_basis_src); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_efluid_src, src_energy * qps1d.w[l1] * dx, efluid_basis_src); CHKERRQ(ierr);
        }

        /* Map local vectors to global vectors */
        ierr = Local2Global(elem); CHKERRQ(ierr);

        /* set vector entries equal to zero (while maintaining structure) */
        ierr = InitializeLocalRHSF(); CHKERRQ(ierr);
    }
    /* Solve nonlinear system of equations over elem */
    ierr = NLSolve(); CHKERRQ(ierr);
    /* print global solution */
    PetscScalar tmp_vel, tmp_rho, tmp_efluid;
    PetscPrintf(PETSC_COMM_WORLD, "cm\tvelocity\tdensity\n");
    for (int index = 0; index < info.nnodes; index++){
        VecGetValues(velocity, 1, &index, &tmp_vel);
        VecGetValues(density, 1, &index, &tmp_rho);
        VecGetValues(energy, 1, &index, &tmp_efluid);
        PetscPrintf(PETSC_COMM_WORLD, "%.4e\t% .8e\t% .8e\t% .8e\t% .8e\t% .8e\t% .8e\n", info.xnod[index], tmp_vel, tmp_rho, tmp_efluid, u(info.xnod[index]), rho(info.xnod[index]), efluid(info.xnod[index]) );
        if (index % 2 == 1){
            PetscPrintf(PETSC_COMM_WORLD,"\n");
        }
    }
    /* Check numerical solution */
    // ierr = VelRho_L2Error(); CHKERRQ(ierr);

    // all petsc based functions need to end with PetscFinalize()
    ierr = PetscFinalize();
    return ierr;
}

PetscErrorCode NonLinear::EvalBasis(const double x, const int ord){
    /* 
     *  - evaluates basis functions for FEM integrals. 
     *  - basis functions are defined on [-1,1]
     */
    if (ord == 2){
        PetscInt idx[2] = {0,1};
        PetscScalar b1, b1p, b2, b2p;
        // basis function evaluation
        b1 = 0.5 * (1.0 - x);
        b2 = 0.5 * (1.0 + x);
        b1p =-0.5;
        b2p = 0.5;
        PetscScalar src[] = {b1, b2};
        ierr = VecSetValues(mass_basis_src, 2, idx, src, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(momen_basis_src, 2, idx, src, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(efluid_basis_src, 2, idx, src, INSERT_VALUES); CHKERRQ(ierr);
    }
    else{
        PetscPrintf(PETSC_COMM_WORLD, "\n\nOrder = %d shape function does not exist.\n\n",ord);
        exit(1);
    }
    return ierr;
}

PetscErrorCode NonLinear::VelRho_L2Error(){
    /*
     *  Check numerical solution for velocity and density
     */
    PetscPrintf(PETSC_COMM_WORLD,"\n\nneed to add in energy error information\n\n");
    PetscScalar uval, duval, rval, drval, eval, deval;
    PetscScalar uhval, duhval, rhval, drhval, ehval, dehval;
    PetscInt mynum;
    PetscScalar tmp_vel, tmp_rho, tmp_em;
    PetscScalar tmp_b, tmp_bp;
    l2Err_Vel = 0.0;
    l2Err_Rho = 0.0;
    l2Err_Em = 0.0;
    h1Err_Vel = 0.0;
    h1Err_Rho = 0.0;
    h1Err_Em = 0.0;

    QuadParams1D qps1d;
    get1D_QPs(info.maxord, qps1d);
    
    PetscShapeFunction1D shape1d;
    ierr = VecCreate(PETSC_COMM_WORLD, &shape1d.psi); CHKERRQ(ierr);
    ierr = VecCreate(PETSC_COMM_WORLD, &shape1d.dpsi); CHKERRQ(ierr);
    ierr = VecSetSizes(shape1d.psi, PETSC_DECIDE, info.order[0]); CHKERRQ(ierr);
    ierr = VecSetSizes(shape1d.dpsi, PETSC_DECIDE, info.order[0]); CHKERRQ(ierr);
    ierr = VecSetFromOptions(shape1d.psi); CHKERRQ(ierr);
    ierr = VecSetFromOptions(shape1d.dpsi); CHKERRQ(ierr);

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
            /* evaluate known functions */
            uval = u(x);
            duval = up(x);
            rval = rho(x);
            drval = rhop(x);
            eval = efluid(x);
            deval = efluidp(x);
            uhval = 0.0;
            duhval = 0.0;
            rhval = 0.0;
            drhval = 0.0;
            ehval = 0.0;
            dehval = 0.0;
            for (int k = 0; k < info.order[0]; k++){
                mynum = info.nod[elem][k];
                ierr = VecGetValues(velocity, 1, &mynum, &tmp_vel); CHKERRQ(ierr);
                ierr = VecGetValues(density, 1, &mynum, &tmp_rho); CHKERRQ(ierr);
                ierr = VecGetValues(energy, 1, &mynum, &tmp_em); CHKERRQ(ierr);
                ierr = VecGetValues(shape1d.psi, 1, &k, &tmp_b); CHKERRQ(ierr);
                ierr = VecGetValues(shape1d.dpsi, 1, &k, &tmp_bp); CHKERRQ(ierr);
                uhval += tmp_vel * tmp_b;
                rhval += tmp_rho * tmp_b;
                ehval += tmp_em * tmp_b;
                duhval += tmp_vel * tmp_bp / dx;
                drhval += tmp_rho * tmp_bp / dx;
                dehval += tmp_em * tmp_bp / dx;
            }
            l2Err_Vel += pow(uval - uhval, 2.0) * qps1d.w[l1] * dx;
            l2Err_Rho += pow(rval - rhval, 2.0) * qps1d.w[l1] * dx;
            l2Err_Em += pow(eval - ehval, 2.0) * qps1d.w[l1] * dx;
            h1Err_Vel += pow(duval - duhval, 2.0) * qps1d.w[l1] * dx;
            h1Err_Rho += pow(drval - drhval, 2.0) * qps1d.w[l1] * dx;
            h1Err_Em += pow(deval - dehval, 2.0) * qps1d.w[l1] * dx;
        }
    }
    l2Err_Vel = sqrt(l2Err_Vel);
    l2Err_Rho = sqrt(l2Err_Rho);
    l2Err_Em = sqrt(l2Err_Em);
    h1Err_Vel = sqrt(l2Err_Vel + h1Err_Vel);
    h1Err_Rho = sqrt(l2Err_Rho + h1Err_Rho);
    h1Err_Em = sqrt(l2Err_Em + h1Err_Em);
    PetscPrintf(PETSC_COMM_WORLD, "%.8e\t%.8e\t%.8e\t%.8e\t\t%.8e\t%.8e\t%.8e\n", info.hel, l2Err_Vel, l2Err_Rho, l2Err_Em, h1Err_Vel, h1Err_Rho, h1Err_Em);
    return ierr;
}

PetscErrorCode NonLinear::InitializeLocalRHSF(){
    /*
    used to initialize local FE vecters. Sets everything to 0.0. 
    */
    ierr = VecSet(ctx.loc_mass_src, 0.0); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_mass_src); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_mass_src); CHKERRQ(ierr);
    // conservation of momentum
    ierr = VecSet(ctx.loc_momen_src, 0.0); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_momen_src); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_momen_src); CHKERRQ(ierr);
    // conservation of fluid energy
    ierr = VecSet(ctx.loc_efluid_src, 0.0); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_efluid_src); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_efluid_src); CHKERRQ(ierr);
    return ierr;
}

PetscErrorCode NonLinear::Local2Global(const int el){
    /*
     *  Maps local FE integration information to global system of nonlinear equations
     */
    PetscInt glo_idx[2];
    glo_idx[0] = el*info.order[el];
    glo_idx[1] = glo_idx[0] + (info.order[el]-1);
    PetscInt loc_idx[2] = {0,1};
    PetscScalar loc_val[2];
    ierr = VecGetValues(ctx.loc_mass_src, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_mass_src, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecGetValues(ctx.loc_momen_src, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_momen_src, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecGetValues(ctx.loc_efluid_src, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_efluid_src, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    return ierr;
}

PetscErrorCode Monitor(SNES snes, PetscInt its, PetscReal fnorm, void *ctx){
    NonLinear *user = (NonLinear *)ctx;
    Vec soln;
    PetscScalar *tmp_vel, *tmp_rho, *tmp_efluid;
    tmp_vel = (PetscScalar *)malloc(user->info.nnodes * sizeof(PetscScalar));
    tmp_rho = (PetscScalar *)malloc(user->info.nnodes * sizeof(PetscScalar));
    tmp_efluid = (PetscScalar *)malloc(user->info.nnodes * sizeof(PetscScalar));

    PetscInt *idu, *idr, *idem;
    idu = (PetscInt*) malloc(user->info.nnodes * sizeof(PetscInt));
    idr = (PetscInt *)malloc(user->info.nnodes * sizeof(PetscInt));
    idem = (PetscInt *)malloc(user->info.nnodes * sizeof(PetscInt));
    for (int i=0; i<user->info.nnodes; i++){
        idu[i] = i;
        idr[i] = i + user->info.nnodes;
        idem[i] = i + 2*user->info.nnodes;
    }

    PetscPrintf(PETSC_COMM_WORLD,"iter = %D, SNES Function norm %g\n",its,(double)fnorm);
    SNESGetSolution(snes, &soln);
    VecGetValues(soln, user->info.nnodes, idu, tmp_vel);
    VecGetValues(soln, user->info.nnodes, idr, tmp_rho);
    VecGetValues(soln, user->info.nnodes, idem, tmp_efluid);
    VecSetValues(user->velocity, user->info.nnodes, idu, tmp_vel, INSERT_VALUES);
    VecSetValues(user->density, user->info.nnodes, idu, tmp_rho, INSERT_VALUES);
    VecSetValues(user->energy, user->info.nnodes, idu, tmp_efluid, INSERT_VALUES);
    VecView(user->velocity, user->monCTX.viewer1);
    VecView(user->density, user->monCTX.viewer2);
    VecView(user->energy, user->monCTX.viewer3);
    // for (int index = 0; index < user->info.nnodes; index++){
    //     PetscPrintf(PETSC_COMM_WORLD, "%.4e\t% .8e\t% .8e\t% .8e\n", user->info.xnod[index], tmp_vel[index], tmp_rho[index], tmp_efluid[index]);
    //     if (index % 2 == 1){
    //         PetscPrintf(PETSC_COMM_WORLD,"\n");
    //     }
    // }

    return 0;
}

PetscErrorCode NonLinear::NLSolve(){
    /*
     *  create nonlinear solver context and solve equations
     */
    /* Set initial guess */
    for (PetscInt i=0; i<info.nnodes; i++){
        ierr = VecSetValue(soln, i, u(0.0), INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValue(soln, i+info.nnodes, rho(0.0), INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValue(soln, i+2*info.nnodes, efluid(0.0), INSERT_VALUES); CHKERRQ(ierr);
    }
    /* Solve nonlinear system */
    ierr = SNESSolve(snes, NULL, soln); CHKERRQ(ierr);
    ierr = SNESView(snes, PETSC_VIEWER_STDOUT_WORLD);
    
    /* Map computed solution to solution vectors */
    PetscScalar *tmp_vel, *tmp_rho, *tmp_efluid;
    tmp_vel = (PetscScalar*) malloc(info.nnodes * sizeof(PetscScalar));
    tmp_rho = (PetscScalar*) malloc(info.nnodes * sizeof(PetscScalar));
    tmp_efluid = (PetscScalar*) malloc(info.nnodes * sizeof(PetscScalar));

    PetscInt *idu, *idr, *idem;
    idu = (PetscInt*) malloc(info.nnodes * sizeof(PetscInt));
    idr = (PetscInt *)malloc(info.nnodes * sizeof(PetscInt));
    idem = (PetscInt *)malloc(info.nnodes * sizeof(PetscInt));
    for (int i=0; i<info.nnodes; i++){
        idu[i] = i;
        idr[i] = i + info.nnodes;
        idem[i] = i + 2*info.nnodes;
    }
    ierr = VecGetValues(soln, info.nnodes, idu, tmp_vel); CHKERRQ(ierr);
    ierr = VecGetValues(soln, info.nnodes, idr, tmp_rho); CHKERRQ(ierr);
    ierr = VecGetValues(soln, info.nnodes, idem, tmp_efluid); CHKERRQ(ierr);
    ierr = VecSetValues(velocity, info.nnodes, idu, tmp_vel, INSERT_VALUES);
    ierr = VecSetValues(density, info.nnodes, idu, tmp_rho, INSERT_VALUES);
    ierr = VecSetValues(energy, info.nnodes, idu, tmp_efluid, INSERT_VALUES);

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

    NonLinear *user = (NonLinear *)ctx;
    PetscErrorCode ierr;
    const PetscScalar *xx;
    PetscScalar       *ff;
    const int nn {user->info.nnodes};

    PetscScalar mass_upwind, *mass_src;
    mass_src = (PetscScalar *) malloc(nn * sizeof(PetscScalar));
    PetscScalar momen_upwind, *momen_src;
    momen_src = (PetscScalar *)malloc(nn * sizeof(PetscScalar));
    PetscScalar efluid_upwind, *efluid_src;
    efluid_src = (PetscScalar *)malloc(nn * sizeof(PetscScalar));

    // assign petsc Vec to c array for use
    PetscInt *idx;
    idx = (PetscInt*) malloc(nn * sizeof(PetscInt));
    for (int i=0; i<nn; i++){
        idx[i] = i;
    }
    ierr = VecGetValues(user->ctx.glo_mass_src, nn, idx, mass_src); CHKERRQ(ierr);
    ierr = VecGetValues(user->ctx.glo_momen_src, nn, idx, momen_src); CHKERRQ(ierr);
    ierr = VecGetValues(user->ctx.glo_efluid_src, nn, idx, efluid_src); CHKERRQ(ierr);

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
    /* loop over elements by looping over node 1 of each element
     * e.g. with liner FEs:
     *    cell 1 -> i = 0
     *    cell 2 -> i = 2
     */
    const double l{user->info.bounds[0]};
    for (PetscInt i=0; i<nn; i+=user->info.order[0]){
        if (i == 0){
            mass_upwind = user->rho(l) * user->u(l);
            momen_upwind = user->rho(l) * pow(user->u(l),2.0) + user->ctx.gamma_s * user->efluid(l);
            efluid_upwind = 0.5 * user->rho(l) * pow(user->u(l),3.0) + (1.0 + user->ctx.gamma_s)* user->u(l) * user->efluid(l);
        }
        else{
            mass_upwind = xx[nn+i-1] * xx[i-1];
            momen_upwind = xx[nn+i-1] * pow(xx[i-1],2.0) + user->ctx.gamma_s * xx[2*nn+i-1];
            efluid_upwind = 0.5 * xx[nn+i-1] * pow(xx[i-1],3.0) + (1.0 + user->ctx.gamma_s) * xx[i-1] * xx[2*nn+i-1];
        }
        // Conservation of mass
        ff[i]   =   xx[nn+i]*xx[i]     * (1.0/3.0)
                  + xx[nn+i+1]*xx[i]   * (1.0/6.0)
                  + xx[nn+i]*xx[i+1]   * (1.0/6.0)
                  + xx[nn+i+1]*xx[i+1] * (1.0/3.0)
                  - mass_upwind
                  - mass_src[i];
        ff[i+1] = - xx[nn+i]*xx[i]     * (1.0/3.0)
                  - xx[nn+i+1]*xx[i]   * (1.0/6.0)
                  - xx[nn+i]*xx[i+1]   * (1.0/6.0)
                  - xx[nn+i+1]*xx[i+1] * (1.0/3.0)
                  + xx[nn+i+1]*xx[i+1]
                  - mass_src[i+1];
        // Conservation of momentum
        ff[nn+i]  =   xx[nn+i]*pow(xx[i],2.0)     * (1.0/4.0)
                    + xx[nn+i+1]*pow(xx[i],2.0)   * (1.0/12.0)
                    + xx[nn+i]*xx[i]*xx[i+1]      * (1.0/6.0)
                    + xx[nn+i+1]*xx[i]*xx[i+1]    * (1.0/6.0)
                    + xx[nn+i]*pow(xx[i+1],2.0)   * (1.0/12.0)
                    + xx[nn+i+1]*pow(xx[i+1],2.0) * (1.0/4.0)
                    + xx[2*nn+i]                  * (user->ctx.gamma_s/2.0)
                    + xx[2*nn+i+1]                * (user->ctx.gamma_s/2.0)
                    - momen_upwind
                    - momen_src[i];
        ff[nn+i+1]= - xx[nn+i]*pow(xx[i],2.0)     * (1.0/4.0)
                    - xx[nn+i+1]*pow(xx[i],2.0)   * (1.0/12.0)
                    - xx[nn+i]*xx[i]*xx[i+1]      * (1.0/6.0)
                    - xx[nn+i+1]*xx[i]*xx[i+1]    * (1.0/6.0)
                    - xx[nn+i]*pow(xx[i+1],2.0)   * (1.0/12.0)
                    - xx[nn+i+1]*pow(xx[i+1],2.0) * (1.0/4.0)
                    - xx[2*nn+i]                  * (user->ctx.gamma_s/2.0)
                    - xx[2*nn+i+1]                * (user->ctx.gamma_s/2.0)
                    + xx[nn+i+1]*pow(xx[i+1],2.0) + (user->ctx.gamma_s*xx[2*nn+i+1])
                    - momen_src[i+1];
        // Conservation of fluid energy
        ff[2*nn+i]  =   xx[nn+i]*pow(xx[i],3.0)           * (1.0/10.0)
                      + xx[nn+i+1]*pow(xx[i],3.0)         * (1.0/40.0)
                      + xx[nn+i]*pow(xx[i],2.0)*xx[i+1]   * (3.0/40.0)
                      + xx[nn+i+1]*pow(xx[i],2.0)*xx[i+1] * (1.0/20.0)
                      + xx[nn+i]*xx[i]*pow(xx[i+1],2.0)   * (1.0/20.0)
                      + xx[nn+i+1]*xx[i]*pow(xx[i+1],2.0) * (3.0/40.0)
                      + xx[nn+i]*pow(xx[i+1],3.0)         * (1.0/40.0)
                      + xx[nn+i+1]*pow(xx[i+1],3.0)       * (1.0/10.0)
                      + xx[2*nn+i]*xx[i]                  * (1.0/3.0) * (1.0+user->ctx.gamma_s)
                      + xx[2*nn+i+1]*xx[i]                * (1.0/6.0) * (1.0+user->ctx.gamma_s)
                      + xx[2*nn+i]*xx[i+1]                * (1.0/6.0) * (1.0+user->ctx.gamma_s)
                      + xx[2*nn+i+1]*xx[i+1]              * (1.0/3.0) * (1.0+user->ctx.gamma_s)
                      - efluid_upwind
                      - efluid_src[i];
        ff[2*nn+i+1]= - xx[nn+i]*pow(xx[i],3.0)           * (1.0/10.0)
                      - xx[nn+i+1]*pow(xx[i],3.0)         * (1.0/40.0)
                      - xx[nn+i]*pow(xx[i],2.0)*xx[i+1]   * (3.0/40.0)
                      - xx[nn+i+1]*pow(xx[i],2.0)*xx[i+1] * (1.0/20.0)
                      - xx[nn+i]*xx[i]*pow(xx[i+1],2.0)   * (1.0/20.0)
                      - xx[nn+i+1]*xx[i]*pow(xx[i+1],2.0) * (3.0/40.0)
                      - xx[nn+i]*pow(xx[i+1],3.0)         * (1.0/40.0)
                      - xx[nn+i+1]*pow(xx[i+1],3.0)       * (1.0/10.0)
                      - xx[2*nn+i]*xx[i]                  * (1.0/3.0) * (1.0+user->ctx.gamma_s)
                      - xx[2*nn+i+1]*xx[i]                * (1.0/6.0) * (1.0+user->ctx.gamma_s)
                      - xx[2*nn+i]*xx[i+1]                * (1.0/6.0) * (1.0+user->ctx.gamma_s)
                      - xx[2*nn+i+1]*xx[i+1]              * (1.0/3.0) * (1.0+user->ctx.gamma_s)
                      + 0.5*xx[nn+i+1]*pow(xx[i+1],3.0) + (1.0+user->ctx.gamma_s)*(xx[i+1]*xx[2*nn+i+1])
                      - efluid_src[i+1];
    }
    /* Restore vectors */
    ierr = VecRestoreArrayRead(x,&xx);CHKERRQ(ierr);
    ierr = VecRestoreArray(f,&ff);CHKERRQ(ierr);

    return 0;
}

PetscErrorCode FormJacobian(SNES snes, Vec x, Mat jac, Mat B, void *ctx) {
    /*
     *  compute Jacobia entries and insert into matrix
     */
    PetscErrorCode ierr;
    NonLinear *user = (NonLinear *)ctx;
    const PetscScalar *xx;
    const int nn{user->N};
    const int n{user->info.nnodes};
    // const int nels{user->info.nels};

    PetscScalar *A;
    A = (PetscScalar *)calloc(nn*nn, sizeof(PetscScalar));

    // assign petsc Vec to c array for use
    PetscInt *idx;
    idx = (PetscInt*) malloc(nn*nn * sizeof(PetscInt));
    for (int i=0; i<nn*nn; i++){
        idx[i] = i;
    }

    /* Get pointer to vector data */
    ierr = VecGetArrayRead(x,&xx);CHKERRQ(ierr);

    PetscInt j, k;
    /* Compute Jacobian entries */
    /* loop over conservation of mass */
    k = 0;
    for (PetscInt i = 0; i < nn*(1*n); i+=nn){
        // over f_i
        A[i] = xx[n+k]/3.0 + xx[n+k+1]/6.0;  
        A[i+1] = xx[n+k]/6.0 + xx[n+k+1]/3.0;
        A[n+i] = xx[k]/3.0 + xx[k+1]/6.0;    
        A[n+i+1] = xx[k]/6.0 + xx[k+1]/3.0;  
        A[2*n+i] = 0.0;
        A[2*n+i+1] = 0.0;
        if (k > 0) {
            A[i-1] = -xx[n+k-1];
            A[n+i-1] = -xx[k-1];
        }
        // over f_{i+1}
        j = i;
        i += nn;
        A[i] = -A[j];
        A[i+1] = -A[j+1] + xx[n+k+1];
        A[n+i] = -A[n+j];
        A[n+i+1] = -A[n+j+1] + xx[k+1];
        A[2*n+i] = 0.0;
        A[2*n+i+1] = 0.0;
        i += 2;
        k += 2;
    }
    /* Conservation of Momentum */
    k = 0;
    for (PetscInt i = nn*(1*n); i < nn*(2*n); i+=nn){
        // over f_{6*n+i}
        A[i] = xx[n+k]*xx[k]/2.0 + xx[n+k+1]*xx[k]/6.0 + xx[n+k]*xx[k+1]/6.0 + xx[n+k+1]*xx[k+1]/6.0;  
        A[i+1] = xx[n+k]*xx[k]/6.0 + xx[n+k+1]*xx[k]/6.0 + xx[n+k]*xx[k+1]/6.0 + xx[n+k+1]*xx[k+1]/2.0;
        A[n+i] = pow(xx[k],2.0)/4.0 + xx[k]*xx[k+1]/6.0 + pow(xx[k+1],2.0)/12.0;                       
        A[n+i+1] = pow(xx[k],2.0)/12.0 + xx[k]*xx[k+1]/6.0 + pow(xx[k+1],2.0)/4.0;                     
        A[2*n+i] = user->ctx.gamma_s/2.0;
        A[2*n+i+1] = user->ctx.gamma_s/2.0;
        if (k > 0) {
            A[i-1] = -2.0*xx[k-1]*xx[n+k-1];
            A[n+i-1] = -pow(xx[k-1],2.0);
            A[2*n+i-1] = -user->ctx.gamma_s;
        }
        // over f_{6*n+i+1}
        j = i;
        i += nn;
        A[i] = -A[j];
        A[i+1] = -A[j+1] + 2.0*xx[k+1]*xx[n+k+1];
        A[n+i] = -A[n+j];
        A[n+i+1] = -A[n+j+1] + pow(xx[k+1],2.0);
        A[2*n+i] = -A[2*n+j];
        A[2*n+i+1] = -A[2*n+j+1] + user->ctx.gamma_s;
        i += 2;
        k += 2;
    }
    /* Conservation of Energy */
    k = 0;
    for (PetscInt i = nn*(2*n); i < nn*(3*n); i+=nn){
        // over f_{12*n+i}
        A[i] = (2.0*xx[2*n+k] + xx[2*n+k+1])*(1.0+user->ctx.gamma_s)/6.0 + 3.0/10.0*xx[n+k]*pow(xx[k],2.0) + 3.0/40.0*xx[n+k+1]*pow(xx[k],2.0) + 3.0/20.0*xx[n+k]*xx[k]*xx[k+1] + 1.0/10.0*xx[n+k+1]*xx[k]*xx[k+1] + 1.0/20.0*xx[n+k]*pow(xx[k+1],2.0) + 3.0/40.0*xx[n+k+1]*pow(xx[k+1],2.0);
        A[i+1] = (xx[2*n+k] + 2.0*xx[2*n+k+1])*(1.0+user->ctx.gamma_s)/6.0 + 3.0/40.0*xx[n+k]*pow(xx[k],2.0) + 1.0/20.0*xx[n+k+1]*pow(xx[k],2.0) + 1.0/10.0*xx[n+k]*xx[k]*xx[k+1] + 3.0/20.0*xx[n+k+1]*xx[k]*xx[k+1] + 3.0/40.0*xx[n+k]*pow(xx[k+1],2.0) + 3.0/10.0*xx[n+k+1]*pow(xx[k+1],2.0);
        A[n+i] = pow(xx[k],3.0)/10.0 + 3.0/40.0*pow(xx[k],2.0)*xx[k+1] + 1.0/20.0*xx[k]*pow(xx[k+1],2.0) + 1.0/40.0*pow(xx[k+1],3.0);
        A[n+i+1] = pow(xx[k],3.0)/40.0 + 1.0/20.0*pow(xx[k],2.0)*xx[k+1] + 3.0/40.0*xx[k]*pow(xx[k+1],2.0) + 1.0/10.0*pow(xx[k+1],3.0);
        A[2*n+i] = (2.0*xx[k] + xx[k+1])*(1.0+user->ctx.gamma_s)/6.0;
        A[2*n+i+1] = (xx[k] + 2.0*xx[k+1])*(1.0+user->ctx.gamma_s)/6.0;
        if (k > 0) {
            A[i-1] = -(0.5*3.0*pow(xx[k-1],2.0)*xx[n+k-1] + (1.0+user->ctx.gamma_s)*xx[2*n+k-1]);
            A[n+i-1] = -0.5*pow(xx[k-1],3.0);
            A[2*n+i-1] = -(1.0+user->ctx.gamma_s)*xx[k-1];
        }
        // over f_{12*n+i+1}
        j = i;
        i += nn;
        A[i] = -A[j];
        A[i+1] = -A[j+1] + (3.0*pow(xx[k+1],2.0)*xx[n+k+1])/2.0 + (1.0 + user->ctx.gamma_s) * xx[2*n+k+1];
        A[n+i] = -A[n+j];
        A[n+i+1] = -A[n+j+1] + pow(xx[k+1],3.0)/2.0;
        A[2*n+i] = -A[2*n+j];
        A[2*n+i+1] = -A[2*n+j+1] + (1.0 + user->ctx.gamma_s) * xx[k+1];
        i += 2;
        k += 2;
    }

    /* and insert into matrix B */
    ierr = MatSetValues(B, nn, idx, nn, idx, A, INSERT_VALUES); CHKERRQ(ierr);

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