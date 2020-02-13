#include "Nonlinear.hpp"

PetscErrorCode NonLinear::Initialize_NL_1D(){
    /*
     *  initialize all PETSC data structures for NonLinear::NL_1D
     */
    if (info.nnodes == NAN){
        SETERRQ(PETSC_COMM_WORLD,1,"Please build mesh via SetUpGrid::add_DFEMGrid before executing Initialize_NL_1D()!");
    }
    N = info.nnodes * 4;
    /* ------ Initialize SNES ------ */
    ierr = SNESCreate(PETSC_COMM_WORLD, &snes); CHKERRQ(ierr); // create context 
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
    ierr = SNESSetJacobian(snes, J, J, NULL, NULL); CHKERRQ(ierr); // set Jacobian matrix data strcuture and evaluation routine

    // /*
    //  * nonlinear preconditioning
    //  * - options are passed to this via the "-npc_" prefix
    //  */
    // SNES psnes;                         // "inner" snes instance
    // SNESGetNPC(snes, &psnes);           // extract inner snes
    // // SNESSetType(psnes, SNESSHELL);
    // // SNESShellSetSolve(psnes, );

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
    ierr = VecDuplicate(mass_basis_src, &efluid_basis_i); CHKERRQ(ierr);
    ierr = VecDuplicate(mass_basis_src, &efluid_basis_ii); CHKERRQ(ierr);
    ierr = VecDuplicate(mass_basis_src, &erad_basis_src); CHKERRQ(ierr);
    ierr = VecDuplicate(mass_basis_src, &erad_basis_i); CHKERRQ(ierr);
    ierr = VecDuplicate(mass_basis_src, &erad_basis_ii); CHKERRQ(ierr);
    ierr = VecDuplicate(mass_basis_src, &erad_basis_iii); CHKERRQ(ierr);
    ierr = VecDuplicate(mass_basis_src, &erad_basis_iv); CHKERRQ(ierr);
    
    ierr = VecCreate(PETSC_COMM_WORLD, &ctx.loc_mass_src); CHKERRQ(ierr);
    ierr = VecSetSizes(ctx.loc_mass_src, PETSC_DECIDE, info.order[0]); CHKERRQ(ierr);
    ierr = VecSetFromOptions(ctx.loc_mass_src); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_src, &ctx.loc_momen_src); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_src, &ctx.loc_efluid_src); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_src, &ctx.loc_efluid_i); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_src, &ctx.loc_efluid_ii); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_src, &ctx.loc_erad_src); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_src, &ctx.loc_erad_i); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_src, &ctx.loc_erad_ii); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_src, &ctx.loc_erad_iii); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_src, &ctx.loc_erad_iv); CHKERRQ(ierr);

    // ------ Global FE vec initialization ------
    // source information
    ierr = VecCreate(PETSC_COMM_WORLD, &ctx.glo_mass_src); CHKERRQ(ierr);
    ierr = VecSetSizes(ctx.glo_mass_src, PETSC_DECIDE, info.nnodes); CHKERRQ(ierr);
    ierr = VecSetFromOptions(ctx.glo_mass_src); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_src, &ctx.glo_momen_src); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_src, &ctx.glo_efluid_src); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_src, &ctx.glo_efluid_i); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_src, &ctx.glo_efluid_ii); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_src, &ctx.glo_erad_src); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_src, &ctx.glo_erad_i); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_src, &ctx.glo_erad_ii); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_src, &ctx.glo_erad_iii); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_src, &ctx.glo_erad_iv); CHKERRQ(ierr);
    // solution vectors
    ierr = VecDuplicate(ctx.glo_mass_src, &velocity); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_src, &density); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_src, &FluidEnergy); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_src, &RadEnergy); CHKERRQ(ierr);
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
    ierr = VecSetValue(FluidEnergy, 0, efluid(info.bounds[0]), INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecSetValue(RadEnergy, 0, erad(info.bounds[0]), INSERT_VALUES); CHKERRQ(ierr);
    // assemble vectors
    ierr = VecAssemblyBegin(velocity); CHKERRQ(ierr); 
    ierr = VecAssemblyEnd(velocity); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(density); CHKERRQ(ierr); 
    ierr = VecAssemblyEnd(density); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(FluidEnergy); CHKERRQ(ierr); 
    ierr = VecAssemblyEnd(FluidEnergy); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(RadEnergy); CHKERRQ(ierr); 
    ierr = VecAssemblyEnd(RadEnergy); CHKERRQ(ierr);

    ierr = InitializeLocalRHSF(); CHKERRQ(ierr);
    get1D_QPs(info.maxord, qps1d); // set qps

    /* Integrate over each cell locally and store to global vector */
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
            src_efluid = MMS_Src_EFluid(x);
            src_erad = MMS_Src_ERad(x);

            ierr = VecAXPY(ctx.loc_mass_src, src_mass * qps1d.w[l1] * dx, mass_basis_src); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_momen_src, src_momen * qps1d.w[l1] * dx, momen_basis_src); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_efluid_src, src_efluid * qps1d.w[l1] * dx, efluid_basis_src); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_efluid_i,  c/(3.0*sig_r(x)) * qps1d.w[l1] * dx / pow(dx,2.0), efluid_basis_i); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_efluid_ii, c/(3.0*sig_r(x)) * qps1d.w[l1] * dx / pow(dx,2.0), efluid_basis_ii); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_erad_src, src_erad * qps1d.w[l1] * dx / pow(dx,2.0), erad_basis_src); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_erad_i,  c/(3.0*sig_r(x)) * qps1d.w[l1] * dx / pow(dx,2.0), erad_basis_i); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_erad_ii, c/(3.0*sig_r(x)) * qps1d.w[l1] * dx / pow(dx,2.0), erad_basis_ii); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_erad_iii, c*sig_r(x) * qps1d.w[l1] * dx, erad_basis_iii); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_erad_iv,  c*sig_r(x) * qps1d.w[l1] * dx, erad_basis_iv); CHKERRQ(ierr);
        }

        /* Map local vectors to global vectors */
        ierr = Local2Global(elem); CHKERRQ(ierr);

        /* set vector entries equal to zero (while maintaining structure) */
        ierr = InitializeLocalRHSF(); CHKERRQ(ierr);
    }
    /* Solve nonlinear system of equations over elem */
    ierr = NLSolve(); CHKERRQ(ierr);
    /* print global solution */
    PetscScalar tmp_vel, tmp_rho, tmp_efluid, tmp_erad;
    PetscPrintf(PETSC_COMM_WORLD, "cm\tvelocity\tdensity\n");
    for (int index = 0; index < info.nnodes; index++){
        VecGetValues(velocity, 1, &index, &tmp_vel);
        VecGetValues(density, 1, &index, &tmp_rho);
        VecGetValues(FluidEnergy, 1, &index, &tmp_efluid);
        VecGetValues(RadEnergy, 1, &index, &tmp_erad);
        PetscPrintf(PETSC_COMM_WORLD, "%.4e\t% .8e\t% .8e\t% .8e\t% .8e\t% .8e\t% .8e\t% .8e\t% .8e\n", info.xnod[index], tmp_vel, tmp_rho, tmp_efluid, tmp_erad, u(info.xnod[index]), rho(info.xnod[index]), efluid(info.xnod[index]), erad(info.xnod[index]) );
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
        PetscScalar i[]   = {b1p*b1p, b1p*b2p};
        PetscScalar ii[]  = {b2p*b1p, b2p*b2p};
        PetscScalar iii[] = {b1*b1, b1*b2};
        PetscScalar iv[]  = {b2*b1, b2*b2};
        ierr = VecSetValues(mass_basis_src, 2, idx, src, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(momen_basis_src, 2, idx, src, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(efluid_basis_src, 2, idx, src, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(erad_basis_src, 2, idx, src, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(efluid_basis_i, 2, idx, i, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(efluid_basis_ii, 2, idx, ii, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(erad_basis_i, 2, idx, i, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(erad_basis_ii, 2, idx, ii, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(erad_basis_iii, 2, idx, iii, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(erad_basis_iv, 2, idx, iv, INSERT_VALUES); CHKERRQ(ierr);
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
                ierr = VecGetValues(FluidEnergy, 1, &mynum, &tmp_em); CHKERRQ(ierr);
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
    ierr = VecSet(ctx.loc_efluid_i, 0.0); CHKERRQ(ierr);
    ierr = VecSet(ctx.loc_efluid_ii, 0.0); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_efluid_src); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_efluid_i); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_efluid_ii); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_efluid_src); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_efluid_i); CHKERRQ(ierr);    
    ierr = VecAssemblyEnd(ctx.loc_efluid_ii); CHKERRQ(ierr);    
    // radiation diffusion
    ierr = VecSet(ctx.loc_erad_src, 0.0); CHKERRQ(ierr);
    ierr = VecSet(ctx.loc_erad_i, 0.0); CHKERRQ(ierr);
    ierr = VecSet(ctx.loc_erad_ii, 0.0); CHKERRQ(ierr);
    ierr = VecSet(ctx.loc_erad_iii, 0.0); CHKERRQ(ierr);
    ierr = VecSet(ctx.loc_erad_iv, 0.0); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_erad_src); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_erad_i); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_erad_ii); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_erad_iii); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_erad_iv); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_erad_src); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_erad_i); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_erad_ii); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_erad_iii); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_erad_iv); CHKERRQ(ierr);

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
    ierr = VecGetValues(ctx.loc_erad_src, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_erad_src, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecGetValues(ctx.loc_efluid_i, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_efluid_i, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecGetValues(ctx.loc_efluid_ii, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_efluid_ii, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecGetValues(ctx.loc_erad_i, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_erad_i, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecGetValues(ctx.loc_erad_ii, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_erad_ii, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecGetValues(ctx.loc_erad_iii, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_erad_iii, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecGetValues(ctx.loc_erad_iv, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_erad_iv, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);

    return ierr;
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
        ierr = VecSetValue(soln, i+3*info.nnodes, erad(0.0), INSERT_VALUES); CHKERRQ(ierr);
    }
    /* Solve nonlinear system */
    ierr = SNESSolve(snes, NULL, soln); CHKERRQ(ierr);
    // ierr = SNESView(snes, PETSC_VIEWER_STDOUT_WORLD);
    
    /* Map computed solution to solution vectors */
    PetscScalar *tmp_vel, *tmp_rho, *tmp_efluid, *tmp_erad;
    tmp_vel = (PetscScalar*) malloc(info.nnodes * sizeof(PetscScalar));
    tmp_rho = (PetscScalar*) malloc(info.nnodes * sizeof(PetscScalar));
    tmp_efluid = (PetscScalar*) malloc(info.nnodes * sizeof(PetscScalar));
    tmp_erad = (PetscScalar*) malloc(info.nnodes * sizeof(PetscScalar));

    PetscInt *idu, *idr, *idem, *iderad;
    idu = (PetscInt*) malloc(info.nnodes * sizeof(PetscInt));
    idr = (PetscInt *)malloc(info.nnodes * sizeof(PetscInt));
    idem = (PetscInt *)malloc(info.nnodes * sizeof(PetscInt));
    iderad = (PetscInt *)malloc(info.nnodes * sizeof(PetscInt));
    for (int i=0; i<info.nnodes; i++){
        idu[i] = i;
        idr[i] = i + info.nnodes;
        idem[i] = i + 2*info.nnodes;
        iderad[i] = i + 3*info.nnodes;
    }
    ierr = VecGetValues(soln, info.nnodes, idu, tmp_vel); CHKERRQ(ierr);
    ierr = VecGetValues(soln, info.nnodes, idr, tmp_rho); CHKERRQ(ierr);
    ierr = VecGetValues(soln, info.nnodes, idem, tmp_efluid); CHKERRQ(ierr);
    ierr = VecGetValues(soln, info.nnodes, iderad, tmp_erad); CHKERRQ(ierr);
    ierr = VecSetValues(velocity, info.nnodes, idu, tmp_vel, INSERT_VALUES);
    ierr = VecSetValues(density, info.nnodes, idu, tmp_rho, INSERT_VALUES);
    ierr = VecSetValues(FluidEnergy, info.nnodes, idu, tmp_efluid, INSERT_VALUES);
    ierr = VecSetValues(RadEnergy, info.nnodes, idu, tmp_erad, INSERT_VALUES);

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
    PetscScalar *efluid_src;
    PetscScalar *efluid_i, *efluid_ii;
    efluid_src = (PetscScalar *)malloc(nn * sizeof(PetscScalar));
    efluid_i = (PetscScalar *)malloc(nn * sizeof(PetscScalar));
    efluid_ii = (PetscScalar *)malloc(nn * sizeof(PetscScalar));
    PetscScalar *erad_src;
    PetscScalar *erad_i, *erad_ii, *erad_iii, *erad_iv;
    erad_src = (PetscScalar *)malloc(nn * sizeof(PetscScalar));
    erad_i = (PetscScalar *)malloc(nn * sizeof(PetscScalar));
    erad_ii = (PetscScalar *)malloc(nn * sizeof(PetscScalar));
    erad_iii = (PetscScalar *)malloc(nn * sizeof(PetscScalar));
    erad_iv = (PetscScalar *)malloc(nn * sizeof(PetscScalar));

    // assign petsc Vec to c array for use
    PetscInt *idx;
    idx = (PetscInt*) malloc(nn * sizeof(PetscInt));
    for (int i=0; i<nn; i++){
        idx[i] = i;
    }
    ierr = VecGetValues(user->ctx.glo_mass_src, nn, idx, mass_src); CHKERRQ(ierr);
    ierr = VecGetValues(user->ctx.glo_momen_src, nn, idx, momen_src); CHKERRQ(ierr);
    ierr = VecGetValues(user->ctx.glo_efluid_src, nn, idx, efluid_src); CHKERRQ(ierr);
    ierr = VecGetValues(user->ctx.glo_efluid_i,   nn, idx, efluid_i); CHKERRQ(ierr);
    ierr = VecGetValues(user->ctx.glo_efluid_ii,  nn, idx, efluid_ii); CHKERRQ(ierr);
    ierr = VecGetValues(user->ctx.glo_erad_src, nn, idx, erad_src); CHKERRQ(ierr);
    ierr = VecGetValues(user->ctx.glo_erad_i,   nn, idx, erad_i); CHKERRQ(ierr);
    ierr = VecGetValues(user->ctx.glo_erad_ii,  nn, idx, erad_ii); CHKERRQ(ierr);
    ierr = VecGetValues(user->ctx.glo_erad_iii, nn, idx, erad_iii); CHKERRQ(ierr);
    ierr = VecGetValues(user->ctx.glo_erad_iv,  nn, idx, erad_iv); CHKERRQ(ierr);
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
    const double r{user->info.bounds[1]};
    PetscScalar efluid_upwind[2], erad_upwind[2];
    // PetscPrintf(PETSC_COMM_SELF, "In FormFunction\n");
    for (PetscInt i=0; i<nn; i+=user->info.order[0]){
        if (i == 0){
            // PetscPrintf(PETSC_COMM_SELF, "  \nel = %d, (%d, %d), ( i==%d )\n", i/2, i, i+1, i);
            /* fluid mass and momentum upwind values for left cell edge */
            mass_upwind = user->rho(l) * user->u(l);
            momen_upwind = user->rho(l) * pow(user->u(l),2.0) + user->ctx.gamma_s * user->efluid(l) + user->erad(l)/3.0;
            /* fluid energy upwind values for left and right edge */
            efluid_upwind[0] = -0.5 * user->rho(l) * pow(user->u(l),3.0) - (1.0 + user->ctx.gamma_s) * user->u(l) * user->efluid(l) + 0.0; // + 0.0 -> net current at left edge is zero for reflecting BCs
            efluid_upwind[1] =  0.5 * xx[nn+i+1] * pow(xx[i+1],3.0) + (1.0+user->ctx.gamma_s)*xx[i+1]*xx[2*nn+i+1]
                                + ( xx[3*nn+i+1]/4.0 - user->c/(6.0*user->sig_r(user->info.xnod[i+1])) * (xx[3*nn+i] - xx[3*nn+i+1])/(user->info.xnod[i+1] - user->info.xnod[i]))
                                - ( xx[3*nn+i+2]/4.0 + user->c/(6.0*user->sig_r(user->info.xnod[i+2])) * (xx[3*nn+i+2] - xx[3*nn+i+3])/(user->info.xnod[i+3] - user->info.xnod[i+2])) ;
            /* radiation energy upwind values for left and right edge */
            erad_upwind[0] = 0.0; // net current at left edge is zero for reflecting BCs
            erad_upwind[1] =   ( xx[3*nn+i+1]/4.0 - user->c/(6.0*user->sig_r(user->info.xnod[i+1])) * (xx[3*nn+i] - xx[3*nn+i+1])/(user->info.xnod[i+1] - user->info.xnod[i]))
                             - ( xx[3*nn+i+2]/4.0 + user->c/(6.0*user->sig_r(user->info.xnod[i+2])) * (xx[3*nn+i+2] - xx[3*nn+i+3])/(user->info.xnod[i+3] - user->info.xnod[i+2])) ;
        }
        else if (i == nn-2){
            // PetscPrintf(PETSC_COMM_SELF, "  \nel = %d, (%d, %d), (i == %d)\n", i/2, i, i+1, i);
            /* fluid mass and momentum upwind values for left cell edge */
            mass_upwind = xx[nn+i-1] * xx[i-1];
            momen_upwind = xx[nn+i-1] * pow(xx[i-1],2.0) + user->ctx.gamma_s * xx[2*nn+i-1] + xx[3*nn+i-1]/3.0;
            /* fluid energy upwind values for left and right */
            efluid_upwind[0] = -0.5 * xx[nn+i-1] * pow(xx[i-1],3.0) - (1.0 + user->ctx.gamma_s) * xx[i-1] * xx[2*nn+i-1]
                                - ( xx[3*nn+i-1]/4.0 - user->c/(6.0*user->sig_r(user->info.xnod[i-1])) * (xx[3*nn+i-2] - xx[3*nn+i-1])/(user->info.xnod[i-1] - user->info.xnod[i-2]) )
                                + ( xx[3*nn+i]/4.0   + user->c/(6.0*user->sig_r(user->info.xnod[i]))   * (xx[3*nn+i] - xx[3*nn+i+1])/(user->info.xnod[i+1] - user->info.xnod[i]) );
            efluid_upwind[1] =  0.5 * user->rho(r) * pow(user->u(r),3.0) + (1.0+user->ctx.gamma_s)*user->u(r)*user->efluid(r) + 0.0; // + 0.0 -> net current at right edge is zero for reflecting BCs
            /* radiation energy upwind values for left and right edge */
            erad_upwind[0] = - ( xx[3*nn+i-1]/4.0 - user->c/(6.0*user->sig_r(user->info.xnod[i-1])) * (xx[3*nn+i-2] - xx[3*nn+i-1])/(user->info.xnod[i-1] - user->info.xnod[i-2]) )
                             + ( xx[3*nn+i]/4.0   + user->c/(6.0*user->sig_r(user->info.xnod[i]))   * (xx[3*nn+i] - xx[3*nn+i+1])/(user->info.xnod[i+1] - user->info.xnod[i]) );
            erad_upwind[1] = 0.0; // net current at left edge is zero for reflecting BCs
        }
        else{
            // PetscPrintf(PETSC_COMM_SELF, "  \nel = %d, (%d, %d)\n", i/2, i, i+1);
            /* fluid mass and momentum upwind values for left cell edge */
            mass_upwind = xx[nn+i-1] * xx[i-1];
            momen_upwind = xx[nn+i-1] * pow(xx[i-1],2.0) + user->ctx.gamma_s * xx[2*nn+i-1] + xx[3*nn+i-1]/3.0;
            /* fluid energy upwind values for left and right edge */
            efluid_upwind[0] = -0.5 * xx[nn+i-1] * pow(xx[i-1],3.0) - (1.0 + user->ctx.gamma_s) * xx[i-1] * xx[2*nn+i-1]
                                - ( xx[3*nn+i-1]/4.0 - user->c/(6.0*user->sig_r(user->info.xnod[i-1])) * (xx[3*nn+i-2] - xx[3*nn+i-1])/(user->info.xnod[i-1] - user->info.xnod[i-2]) ) 
                                + ( xx[3*nn+i]/4.0   + user->c/(6.0*user->sig_r(user->info.xnod[i]))   * (xx[3*nn+i] - xx[3*nn+i+1])/(user->info.xnod[i+1] - user->info.xnod[i]) );
            efluid_upwind[1] =  0.5 * xx[nn+i+1] * pow(xx[i+1],3.0) + (1.0+user->ctx.gamma_s)*xx[i+1]*xx[2*nn+i+1]
                                + (xx[3*nn+i+1]/4.0 - user->c/(6.0*user->sig_r(user->info.xnod[i+1])) * (xx[3*nn+i]-xx[3*nn+i+1])/(user->info.xnod[i+1] - user->info.xnod[i]) )
                                - (xx[3*nn+i+2]/4.0 + user->c/(6.0*user->sig_r(user->info.xnod[i+2])) * (xx[3*nn+i+2]-xx[3*nn+i+3])/(user->info.xnod[i+3]-user->info.xnod[i+2]) );
            /* radiation energy upwind values for left and right edge */
            erad_upwind[0] = - ( xx[3*nn+i-1]/4.0 - user->c/(6.0*user->sig_r(user->info.xnod[i-1])) * (xx[3*nn+i-2] - xx[3*nn+i-1])/(user->info.xnod[i-1] - user->info.xnod[i-2]) )
                             + ( xx[3*nn+i]/4.0   + user->c/(6.0*user->sig_r(user->info.xnod[i]))   * (xx[3*nn+i] - xx[3*nn+i+1])/(user->info.xnod[i+1] - user->info.xnod[i]) );
            erad_upwind[1] =   ( xx[3*nn+i+1]/4.0 - user->c/(6.0*user->sig_r(user->info.xnod[i+1])) * (xx[3*nn+i] - xx[3*nn+i+1])/(user->info.xnod[i+1] - user->info.xnod[i]))
                             - ( xx[3*nn+i+2]/4.0 + user->c/(6.0*user->sig_r(user->info.xnod[i+2])) * (xx[3*nn+i+2] - xx[3*nn+i+3])/(user->info.xnod[i+3] - user->info.xnod[i+2])) ;
        }
        // PetscPrintf(PETSC_COMM_SELF, "  mass_upwind = % .8e\n", mass_upwind);
        // PetscPrintf(PETSC_COMM_SELF, "  momen_upwind = % .8e\n", mass_upwind);
        // PetscPrintf(PETSC_COMM_SELF, "  efluid_upwind = (% .8e, % .8e)\n", efluid_upwind[0], efluid_upwind[1]);
        // PetscPrintf(PETSC_COMM_SELF, "  erad_upwind = (% .8e, % .8e)\n", erad_upwind[0], erad_upwind[1]);
        // PetscPrintf(PETSC_COMM_SELF, "  mass_src = (% .8e, % .8e)\n", mass_src[0], mass_src[1]);
        // PetscPrintf(PETSC_COMM_SELF, "  momen_src = (% .8e, % .8e)\n", momen_src[0], momen_src[1]);
        // PetscPrintf(PETSC_COMM_SELF, "  efluid_src = (% .8e, % .8e)\n", efluid_src[0], efluid_src[1]);
        // PetscPrintf(PETSC_COMM_SELF, "  erad_src = (% .8e, % .8e)\n", erad_src[0], erad_src[1]);
        // PetscPrintf(PETSC_COMM_SELF, "  efluid_i = (% .8e, % .8e)\n", efluid_i[0], efluid_i[1]);
        // PetscPrintf(PETSC_COMM_SELF, "  efluid_ii = (% .8e, % .8e)\n", efluid_ii[0], efluid_ii[1]);
        // PetscPrintf(PETSC_COMM_SELF, "  erad_i = (% .8e, % .8e)\n", erad_i[0], erad_i[1]);
        // PetscPrintf(PETSC_COMM_SELF, "  erad_ii = (% .8e, % .8e)\n", erad_ii[0], erad_ii[1]);
        // PetscPrintf(PETSC_COMM_SELF, "  erad_iii = (% .8e, % .8e)\n", erad_iii[0], erad_iii[1]);
        // PetscPrintf(PETSC_COMM_SELF, "  erad_iv = (% .8e, % .8e)\n", erad_iv[0], erad_iv[1]);
        // exit(-1);
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
                    + xx[3*nn+i]                  * (1.0/6.0)
                    + xx[3*nn+i+1]                * (1.0/6.0)
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
                    - xx[3*nn+i]                  * (1.0/6.0)
                    - xx[3*nn+i+1]                * (1.0/6.0)
                    + xx[nn+i+1]*pow(xx[i+1],2.0) + (user->ctx.gamma_s*xx[2*nn+i+1]) + xx[3*nn+i+1]/3.0
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
                      + xx[3*nn+i]*efluid_i[i] + xx[3*nn+i+1]*efluid_i[i+1]
                      + efluid_upwind[0]
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
                      + xx[3*nn+i]*efluid_ii[i] + xx[3*nn+i+1]*efluid_ii[i+1]
                      + efluid_upwind[1]
                      - efluid_src[i+1];
        // radiation energy
        ff[3*nn+i]   =   xx[3*nn+i] * erad_i[i] + xx[3*nn+i+1] * erad_i[i+1]
                       + (xx[3*nn+i] - xx[2*nn+i]) * erad_iii[0] + (xx[3*nn+i+1] - xx[2*nn+i+1]) * erad_iii[1]
                       + erad_upwind[0]
                       - erad_src[i];
        ff[3*nn+i+1] =   xx[3*nn+i] * erad_ii[i] + xx[3*nn+i+1] * erad_ii[i+1]
                       + (xx[3*nn+i] - xx[2*nn+i]) * erad_iv[0] + (xx[3*nn+i+1] - xx[2*nn+i+1]) * erad_iv[1]
                       + erad_upwind[1]
                       - erad_src[i+1];
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
    ApplicationCTX *user = (ApplicationCTX *)ctx;
    const PetscScalar *xx;
    PetscScalar A[36];
    PetscErrorCode ierr;
    PetscInt idx[36] = { 0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
                        10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
                        20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                        30, 31, 32, 33, 34, 35};

    /* Get pointer to vector data */
    ierr = VecGetArrayRead(x,&xx);CHKERRQ(ierr);

    /* Compute Jacobian entries */
    /* Conservation of Mass */
    // over f0
    A[0] = xx[2]/3.0 + xx[3]/6.0;
    A[1] = xx[2]/6.0 + xx[3]/3.0;
    A[2] = xx[0]/3.0 + xx[1]/6.0;
    A[3] = xx[0]/6.0 + xx[1]/3.0;
    A[4] = 0.0;
    A[5] = 0.0;
    // over f1
    A[6] = -A[0];        
    A[7] = -A[1] + xx[3];
    A[8] = -A[2];        
    A[9] = -A[3] + xx[1];
    A[10] = 0.0;
    A[11] = 0.0;
    /* Conservation of Momentum */
    // over f2
    A[12] = 1.0/6.0 * ( xx[3]*(xx[0]+xx[1]) + xx[2]*(3.0*xx[0]+xx[1]) );
    A[13] = 1.0/6.0 * ( xx[2]*(xx[0]+xx[1]) + xx[3]*(xx[0]+3.0*xx[1]) );
    A[14] = 1/12.0 * ( 3.0*pow(xx[0],2.0) + 2.0*xx[0]*xx[1] + pow(xx[1],2.0) );
    A[15] = 1/12.0 * ( pow(xx[0],2.0) + 2.0*xx[0]*xx[1] + 3.0*pow(xx[1],2.0) );
    A[16] = user->gamma_s / 2.0;
    A[17] = user->gamma_s / 2.0;
    // over f3
    A[18] = -A[12];
    A[19] = -A[13] + 2.0*xx[1]*xx[3];
    A[20] = -A[14];
    A[21] = -A[15] + pow(xx[1], 2.0);
    A[22] = -A[16];
    A[23] = -A[17] + user->gamma_s;
    /* Conservation of Energy */
    // over f4
    A[24] = 1.0/120.0 * (40.0*xx[4]*(1.0+user->gamma_s) + 20.0*xx[5]*(1.0+user->gamma_s) + 6.0*xx[2]*(6.0*pow(xx[0],2.0) + 3.0*xx[0]*xx[1] + pow(xx[1],2.0)) + 3.0*xx[3]*(3.0*pow(xx[0],2.0) + 4.0*xx[0]*xx[1] + 3.0*pow(xx[1],2.0)) );
    A[25] = 1.0/120.0 * (20.0*xx[4]*(1.0+user->gamma_s) + 40.0*xx[5]*(1.0+user->gamma_s) + 3.0*xx[2]*(3.0*pow(xx[0],2.0) + 4.0*xx[0]*xx[1] + 3.0*pow(xx[1],2.0)) + 6.0*xx[3]*(pow(xx[0],2.0) + 3.0*xx[0]*xx[1] + 6.0*pow(xx[1],2.0)) );
    A[26] = 1.0/40.0 * (4.0*pow(xx[0],3.0) + 3.0*pow(xx[0],2.0)*xx[1] + 2.0*xx[0]*pow(xx[1],2.0) + pow(xx[1],3.0));
    A[27] = 1.0/40.0 * (pow(xx[0],3.0) + 2.0*pow(xx[0],2.0)*xx[1] + 3.0*xx[0]*pow(xx[1],2.0) + 4.0*pow(xx[1],3.0));
    A[28] = 1.0/6.0 * (1.0+user->gamma_s) * (2.0*xx[0] + xx[1]);
    A[29] = 1.0/6.0 * (1.0+user->gamma_s) * (xx[0] + 2.0*xx[1]);
    // over f5
    A[30] = -A[24];
    A[31] = -A[25] + (3.0*pow(xx[1],2.0)*xx[3])/2.0 + (1.0 + user->gamma_s) * xx[5];
    A[32] = -A[26];
    A[33] = -A[27] + pow(xx[1],3.0)/2.0;
    A[34] = -A[28];
    A[35] = -A[29] + (1.0 + user->gamma_s) * xx[1];

    /* and insert into matrix B */
    ierr = MatSetValues(B, 6, idx, 6, idx, A, INSERT_VALUES); CHKERRQ(ierr);

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