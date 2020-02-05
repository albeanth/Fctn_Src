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
    // conservation of mass
    ierr = VecCreate(PETSC_COMM_WORLD, &mass_basis.i); CHKERRQ(ierr);
    ierr = VecSetSizes(mass_basis.i, PETSC_DECIDE, info.order[0]); CHKERRQ(ierr);
    ierr = VecSetFromOptions(mass_basis.i); CHKERRQ(ierr);
    ierr = VecDuplicate(mass_basis.i, &mass_basis.ii); CHKERRQ(ierr);
    ierr = VecDuplicate(mass_basis.i, &mass_basis.iii); CHKERRQ(ierr);
    ierr = VecDuplicate(mass_basis.i, &mass_basis.iv); CHKERRQ(ierr);
    ierr = VecDuplicate(mass_basis.i, &mass_basis.src); CHKERRQ(ierr);
    // conservation of momentum
    ierr = VecDuplicate(mass_basis.i, &momen_basis.i); CHKERRQ(ierr);
    ierr = VecDuplicate(mass_basis.i, &momen_basis.ii); CHKERRQ(ierr);
    ierr = VecDuplicate(mass_basis.i, &momen_basis.iii); CHKERRQ(ierr);
    ierr = VecDuplicate(mass_basis.i, &momen_basis.iv); CHKERRQ(ierr);
    ierr = VecDuplicate(mass_basis.i, &momen_basis.v); CHKERRQ(ierr);
    ierr = VecDuplicate(mass_basis.i, &momen_basis.vi); CHKERRQ(ierr);
    ierr = VecDuplicate(mass_basis.i, &momen_basis.vii); CHKERRQ(ierr);
    ierr = VecDuplicate(mass_basis.i, &momen_basis.viii); CHKERRQ(ierr);
    ierr = VecDuplicate(mass_basis.i, &momen_basis.src); CHKERRQ(ierr);
    // conservation of fluid energy
    ierr = VecDuplicate(mass_basis.i, &efluid_basis.i); CHKERRQ(ierr);
    ierr = VecDuplicate(mass_basis.i, &efluid_basis.ii); CHKERRQ(ierr);
    ierr = VecDuplicate(mass_basis.i, &efluid_basis.iii); CHKERRQ(ierr);
    ierr = VecDuplicate(mass_basis.i, &efluid_basis.iv); CHKERRQ(ierr);
    ierr = VecDuplicate(mass_basis.i, &efluid_basis.v); CHKERRQ(ierr);
    ierr = VecDuplicate(mass_basis.i, &efluid_basis.vi); CHKERRQ(ierr);
    ierr = VecDuplicate(mass_basis.i, &efluid_basis.vii); CHKERRQ(ierr);
    ierr = VecDuplicate(mass_basis.i, &efluid_basis.viii); CHKERRQ(ierr);
    ierr = VecDuplicate(mass_basis.i, &efluid_basis.ix); CHKERRQ(ierr);
    ierr = VecDuplicate(mass_basis.i, &efluid_basis.x); CHKERRQ(ierr);
    ierr = VecDuplicate(mass_basis.i, &efluid_basis.xi); CHKERRQ(ierr);
    ierr = VecDuplicate(mass_basis.i, &efluid_basis.xii); CHKERRQ(ierr);
    ierr = VecDuplicate(mass_basis.i, &efluid_basis.src); CHKERRQ(ierr);

    // ------ local FE vec initialization ------
    // conservation of mass
    ierr = VecCreate(PETSC_COMM_WORLD, &ctx.loc_mass_i); CHKERRQ(ierr);
    ierr = VecSetSizes(ctx.loc_mass_i, PETSC_DECIDE, info.order[0]); CHKERRQ(ierr);
    ierr = VecSetFromOptions(ctx.loc_mass_i); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_i, &ctx.loc_mass_ii); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_i, &ctx.loc_mass_iii); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_i, &ctx.loc_mass_iv); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_i, &ctx.loc_mass_src); CHKERRQ(ierr);
    // conservation of momentum
    ierr = VecDuplicate(ctx.loc_mass_i, &ctx.loc_momen_i); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_i, &ctx.loc_momen_ii); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_i, &ctx.loc_momen_iii); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_i, &ctx.loc_momen_iv); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_i, &ctx.loc_momen_v); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_i, &ctx.loc_momen_vi); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_i, &ctx.loc_momen_vii); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_i, &ctx.loc_momen_viii); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_i, &ctx.loc_momen_src); CHKERRQ(ierr);
    // conservation of fluid energy
    ierr = VecDuplicate(ctx.loc_mass_i, &ctx.loc_efluid_i); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_i, &ctx.loc_efluid_ii); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_i, &ctx.loc_efluid_iii); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_i, &ctx.loc_efluid_iv); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_i, &ctx.loc_efluid_v); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_i, &ctx.loc_efluid_vi); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_i, &ctx.loc_efluid_vii); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_i, &ctx.loc_efluid_viii); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_i, &ctx.loc_efluid_ix); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_i, &ctx.loc_efluid_x); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_i, &ctx.loc_efluid_xi); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_i, &ctx.loc_efluid_xii); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.loc_mass_i, &ctx.loc_efluid_src); CHKERRQ(ierr);

    // ------ Global FE vec initialization ------
    // conservation of mass
    ierr = VecCreate(PETSC_COMM_WORLD, &ctx.glo_mass_i); CHKERRQ(ierr);
    ierr = VecSetSizes(ctx.glo_mass_i, PETSC_DECIDE, info.nnodes); CHKERRQ(ierr);
    ierr = VecSetFromOptions(ctx.glo_mass_i); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_i, &ctx.glo_mass_ii); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_i, &ctx.glo_mass_iii); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_i, &ctx.glo_mass_iv); CHKERRQ(ierr);
    // conservation of momentum
    ierr = VecDuplicate(ctx.glo_mass_i, &ctx.glo_momen_i); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_i, &ctx.glo_momen_ii); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_i, &ctx.glo_momen_iii); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_i, &ctx.glo_momen_iv); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_i, &ctx.glo_momen_v); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_i, &ctx.glo_momen_vi); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_i, &ctx.glo_momen_vii); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_i, &ctx.glo_momen_viii); CHKERRQ(ierr);
    // conservation of fluid energy
    ierr = VecDuplicate(ctx.glo_mass_i, &ctx.glo_efluid_i); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_i, &ctx.glo_efluid_ii); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_i, &ctx.glo_efluid_iii); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_i, &ctx.glo_efluid_iv); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_i, &ctx.glo_efluid_v); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_i, &ctx.glo_efluid_vi); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_i, &ctx.glo_efluid_vii); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_i, &ctx.glo_efluid_viii); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_i, &ctx.glo_efluid_ix); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_i, &ctx.glo_efluid_x); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_i, &ctx.glo_efluid_xi); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_i, &ctx.glo_efluid_xii); CHKERRQ(ierr);
    // source information
    ierr = VecDuplicate(ctx.glo_mass_i, &ctx.glo_mass_src); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_i, &ctx.glo_momen_src); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_i, &ctx.glo_efluid_src); CHKERRQ(ierr);
    // upwind vectors
    ierr = VecDuplicate(ctx.glo_mass_i, &ctx.mass_upwind); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_i, &ctx.momen_upwind); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_i, &ctx.efluid_upwind); CHKERRQ(ierr);
    // solution vectors
    ierr = VecDuplicate(ctx.glo_mass_i, &velocity); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_i, &density); CHKERRQ(ierr);
    ierr = VecDuplicate(ctx.glo_mass_i, &energy); CHKERRQ(ierr);
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

    PetscScalar tmp_vel, tmp_rho, tmp_efluid;
    /* ------------------------------------------------------ */
    // // uncomment for viewing results of local integrations
    // PetscScalar aa[2], bb[2], cc[2], dd[2], ee[2], ff[2], gg[2], hh[2];
    // PetscInt tmpIdx[2] = {0, 1};
    /* ------------------------------------------------------ */
    /* Sweep over elements and solve */
    for (int elem = 0; elem < info.nels; elem ++){
        xL = info.xnod[ info.nod[elem][0] ];
        xR = info.xnod[ info.nod[elem][info.order[elem]-1] ];
        dx = (xR-xL)/2.0;
        loc_to_glob_map[0] = elem*info.order[elem];
        loc_to_glob_map[1] = loc_to_glob_map[0] + (info.order[elem] - 1);

        for (int l1 = 0; l1 < qps1d.nw; l1++){
            /* map from ref elem to real elem */
            x = xL + (1.0 + qps1d.xw[l1]) * dx;
            /* evaluate basis functions */
            ierr = EvalBasis(qps1d.xw[l1], info.order[elem]); CHKERRQ(ierr);
            /* evaluate known functions and integrate */
            src_mass = MMS_Src_Mass(x);
            src_momen = MMS_Src_Momentum(x);
            src_energy = MMS_Src_Energy(x);
            // conservation of mass
            ierr = VecAXPY(ctx.loc_mass_i, qps1d.w[l1] * dx, mass_basis.i); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_mass_ii, qps1d.w[l1] * dx, mass_basis.ii); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_mass_iii, qps1d.w[l1] * dx, mass_basis.iii); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_mass_iv, qps1d.w[l1] * dx, mass_basis.iv); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_mass_src, src_mass * qps1d.w[l1] * dx, mass_basis.src); CHKERRQ(ierr);
            // // conservation of momentum
            ierr = VecAXPY(ctx.loc_momen_i, qps1d.w[l1] * dx, momen_basis.i); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_momen_ii, qps1d.w[l1] * dx, momen_basis.ii); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_momen_iii, 2.0 * qps1d.w[l1] * dx, momen_basis.iii); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_momen_iv, 2.0 * qps1d.w[l1] * dx, momen_basis.iv); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_momen_v, qps1d.w[l1] * dx, momen_basis.v); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_momen_vi, qps1d.w[l1] * dx, momen_basis.vi); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_momen_vii, qps1d.w[l1] * dx, momen_basis.vii); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_momen_viii, qps1d.w[l1] * dx, momen_basis.vii); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_momen_src, src_momen * qps1d.w[l1] * dx, momen_basis.src); CHKERRQ(ierr);
            // // conservation of energy
            ierr = VecAXPY(ctx.loc_efluid_i, qps1d.w[l1] * dx, efluid_basis.i); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_efluid_ii, qps1d.w[l1] * dx, efluid_basis.ii); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_efluid_iii, 3.0 * qps1d.w[l1] * dx, efluid_basis.iii); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_efluid_iv, 3.0 * qps1d.w[l1] * dx, efluid_basis.iv); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_efluid_v, 3.0 * qps1d.w[l1] * dx, efluid_basis.v); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_efluid_vi, 3.0 * qps1d.w[l1] * dx, efluid_basis.vi); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_efluid_vii, qps1d.w[l1] * dx, efluid_basis.vii); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_efluid_viii, qps1d.w[l1] * dx, efluid_basis.viii); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_efluid_ix, qps1d.w[l1] * dx, efluid_basis.ix); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_efluid_x, qps1d.w[l1] * dx, efluid_basis.x); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_efluid_xi, qps1d.w[l1] * dx, efluid_basis.xi); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_efluid_xii, qps1d.w[l1] * dx, efluid_basis.xii); CHKERRQ(ierr);
            ierr = VecAXPY(ctx.loc_efluid_src, src_energy * qps1d.w[l1] * dx, efluid_basis.src); CHKERRQ(ierr);
        }
        /* ------------------------------------------------------ */
        // // uncomment for viewing results of local integrations
        // VecGetValues(ctx.loc_mass_i, 2, tmpIdx, aa);
        // VecGetValues(ctx.loc_mass_ii, 2, tmpIdx, bb);
        // VecGetValues(ctx.loc_mass_iv, 2, tmpIdx, cc);
        // PetscPrintf(PETSC_COMM_WORLD, "% .8e\t% .8e\t% .8e\n", aa[0], bb[0], cc[0]);
        // PetscPrintf(PETSC_COMM_WORLD, "% .8e\t% .8e\t% .8e\n\n", aa[1], bb[1], cc[1]);
        // VecGetValues(ctx.loc_momen_i, 2, tmpIdx, aa);
        // VecGetValues(ctx.loc_momen_ii, 2, tmpIdx, bb);
        // VecGetValues(ctx.loc_momen_iv, 2, tmpIdx, cc);
        // VecGetValues(ctx.loc_momen_vi, 2, tmpIdx, dd);
        // VecGetValues(ctx.loc_momen_vii, 2, tmpIdx, ee);
        // VecGetValues(ctx.loc_momen_viii, 2, tmpIdx, ff);
        // PetscPrintf(PETSC_COMM_WORLD, "% .8e\t% .8e\t% .8e\t% .8e\t% .8e\t% .8e\n", aa[0], bb[0], cc[0], dd[0], ee[0], ff[0]);
        // PetscPrintf(PETSC_COMM_WORLD, "% .8e\t% .8e\t% .8e\t% .8e\t% .8e\t% .8e\n\n", aa[1], bb[1], cc[1], dd[1], ee[1], ff[1]);
        // VecGetValues(ctx.loc_efluid_i, 2, tmpIdx, aa);
        // VecGetValues(ctx.loc_efluid_ii, 2, tmpIdx, bb);
        // VecGetValues(ctx.loc_efluid_iv, 2, tmpIdx, cc);
        // VecGetValues(ctx.loc_efluid_vi, 2, tmpIdx, dd);
        // VecGetValues(ctx.loc_efluid_viii, 2, tmpIdx, ee);
        // VecGetValues(ctx.loc_efluid_ix, 2, tmpIdx, ff);
        // VecGetValues(ctx.loc_efluid_x, 2, tmpIdx, gg);
        // VecGetValues(ctx.loc_efluid_xii, 2, tmpIdx, hh);
        // PetscPrintf(PETSC_COMM_WORLD, "% .8e\t% .8e\t% .8e\t% .8e\t% .8e\t% .8e\t% .8e\t% .8e\n", aa[0], bb[0], cc[0], dd[0], ee[0], ff[0], gg[0], hh[0]);
        // PetscPrintf(PETSC_COMM_WORLD, "% .8e\t% .8e\t% .8e\t% .8e\t% .8e\t% .8e\t% .8e\t% .8e\n\n", aa[1], bb[1], cc[1], dd[1], ee[1], ff[1], gg[1], hh[1]);
        // exit(-1);
        /* ------------------------------------------------------ */

        /* get upwind info */
        if (elem == 0) { // use BC info
          tmp_vel = u(info.bounds[0]);
          tmp_rho = rho(info.bounds[0]);
          tmp_efluid = efluid(info.bounds[0]);
        } 
        else { // use previous element info
          index = elem * info.order[elem] - 1;
          VecGetValues(velocity, 1, &index, &tmp_vel);
          VecGetValues(density, 1, &index, &tmp_rho);
          VecGetValues(energy, 1, &index, &tmp_efluid);
        }
        // compute upwind values for mass and momentum and store them in their global matrices
        VecSetValue(ctx.mass_upwind, elem*info.order[elem], tmp_vel * tmp_rho, INSERT_VALUES);
        VecSetValue(ctx.momen_upwind, elem*info.order[elem], tmp_rho * pow(tmp_vel, 2.0) + ctx.gamma_s * tmp_efluid, INSERT_VALUES);
        VecSetValue(ctx.efluid_upwind, elem*info.order[elem], 1.0/2.0 * tmp_rho * pow(tmp_vel, 3.0) + (1.0+ctx.gamma_s) * tmp_vel * tmp_efluid, INSERT_VALUES);

        /* Map local vectors to global vectors */
        ierr = Local2Global(elem); CHKERRQ(ierr);

        /* set vector entries equal to zero (while maintaining structure) */
        ierr = InitializeLocalRHSF(); CHKERRQ(ierr);
    }
    /* print global solution */
    PetscPrintf(PETSC_COMM_WORLD, "cm\tvelocity\tdensity\n");
    for (int index = 0; index < info.nnodes; index++){
        VecGetValues(velocity, 1, &index, &tmp_vel);
        VecGetValues(density, 1, &index, &tmp_rho);
        VecGetValues(energy, 1, &index, &tmp_efluid);
        PetscPrintf(PETSC_COMM_WORLD, "%.4e\t% .8e\t% .8e\t% .8e\t% .8e\t% .8e\t% .8e\n", info.xnod[index], tmp_vel, tmp_rho, tmp_efluid, u(info.xnod[index]), rho(info.xnod[index]), efluid(info.xnod[index]) );
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
        // conservation of mass integral evaluation
        PetscScalar i[2], ii[2], iii[2], iv[2];
        i[0] = pow(b1,2.0)*b1p;
        i[1] = pow(b1,2.0)*b2p;
        ii[0] = b1*b2*b1p;
        ii[1] = b1*b2*b2p;
        iii[0] = ii[0];
        iii[1] = ii[1];
        iv[0] = pow(b2,2.0)*b1p;
        iv[1] = pow(b2,2.0)*b2p;
        ierr = VecSetValues(mass_basis.i, 2, idx, i, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(mass_basis.ii, 2, idx, ii, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(mass_basis.iii, 2, idx, iii, INSERT_VALUES); CHKERRQ(ierr); 
        ierr = VecSetValues(mass_basis.iv, 2, idx, iv, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(mass_basis.src, 2, idx, src, INSERT_VALUES); CHKERRQ(ierr);
        // conservation of momentum integral evaluation
        PetscScalar v[2], vi[2], vii[2], viii[2];
        i[0] = pow(b1,3.0)*b1p;
        i[1] = pow(b1,3.0)*b2p;
        ii[0] = pow(b1,2.0)*b2*b1p;
        ii[1] = pow(b1,2.0)*b2*b2p;
        iii[0] = ii[0];
        iii[1] = ii[1];
        iv[0] = b1*pow(b2,2.0)*b1p;
        iv[1] = b1*pow(b2,2.0)*b2p;
        v[0] = iii[0];
        v[1] = iii[1];
        vi[0] = pow(b2,3.0)*b1p;
        vi[1] = pow(b2,3.0)*b2p;
        vii[0] = b1*b1p;
        vii[1] = b1*b2p;
        viii[0] = b2*b1p;
        viii[1] = b2*b2p;
        ierr = VecSetValues(momen_basis.i, 2, idx, i, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(momen_basis.ii, 2, idx, ii, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(momen_basis.iii, 2, idx, iii, INSERT_VALUES); CHKERRQ(ierr); 
        ierr = VecSetValues(momen_basis.iv, 2, idx, iv, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(momen_basis.v, 2, idx, v, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(momen_basis.vi, 2, idx, vi, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(momen_basis.vii, 2, idx, vii, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(momen_basis.viii, 2, idx, viii, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(momen_basis.src, 2, idx, src, INSERT_VALUES); CHKERRQ(ierr);
        // conservation of energy integral evaluation
        PetscScalar ix[2], x[2], xi[2], xii[2];
        i[0] = pow(b1,4.0)*b1p;
        i[1] = pow(b1,4.0)*b2p;
        ii[0] = pow(b1,3.0)*b2*b1p;
        ii[1] = pow(b1,3.0)*b2*b2p;
        iii[0] = ii[0];
        iii[1] = ii[1];
        iv[0] = pow(b1,2.0)*pow(b2,2.0)*b1p;
        iv[1] = pow(b1,2.0)*pow(b2,2.0)*b2p;
        v[0] = iv[0];
        v[1] = iv[1];
        vi[0] = b1*pow(b2,3.0)*b1p;
        vi[1] = b1*pow(b2,3.0)*b2p;
        vii[0] = vi[0];
        vii[1] = vi[1];
        viii[0] = pow(b2,4.0)*b1p;
        viii[1] = pow(b2,4.0)*b2p;
        ix[0] = pow(b1,2.0)*b1p;
        ix[1] = pow(b1,2.0)*b2p;
        x[0] = b1*b2*b1p;
        x[1] = b1*b2*b2p;
        xi[0] = x[0];
        xi[1] = x[1];
        xii[0] = pow(b2,2.0)*b1p;
        xii[1] = pow(b2,2.0)*b2p;
        ierr = VecSetValues(efluid_basis.i, 2, idx, i, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(efluid_basis.ii, 2, idx, ii, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(efluid_basis.iii, 2, idx, iii, INSERT_VALUES); CHKERRQ(ierr); 
        ierr = VecSetValues(efluid_basis.iv, 2, idx, iv, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(efluid_basis.v, 2, idx, v, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(efluid_basis.vi, 2, idx, vi, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(efluid_basis.vii, 2, idx, vii, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(efluid_basis.viii, 2, idx, viii, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(efluid_basis.ix, 2, idx, ix, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(efluid_basis.x, 2, idx, x, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(efluid_basis.xi, 2, idx, xi, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(efluid_basis.xii, 2, idx, xii, INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValues(efluid_basis.src, 2, idx, src, INSERT_VALUES); CHKERRQ(ierr);
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
    // conservation of mass
    ierr = VecSet(ctx.loc_mass_i, 0.0); CHKERRQ(ierr);
    ierr = VecSet(ctx.loc_mass_ii, 0.0); CHKERRQ(ierr);
    ierr = VecSet(ctx.loc_mass_iii, 0.0); CHKERRQ(ierr);
    ierr = VecSet(ctx.loc_mass_iv, 0.0); CHKERRQ(ierr);
    ierr = VecSet(ctx.loc_mass_src, 0.0); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_mass_i); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_mass_ii); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_mass_iii); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_mass_iv); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_mass_src); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_mass_i); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_mass_ii); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_mass_iii); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_mass_iv); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_mass_src); CHKERRQ(ierr);
    // conservation of momentum
    ierr = VecSet(ctx.loc_momen_i, 0.0); CHKERRQ(ierr);
    ierr = VecSet(ctx.loc_momen_ii, 0.0); CHKERRQ(ierr);
    ierr = VecSet(ctx.loc_momen_iii, 0.0); CHKERRQ(ierr);
    ierr = VecSet(ctx.loc_momen_iv, 0.0); CHKERRQ(ierr);
    ierr = VecSet(ctx.loc_momen_v, 0.0); CHKERRQ(ierr);
    ierr = VecSet(ctx.loc_momen_vi, 0.0); CHKERRQ(ierr);
    ierr = VecSet(ctx.loc_momen_vii, 0.0); CHKERRQ(ierr);
    ierr = VecSet(ctx.loc_momen_viii, 0.0); CHKERRQ(ierr);
    ierr = VecSet(ctx.loc_momen_src, 0.0); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_momen_i); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_momen_ii); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_momen_iii); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_momen_iv); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_momen_v); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_momen_vi); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_momen_vii); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_momen_viii); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_momen_src); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_momen_i); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_momen_ii); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_momen_iii); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_momen_iv); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_momen_v); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_momen_vi); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_momen_vii); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_momen_viii); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_momen_src); CHKERRQ(ierr);
    // conservation of fluid energy
    ierr = VecSet(ctx.loc_efluid_i, 0.0); CHKERRQ(ierr);
    ierr = VecSet(ctx.loc_efluid_ii, 0.0); CHKERRQ(ierr);
    ierr = VecSet(ctx.loc_efluid_iii, 0.0); CHKERRQ(ierr);
    ierr = VecSet(ctx.loc_efluid_iv, 0.0); CHKERRQ(ierr);
    ierr = VecSet(ctx.loc_efluid_v, 0.0); CHKERRQ(ierr);
    ierr = VecSet(ctx.loc_efluid_vi, 0.0); CHKERRQ(ierr);
    ierr = VecSet(ctx.loc_efluid_vii, 0.0); CHKERRQ(ierr);
    ierr = VecSet(ctx.loc_efluid_viii, 0.0); CHKERRQ(ierr);
    ierr = VecSet(ctx.loc_efluid_ix, 0.0); CHKERRQ(ierr);
    ierr = VecSet(ctx.loc_efluid_x, 0.0); CHKERRQ(ierr);
    ierr = VecSet(ctx.loc_efluid_xi, 0.0); CHKERRQ(ierr);
    ierr = VecSet(ctx.loc_efluid_xii, 0.0); CHKERRQ(ierr);
    ierr = VecSet(ctx.loc_efluid_src, 0.0); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_efluid_i); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_efluid_ii); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_efluid_iii); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_efluid_iv); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_efluid_v); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_efluid_vi); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_efluid_vii); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_efluid_viii); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_efluid_ix); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_efluid_x); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_efluid_xi); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_efluid_xii); CHKERRQ(ierr);
    ierr = VecAssemblyBegin(ctx.loc_efluid_src); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_efluid_i); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_efluid_ii); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_efluid_iii); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_efluid_iv); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_efluid_v); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_efluid_vi); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_efluid_vii); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_efluid_viii); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_efluid_ix); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_efluid_x); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_efluid_xi); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(ctx.loc_efluid_xii); CHKERRQ(ierr);
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
    // conservation of masss
    ierr = VecGetValues(ctx.loc_mass_i, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_mass_i, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecGetValues(ctx.loc_mass_ii, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_mass_ii, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecGetValues(ctx.loc_mass_iii, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_mass_iii, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecGetValues(ctx.loc_mass_iv, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_mass_iv, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecGetValues(ctx.loc_mass_src, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_mass_src, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    // conservation of momentum
    ierr = VecGetValues(ctx.loc_momen_i, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_momen_i, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecGetValues(ctx.loc_momen_ii, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_momen_ii, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecGetValues(ctx.loc_momen_iii, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_momen_iii, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecGetValues(ctx.loc_momen_iv, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_momen_iv, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecGetValues(ctx.loc_momen_v, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_momen_v, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecGetValues(ctx.loc_momen_vi, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_momen_vi, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecGetValues(ctx.loc_momen_vii, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_momen_vii, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecGetValues(ctx.loc_momen_viii, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_momen_viii, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecGetValues(ctx.loc_momen_src, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_momen_src, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    // conservation of energy
    ierr = VecGetValues(ctx.loc_efluid_i, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_efluid_i, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecGetValues(ctx.loc_efluid_ii, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_efluid_ii, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecGetValues(ctx.loc_efluid_iii, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_efluid_iii, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecGetValues(ctx.loc_efluid_iv, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_efluid_iv, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecGetValues(ctx.loc_efluid_v, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_efluid_v, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecGetValues(ctx.loc_efluid_vi, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_efluid_vi, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecGetValues(ctx.loc_efluid_vii, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_efluid_vii, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecGetValues(ctx.loc_efluid_viii, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_efluid_viii, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecGetValues(ctx.loc_efluid_ix, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_efluid_ix, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecGetValues(ctx.loc_efluid_x, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_efluid_x, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecGetValues(ctx.loc_efluid_xi, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_efluid_xi, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecGetValues(ctx.loc_efluid_xii, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_efluid_xii, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecGetValues(ctx.loc_efluid_src, 2, loc_idx, loc_val); CHKERRQ(ierr);
    ierr = VecSetValues(ctx.glo_efluid_src, 2, glo_idx, loc_val, INSERT_VALUES); CHKERRQ(ierr);
    return ierr;
}

    /*
     *  create nonlinear solver context and solve equations
     */
    /* Set initial guess */
    // ierr = VecSet(x, 1.0);CHKERRQ(ierr);
    // Mass
    ierr = VecSetValue(soln, 0, vL, INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(soln, 1, vR, INSERT_VALUES);CHKERRQ(ierr);
    // Momentum
    ierr = VecSetValue(soln, 2, rL, INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(soln, 3, rR, INSERT_VALUES);CHKERRQ(ierr);
    // Energy
    ierr = VecSetValue(soln, 4, eL, INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(soln, 5, eR, INSERT_VALUES);CHKERRQ(ierr);

    /* Solve nonlinear system */
    ierr = SNESSolve(snes, NULL, soln); CHKERRQ(ierr);
    // ierr = SNESView(snes, PETSC_VIEWER_STDOUT_WORLD);

    /* Map solution to global solution */
    PetscScalar value[6];             // array of values computed from NL solve
    PetscInt idx[6] = {0,1,2,3,4,5};  // indices of values to pull from NL solve
    PetscInt el = elem * 2;           // 2 equals num of unknowns per cell
    ierr = VecGetValues(soln, 6, idx, value); CHKERRQ(ierr);
    ierr = VecSetValue(velocity, el, value[0], INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecSetValue(velocity, el+1, value[1], INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecSetValue(density, el, value[2], INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecSetValue(density, el+1, value[3], INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecSetValue(energy, el, value[4], INSERT_VALUES); CHKERRQ(ierr);
    ierr = VecSetValue(energy, el+1, value[5], INSERT_VALUES); CHKERRQ(ierr);
    // PetscPrintf(PETSC_COMM_WORLD, "%.4e\t% .8e\t% .8e\t% .8e\n", info.xnod[el], value[0], value[2], value[4]);
    // PetscPrintf(PETSC_COMM_WORLD, "%.4e\t% .8e\t% .8e\t% .8e\n", info.xnod[el+1], value[1], value[3], value[5]);
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

    NonLinear *user = (NonLinear *)ctx;
    PetscErrorCode ierr;
    const PetscScalar *xx;
    PetscScalar       *ff;
    const int N {user->info.nnodes};

    PetscScalar *mass_i, *mass_ii, *mass_iii, *mass_iv;
    PetscScalar *mass_src;
    mass_i = (PetscScalar*) malloc(N * sizeof(PetscScalar));
    mass_ii = (PetscScalar *)malloc(N * sizeof(PetscScalar));
    mass_iii = (PetscScalar *)malloc(N * sizeof(PetscScalar));
    mass_iv = (PetscScalar *) malloc(N * sizeof(PetscScalar));
    mass_src = (PetscScalar *) malloc(N * sizeof(PetscScalar));
    PetscScalar *momen_i, *momen_ii, *momen_iii, *momen_iv;
    PetscScalar *momen_v, *momen_vi, *momen_vii, *momen_viii;
    PetscScalar *momen_src;
    momen_i = (PetscScalar *)malloc(N * sizeof(PetscScalar));
    momen_ii = (PetscScalar *)malloc(N * sizeof(PetscScalar));
    momen_iii = (PetscScalar *)malloc(N * sizeof(PetscScalar));
    momen_iv = (PetscScalar *)malloc(N * sizeof(PetscScalar));
    momen_v = (PetscScalar *)malloc(N * sizeof(PetscScalar));
    momen_vi = (PetscScalar *)malloc(N * sizeof(PetscScalar));
    momen_vii = (PetscScalar *)malloc(N * sizeof(PetscScalar));
    momen_viii = (PetscScalar *)malloc(N * sizeof(PetscScalar));
    momen_src = (PetscScalar *)malloc(N * sizeof(PetscScalar));
    PetscScalar *efluid_i, *efluid_ii, *efluid_iii, *efluid_iv;
    PetscScalar *efluid_v, *efluid_vi, *efluid_vii, *efluid_viii;
    PetscScalar *efluid_ix, *efluid_x, *efluid_xi, *efluid_xii;
    PetscScalar *efluid_src;
    efluid_i = (PetscScalar *)malloc(N * sizeof(PetscScalar));
    efluid_ii = (PetscScalar *)malloc(N * sizeof(PetscScalar));
    efluid_iii = (PetscScalar *)malloc(N * sizeof(PetscScalar));
    efluid_iv = (PetscScalar *)malloc(N * sizeof(PetscScalar));
    efluid_v = (PetscScalar *)malloc(N * sizeof(PetscScalar));
    efluid_vi = (PetscScalar *)malloc(N * sizeof(PetscScalar));
    efluid_vii = (PetscScalar *)malloc(N * sizeof(PetscScalar));
    efluid_viii = (PetscScalar *)malloc(N * sizeof(PetscScalar));
    efluid_ix = (PetscScalar *)malloc(N * sizeof(PetscScalar));
    efluid_x = (PetscScalar *)malloc(N * sizeof(PetscScalar));
    efluid_xi = (PetscScalar *)malloc(N * sizeof(PetscScalar));
    efluid_xii = (PetscScalar *)malloc(N * sizeof(PetscScalar));
    efluid_src = (PetscScalar *)malloc(N * sizeof(PetscScalar));

    // assign petsc Vec to c array for use
    PetscInt *idx;
    idx = (PetscInt*) malloc(N * sizeof(PetscInt));
    for (int i=0; i<N; i++){
        idx[i] = i;
    }
    // conservation of mass
    ierr = VecGetValues(user->ctx.glo_mass_i, N, idx, mass_i); CHKERRQ(ierr);
    ierr = VecGetValues(user->ctx.glo_mass_ii, N, idx, mass_ii); CHKERRQ(ierr);
    ierr = VecGetValues(user->ctx.glo_mass_iii, N, idx, mass_iii); CHKERRQ(ierr);
    ierr = VecGetValues(user->ctx.glo_mass_iv, N, idx, mass_iv); CHKERRQ(ierr);
    // conservation of momentum
    ierr = VecGetValues(user->ctx.glo_momen_i, N, idx, momen_i); CHKERRQ(ierr);
    ierr = VecGetValues(user->ctx.glo_momen_ii, N, idx, momen_ii); CHKERRQ(ierr);
    ierr = VecGetValues(user->ctx.glo_momen_iii, N, idx, momen_iii); CHKERRQ(ierr);
    ierr = VecGetValues(user->ctx.glo_momen_iv, N, idx, momen_iv); CHKERRQ(ierr);
    ierr = VecGetValues(user->ctx.glo_momen_v, N, idx, momen_v); CHKERRQ(ierr);
    ierr = VecGetValues(user->ctx.glo_momen_vi, N, idx, momen_vi); CHKERRQ(ierr);
    ierr = VecGetValues(user->ctx.glo_momen_vii, N, idx, momen_vii); CHKERRQ(ierr);
    ierr = VecGetValues(user->ctx.glo_momen_viii, N, idx, momen_viii); CHKERRQ(ierr);
    // conservation of fluid energy
    ierr = VecGetValues(user->ctx.glo_efluid_i, N, idx, efluid_i); CHKERRQ(ierr);
    ierr = VecGetValues(user->ctx.glo_efluid_ii, N, idx, efluid_ii); CHKERRQ(ierr);
    ierr = VecGetValues(user->ctx.glo_efluid_iii, N, idx, efluid_iii); CHKERRQ(ierr);
    ierr = VecGetValues(user->ctx.glo_efluid_iv, N, idx, efluid_iv); CHKERRQ(ierr);
    ierr = VecGetValues(user->ctx.glo_efluid_v, N, idx, efluid_v); CHKERRQ(ierr);
    ierr = VecGetValues(user->ctx.glo_efluid_vi, N, idx, efluid_vi); CHKERRQ(ierr);
    ierr = VecGetValues(user->ctx.glo_efluid_vii, N, idx, efluid_vii); CHKERRQ(ierr);
    ierr = VecGetValues(user->ctx.glo_efluid_viii, N, idx, efluid_viii); CHKERRQ(ierr);
    ierr = VecGetValues(user->ctx.glo_efluid_ix, N, idx, efluid_ix); CHKERRQ(ierr);
    ierr = VecGetValues(user->ctx.glo_efluid_x, N, idx, efluid_x); CHKERRQ(ierr);
    ierr = VecGetValues(user->ctx.glo_efluid_xi, N, idx, efluid_xi); CHKERRQ(ierr);
    ierr = VecGetValues(user->ctx.glo_efluid_xii, N, idx, efluid_xii); CHKERRQ(ierr);

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
    // Conservation of mass
    
    // Conservation of Momentum
    
    // Conservation of Energy
    
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