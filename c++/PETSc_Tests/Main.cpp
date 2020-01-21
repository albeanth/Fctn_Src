#include "BVP.hpp"
#include "Nonlinear.hpp"
#include <iostream>
#include <vector>

int main(int argc,char **args){
    const std::vector<int> NumOfElem {32};
    /*
     * test #'s < 100
     *   CFEM using petsc ksp libary tests
     * test #'s \geq 100
     *   nonlinear solver using petsc snes library tests
     */
    const int selection {101};
    const bool hetero {false};
    const std::vector<double> Bnds {0.0, 1.0};

    /* vectors to store errors */
    // std::vector<double> h;
    // std::vector<double> L2Error_Vel;
    // std::vector<double> L2Error_Rho;
    // std::vector<double> H1Error_Vel;
    // std::vector<double> H1Error_Rho;

    if (selection < 100){
        BVP CFEM(selection, hetero, Bnds);
        for (int i=0; i<NumOfElem.size(); i++){
            CFEM.add_CFEMGrid(NumOfElem[i], 1);
            CFEM.CFEM_1D(argc, args);
        }
    }
    else{
        NonLinear NLTest(selection, hetero, Bnds);
        for (int i = 0; i < NumOfElem.size(); i++) {
            NLTest.add_DFEMGrid(NumOfElem[i], 1);
            NLTest.NL_1D(argc, args);
            // L2Error_Vel[i] = NLTest.l2Err_Vel;
            // L2Error_Rho[i] = NLTest.l2Err_Rho;
            // H1Error_Vel[i] = NLTest.h1Err_Vel;
            // H1Error_Rho[i] = NLTest.h1Err_Rho;
            // h[i] = NLTest.info.hel;
        }
    }

    // printf("h\tL2Error");
    // for (int i=0; i<NumOfElem.size(); i++){
    //     printf("%.4e\t%.8e\t%.8e\n", h[i], L2Error_Vel[i], L2Error_Rho[i]);
    // }

    return 0;
}