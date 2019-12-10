#include "BVP.hpp"
#include "Nonlinear.hpp"
#include <iostream>
#include <vector>

int main(int argc,char **args){
    const std::vector<int> NumOfElem {64};
    /*
     * test #'s < 100
     *   CFEM using petsc ksp libary tests
     * test #'s \geq 100
     *   nonlinear solver using petsc snes library tests
     */
    const int selection {100};
    const bool hetero {false};
    const std::vector<double> Bnds {0.0, 1.0};

    if (selection < 100){
        std::vector<double> L2Error (NumOfElem.size(), 0.0);
        std::vector<double> H1Error (NumOfElem.size(), 0.0);
        BVP CFEM(selection, hetero);
        for (int i=0; i<NumOfElem.size(); i++){
            CFEM.add_CFEMGrid(NumOfElem[i], 1, Bnds);
            CFEM.CFEM_1D(argc, args);
        }
    }
    else{
        NonLinear NLTest(selection, hetero);
        NLTest.NL_1D(argc, args);
    }

    return 0;
}