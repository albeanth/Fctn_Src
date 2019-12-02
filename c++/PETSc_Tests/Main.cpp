#include "BVP.hpp"
#include <iostream>
#include <vector>

int main(int argc,char **args){
    const std::vector<int> NumOfElem {32};
    const int selection {0};
    const bool hetero {false};
    const std::vector<double> Bnds {0.0, 1.0};

    std::vector<double> L2Error (NumOfElem.size(), 0.0);
    std::vector<double> H1Error (NumOfElem.size(), 0.0);
    BVP CFEM(selection, hetero);
    for (int i=0; i<NumOfElem.size(); i++){
        CFEM.add_CFEMGrid(NumOfElem[i], 1, Bnds);
        CFEM.CFEM_1D(argc, args);
    }


    return 0;
}