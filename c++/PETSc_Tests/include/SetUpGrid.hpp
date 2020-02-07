// INCLUDE GUARD
#ifndef __SetUpGrid_H_INCLUDED__
#define __SetUpGrid_H_INCLUDED__

// headers for dependents (what this class needs to function)
#include "UtilityFunctions.hpp" // <- parent class

class SetUpGrid : public UtilityFunctions{
  /*
  This class defines the characteristics of a 1D mesh (uniform grid).
  The following attributes are considered:
      1) nels   = number of elements
      2) maxord = maximum order of elements
      3) order  = list of polynomial order for corresponsing elements
      4) nod    = list of global numbering of nodes
      5) xnod   = list of node coordinates
      6) nnodes = notal number of nodes
  */
  public:
    SetUpGrid(const std::vector<double> &Bnds){
      info.bounds = Bnds;
      info.nnodes = NAN;
    };
    // Member variables
    FEMGrid info;
    // Member Functions
    void add_CFEMGrid(const int nels, const int myorder);
    void add_DFEMGrid(const int nels, const int myorder, const double delta = std::numeric_limits<double>::epsilon());
};

#endif
