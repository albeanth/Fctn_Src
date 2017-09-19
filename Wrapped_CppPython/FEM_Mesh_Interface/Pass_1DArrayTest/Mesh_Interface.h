#ifndef Mesh_Interface
#define Mesh_Interface_h 1

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>

/*
This interface allows for python generated 1D FEM meshes to be read into c++ via swig.
*/

class GaussianIntegration{

  public:
    // double GaussInt_1D(double xMin, double xMax, int xNode);
    void ImportMeshGridInfo(std::vector<double> v,std::vector<double> w,int z);

  private:
    // double TestFun1D(double x);
    // std::vector<double> linspace(double a, double b, int NumElem);
};

#endif
