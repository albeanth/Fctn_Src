// Module Name
%module FEMInt

%{
#include "FEM_Integration.h"
%}

%include "std_vector.i"

namespace std {
%template(LineInt)     vector < int >;
%template(LineDouble)  vector < double >;
%template(ArrayInt)    vector < vector < int> >;
%template(ArrayDouble) vector < vector < double> >;
}

%include "FEM_Integration.h"
