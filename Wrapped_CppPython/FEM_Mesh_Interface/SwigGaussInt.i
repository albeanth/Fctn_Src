// Module Name
%module GaussInt_XTE

%{
#include "Mesh_InterfaceXTE.h"
%}

%include "std_vector.i"

namespace std {
%template(LineInt)     vector < int >;
%template(LineDouble)  vector < double >;
%template(ArrayInt)    vector < vector < int> >;
%template(ArrayDouble) vector < vector < double> >;
}

%include "Mesh_InterfaceXTE.h"
