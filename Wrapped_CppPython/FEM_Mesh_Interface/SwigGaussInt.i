// Module Name
%module GaussInt

%{
#include "Mesh_InterfaceXYT.h"
%}

%include "std_vector.i"

namespace std {
%template(LineInt)     vector < int >;
%template(LineDouble)  vector < double >;
%template(ArrayInt)    vector < vector < int> >;
%template(ArrayDouble) vector < vector < double> >;
}

%include "Mesh_InterfaceXYT.h"
