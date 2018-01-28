// Module Name
%module GaussInt_XT

%{
#include "Mesh_InterfaceXT.h"
%}

%include "std_vector.i"

namespace std {
%template(LineInt)     vector < int >;
%template(LineDouble)  vector < double >;
%template(ArrayInt)    vector < vector < int> >;
%template(ArrayDouble) vector < vector < double> >;
}

%include "Mesh_InterfaceXT.h"
