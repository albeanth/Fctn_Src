// Module Name
%module GaussInt

%{
#include "Mesh_Interface.h"
%}

%include "std_vector.i"

namespace std {
%template(Line)  vector < double >;
}

%include "Mesh_Interface.h"
