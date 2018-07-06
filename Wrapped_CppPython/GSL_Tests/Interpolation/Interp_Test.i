// Module Name
%module InterpTest

%{
#include "Interp_Test.h"
%}

%include "std_vector.i"
%include "std_string.i"

namespace std {
%template(LineInt)     vector < int >;
%template(LineDouble)  vector < double >;
%template(ArrayInt)    vector < vector < int> >;
%template(ArrayDouble) vector < vector < double> >;
}

%include "Interp_Test.h"
