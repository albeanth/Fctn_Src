// Module Name
%module GaussInt_ThWave

%{
#include "ThermalWave_Interface.h"
%}

%include "std_vector.i"

namespace std {
%template(LineInt)     vector < int >;
%template(LineDouble)  vector < double >;
%template(ArrayInt)    vector < vector < int> >;
%template(ArrayDouble) vector < vector < double> >;
}

%include "ThermalWave_Interface.h"
