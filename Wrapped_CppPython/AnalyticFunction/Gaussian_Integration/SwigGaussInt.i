
// Module Name
%module GaussInt

// ************************************************************
// Module Includes
// ************************************************************

// These are copied directly to the .cxx file and are not parsed
// by SWIG.  Include include files or definitions that are required
// for the module to build correctly.
%{

#include "GaussianIntegration.h"

// Namespaces
using namespace std;

%}


%include "GaussianIntegration.h"
