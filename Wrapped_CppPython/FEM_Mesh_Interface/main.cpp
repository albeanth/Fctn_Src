#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "Mesh_Interface.h"
using namespace std;


int main()
{
  // Get instance of class
  GaussianIntegration mytest;

  // Initialize dummy mesh data
  char flag;
  int nels;
  vector<int> order;
  vector<vector<int> > nod;
  vector<double> xnod;
  int maxord;

  flag = 'x';
  nels = 100;
  for (int i=0; i<nels; i++){
    order.push_back(2);
  }
  for (int i=0; i<nels; i++ ){
    nod.push_back(vector<int>() );
    for (int j=0; j<2; j++){
      nod[i].push_back(i+j);
    }
  }

  for (int i=0; i<nels; i++){
    xnod.push_back(2);
  }
  maxord = 2;

  // Initialize dummy FE coeffs
  vector<vector<double> > X; vector<vector<double> > Y; vector<vector<double> > T;
  for (int i=0; i<2; i++ ){
    X.push_back(vector<double>() );
    Y.push_back(vector<double>() );
    T.push_back(vector<double>() );
    for (int j=0; j<nels; j++){
      X[i].push_back((double)rand());
      Y[i].push_back((double)rand());
      T[i].push_back((double)rand());
    }
  }

  // test source integrator
  vector<double> v;
  v = mytest.Source_Integrate(flag,
          nels, order, nod, xnod, maxord,
          nels, order, nod, xnod, maxord,
          nels, order, nod, xnod, maxord);

  // test error integrator
  vector<double> error;
  error = mytest.Error_Integrate3D(X, Y, T,
            nels, order, nod, xnod, maxord,
            nels, order, nod, xnod, maxord,
            nels, order, nod, xnod, maxord);
}
