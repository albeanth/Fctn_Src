#include "Mesh_Interface.h"
#include <stdlib.h>

using namespace std;

// void GaussianIntegration::TestFun(std::vector<double> v,std::vector<double> w,int z){
//   cout << "print v" << endl;
//   for (int i=0; i<v.size(); i++){
//     cout << v[i] << endl;
//   }
//   cout << "print w" << endl;
//   for (int i=0; i<w.size(); i++){
//     cout << w[i] << endl;
//   }
//   cout << "print z" << endl << z;
// }

void GaussianIntegration::ImportMeshGridInfo(std::vector<double> v,std::vector<double> w,int z){
  cout << "print v" << endl;
  for (int i=0; i<v.size(); i++){
    cout << v[i] << endl;
  }
  cout << "print w" << endl;
  for (int i=0; i<w.size(); i++){
    cout << w[i] << endl;
  }
  cout << "print z" << endl << z;
}
