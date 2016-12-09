#include <iostream>
#include <vector>
#include "EigSolve.h"
#include <stdlib.h>
using namespace std;
int main(int argc, char **argv)
{
    int m=4; int n=4;

    matrix2D A; //make A an object of class matrix2D
    A.Init(m,n); // initialize A.mat with data
    A.hess(m,n);

    return 0;
}
