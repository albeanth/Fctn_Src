#include "EigSolve.h"
#include <iostream>
#include <math.h>
#include <stdlib.h>
using namespace std;

void matrix2D::Init(int m, int n)
{
  cout << "\nInitializing A.mat\n";
  mat.resize(m, vector<double>(n));

  for (int i=0; i<m; i++)
  {
      for (int j=0; j<n; j++)
      {
          mat[i][j] = i*m+1+j; // assign values to A
          cout << mat[i][j] << '\t';
      }
      cout << endl;
  }
  cout << endl;
}

double norm( vector<double> Vec, int size)
{
  double tmp = 0.0;
  for (int i=0; i<size; i++)
  {
    tmp += pow(Vec[i],2);
  }
  tmp = sqrt(tmp);
  return tmp;
}

vector<vector<double> > OuterProd(vector<double> x, vector<double> y, vector<vector<double> > prod)
{
  int m; int n;
  m = x.size(); n = y.size();
  // cout << m << ", " << n << endl;
  prod.resize(m, vector<double>(n));
  for (int i=0; i<m; i++)
  {
    for (int j=0; j<n; j++)
    {
      prod[i][j] = x[i]*y[j];
    }
  }
  return prod;
}

vector<vector<double> > MatProd(vector<vector<double> > A, vector<vector<double> > B)
{
  // cout << "In MatProd\n";
  int Am; int An;
  int Bm; int Bn;
  Am = A.size(); An = A[0].size();
  Bm = B.size(); Bn = B[0].size();
  // cout << "A = " << A.size() << ", "<< A[0].size() << endl;
  // cout << "B = " << B.size() << ", "<< B[0].size() << endl;
  vector< vector<double> > prod(Am, vector<double>(Bn));
  // cout << "\n";
  for (int row = 0; row < Am; row++)
  {
    for (int col=0; col < Bn; col++)
    {
      for (int inner=0; inner < An; inner ++)
      {
        prod[row][col] += A[row][inner] * B[inner][col];
      }
      // cout << prod[row][col] << "\t" ;
    }
    // cout << "\n";
  }
  return prod;
}

vector<vector<double> > SubMat(vector<vector<double> > mat, vector<vector<double> > sub, int m, int n)
{
  // cout << "Obtaining Submat\n";
  int p; int q;
  p = sub.size(); q = sub[0].size();
  // cout << p << ", " << q << endl;
  int i; int j;
  i = 0;
  for (int u=m-p; u<m; u++)
  {
    j = 0;
    for (int v=n-q; v<n; v++)
    {
      sub[i][j]=mat[u][v];
      // cout << sub[i][j] << "\t";
      j += 1;
    }
    // cout << "\n";
    i += 1;
  }
  return sub;
}


vector<vector<double> > matrix2D::hess(int m, int n)
{
  cout << "Computing Hessenberg form\n";

  std::vector<double> x;
  std::vector<int> eye;
  int xSize;
  double sign;
  double tmp;
  vector< vector<double> > v;

  for (int k=0; k<n-2; k++)
  {
    cout << "k=" << k << endl;
    x.clear(); //clears out x
    eye.clear(); // clears out eye
    for (int i=0; i<(n-(k+1)); i++) // fill x vec and canonical unit vec, eye
    {
      x.push_back(mat[k+1+i][k]);
      if (i == 0)
      {
        eye.push_back(1);
      }
      eye.push_back(0);
    }
    xSize = x.size(); // size of x vec and eye

    if (x[0] >= 0) // get sign of first element of x vec
    {
      sign = 1.0;
    }
    else
    {
      sign = -1.0;
    }

    tmp = norm(x,xSize);
    vector<double> dum(xSize);
    for (int i=0; i<xSize; i++)
    {
      dum[i] = sign*tmp*eye[i]+x[i];
    }

    v.push_back(std::vector<double>(xSize)); // appends another vector
    tmp = norm(dum,xSize);
    for (int i=0; i<xSize; i++)
    {
      v[k][i]=dum[i]/tmp;
      // cout << v[k][i] << endl;
    }

    // start doing matrix transformations
    vector<vector<double> > submat; // submatrix of mat used in transformations
    vector<vector<double> > vprod; // outer product of v vectors
    vector<vector<double> > dumMat; // the submatrix that is subtracted from mat

    vprod = OuterProd(v[k],v[k],vprod); // computes outer product of vectors v
    dumMat.clear();
    dumMat.resize(m-(k+1), vector<double>(n));
    dumMat = SubMat(mat,dumMat,m,n); // obtains submatrix of mat
    dumMat = MatProd(vprod, dumMat); //takes matrix product of outer prod and submat

    int smRow;
    smRow = 0;
    for (int row=k+1; row<m; row++)
    {
      for (int col=k; col<n; col++)
      {
        mat[row][col] -= 2*dumMat[smRow][col];
      }
      smRow+=1;
    }

    dumMat.clear();
    dumMat.resize(m, vector<double>(n-(k+1)));
    dumMat = SubMat(mat,dumMat,m,n); // obtains submatrix of mat
    dumMat = MatProd(dumMat, vprod);

    int smCol;
    for (int row=0; row<m; row++)
    {
      smCol = 0;
      for (int col=k+1; col<n; col++)
      {
        mat[row][col] -= 2*dumMat[row][smCol];
        smCol +=1;
      }
    }
  }

  //print out Hessenberg matrix
  for (int i=0; i<m; i++)
  {
    for (int j=0; j<n; j++)
    {
      printf("%.15e   ",mat[i][j]); //for some reason, printing here with long-double prints nonsense.
    }
    printf("\n");
  }
  return mat;
}
