#include "SetUpGrid.hpp"

void SetUpGrid::add_CFEMGrid(const int nels, const int myorder) {
  /*
   * sets up a discretized continuous finite element domain over the interval 
   * defined by info.Bnds.
   *   INPUT:
   *     myorder =  degree of polynomials per element
   *               -> if one number (scalar) then it is uniform over all elements,
   *               -> if a vector, then it must be as long as the vector of elements
   *     delta   = space for discontinuity between cells. Default value = machine epsilon.
   *   OUTPUT:
   *     mesh object with necessary attributes
   */
  printf("\nSetting up CFEM grid... ");
  info.nels = nels;
  double hel {(info.bounds[1]-info.bounds[0])/info.nels}; // spacing of elements (uniform)
  std::vector<double> xel {arange(info.bounds[0], info.bounds[1], hel)};

  info.order.assign(info.nels, myorder+1);
  info.maxord = *std::max_element(info.order.begin(), info.order.end());

  int tmp = 0;
  for (int elem : info.order){
    tmp += elem - 1;
  }
  info.nnodes = tmp + 1;

  // derive the global indexing of nodes:
  //info.nod(i,1) is the global number of j'th node in element i
  int n = 0;
  for (int k=0; k<info.nels; k++){
    std::vector<int> tmp;
    for (int j=0; j<info.order[k]; j++){
      tmp.push_back(n);
      if (j != info.order[k]-1){
        n+=1;
      }
    }
    info.nod.push_back(tmp);
  }
  // uncomment below to see the way nodes are numbered
  // for (int row=0; row<info.nels; row++){
  //   for (int col=0; col<info.maxord; col++){
  //     printf("%d  ",info.nod[row][col]);
  //   }
  //   printf("\n");
  // }

  // info.xnod , i=1..info.nnodes
  //  -> coordinates of node i
  info.xnod.resize(info.nnodes);
  double h; double hi;
  for (int k=0; k<info.nels-1; k++){
    h = xel[k+1]-xel[k];
    hi = h/(info.order[k]-1);
    for (int j=0; j<info.order[k]; j++){
      info.xnod[info.nod[k][j]] = xel[k] + hi*j;
    }
  }
  int k = info.nels-1;
  h = info.bounds[1]-xel[k];
  hi=h/(info.order[k]-1);
  for (int j=0; j<info.order[k]; j++){
      info.xnod[info.nod[k][j]] = xel[k] + hi*j;
    }
  // for (double elem : info.xnod){
  //   printf("%.4f\n", elem);
  // }
  printf("done!\n\n");

}

void SetUpGrid::add_DFEMGrid(const int nels, const int myorder, const double delta) {
  /*
   * sets up a discretized discontinuous finite element domain over the interval 
   * defined by info.Bnds.
   * - uses a default discontinuity spacing of the order of floating point precision
   *   for the machine calling this function
   *   INPUT:
   *     myorder =  degree of polynomials per element
   *               -> if one number (scalar) then it is uniform over all elements,
   *               -> if a vector, then it must be as long as the vector of elements
   *     delta   = space for discontinuity between cells. Default value = machine epsilon.
   *   OUTPUT:
   *     mesh object with necessary attributes
   */
  printf("\nSetting up DFEM grid... ");
  
  info.nels = nels;
  const double hel {(info.bounds[1]-info.bounds[0])/info.nels}; // spacing of elements (uniform)
  std::vector<double> xel {arange(info.bounds[0],info.bounds[1],hel)};

  info.order.assign(info.nels, myorder+1);
  info.maxord = {*std::max_element(info.order.begin(), info.order.end())};
  info.nnodes = sum(info.order);

  // derive the global indexing of nodes:
  //info.nod(i,1) is the global number of j'th node in element i
  // std::vector<std::vector<int>> info.nod(info.nels, std::vector<int> (info.maxord-1,0));
  int n {0};
  std::vector<int> tmp;
  for (int k=0; k<info.nels; k++){
    tmp.clear();
    for (int j=0; j<info.order[k]; j++){
      tmp.push_back(n);
      n+=1;
    }
    info.nod.push_back(tmp);
  }
  // uncomment below to see the way nodes are numbered
  // printf("\n  'nod'  numbering\n");
  // for (int row=0; row<info.nels; row++){
  //   for (int col=0; col<info.maxord; col++){
  //     printf("    %d  ",info.nod[row][col]);
  //   }
  //   printf("\n");
  // }

  // info.xnod , i=1..info.nnodes
  //  -> coordinates of node i
  info.xnod = zeros(info.nnodes);
  double h, hi;
  for (int k=0; k<info.nels-1; k++){
    h = (xel[k+1]-xel[k]) - (2.0*delta);
    // spacing between each interior node (based on polynomial order of that element)
    hi = h/(info.order[k]-1);
    for (int j=0; j<info.order[k]; j++){
      if (j==0)
        info.xnod[info.nod[k][j]] = xel[k] + delta;
      else
        info.xnod[info.nod[k][j]] = info.xnod[info.nod[k][j-1]] + hi;
    }
  }
  const int k {info.nels-1};
  h = info.bounds[1] - xel[k] - (2*delta);
  hi = h/(info.order[k]-1);
  for (int j=0; j<info.order[k]; j++){
    if (j == 0)
      info.xnod[info.nod[k][j]] = xel[k] + delta;
    else
      info.xnod[info.nod[k][j]] = info.xnod[info.nod[k][j-1]] + hi;
  }
  // printf("\n  Xnod locations\n");
  // for (double elem : info.xnod){
  //   printf("    %.4e\n", elem);
  // }
  // exit(-1);
  printf("done!\n\n");
}

