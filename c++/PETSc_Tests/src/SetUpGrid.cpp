#include "SetUpGrid.hpp"

void SetUpGrid::add_CFEMGrid(const int nels, const int myorder, const std::vector<double> &Bnds){
  // variable declarations
  double a; double b; //bounds
  double hel; double hmax; // element width (constant)
  std::vector<double> xel; // element bounds
  int nnodes; // number of nodes
  int maxord;

  a = Bnds[0];
  b = Bnds[1];

  hel = (b-a)/nels; // spacing of elements (uniform)
  xel = arange(a,b,hel);
  hmax = hel;

  std::vector<int> order(nels,myorder+1);
  maxord = *std::max_element(order.begin(), order.end());

  int tmp = 0;
  for (int elem : order){
    tmp += elem - 1;
  }
  nnodes = tmp + 1;

  // derive the global indexing of nodes:
  //nod(i,1) is the global number of j'th node in element i
  // std::vector<std::vector<int>> nod(nels, std::vector<int> (maxord-1,0));
  std::vector<std::vector<int>> nod;
  int n = 0;
  for (int k=0; k<nels; k++){
    std::vector<int> tmp;
    for (int j=0; j<order[k]; j++){
      tmp.push_back(n);
      if (j != order[k]-1){
        n+=1;
      }
    }
    nod.push_back(tmp);
  }
  // uncomment below to see the way nodes are numbered
  // for (int row=0; row<nels; row++){
  //   for (int col=0; col<maxord; col++){
  //     printf("%d  ",nod[row][col]);
  //   }
  //   printf("\n");
  // }

  // xnod , i=1..nnodes
  //  -> coordinates of node i
  std::vector<double> xnod(nnodes,0);
  double h; double hi;
  for (int k=0; k<nels-1; k++){
    h = xel[k+1]-xel[k];
    hi = h/(order[k]-1);
    for (int j=0; j<order[k]; j++){
      xnod[nod[k][j]] = xel[k] + hi*j;
    }
  }
  int k = nels-1;
  h = b-xel[k];
  hi=h/(order[k]-1);
  for (int j=0; j<order[k]; j++){
      xnod[nod[k][j]] = xel[k] + hi*j;
    }
  // for (double elem : xnod){
  //   printf("%.4f\n", elem);
  // }

  info.nels = nels;
  info.maxord = maxord;
  info.order = order;
  info.nod = nod;
  info.xnod = xnod;
  info.nnodes = nnodes;
  info.bounds = Bnds;
}

void SetUpGrid::add_DFEMGrid(const int nels, const int myorder, const std::vector<double> &Bnds, const double delta){
  /*
   * sets up a discretized discontinuous finite element domain over the interval 
   * defined by Bnds. 
   * - uses a default discontinuity spacing of the order of floating point precision
   *   for the machine calling this function
   *   INPUT:
   *     [a,b]   = bounds of calculation domain
   *     xnel    =  the number or location of elements
   *               -> if one number (scalar) then xnel = number of elements,
   *               -> if a vector, then it contains position of first node of each element
   *     myorder =  degree of polynomials per element
   *               -> if one number (scalar) then it is uniform over all elements,
   *               -> if a vector, then it must be as long as the vector of elements
   *     delta   = space for discontinuity between cells. Default value = machine epsilon.
   *   OUTPUT:
   *     mesh object with necessary attributes
   */
  printf("\nSetting up DFEM grid");
  info.nels = nels;
  // variable declarations
  const double a {Bnds[0]};
  const double b {Bnds[1]};
  const double hel {(b-a)/info.nels}; // spacing of elements (uniform)
  std::vector<double> xel = arange(a,b,hel);

  std::vector<int> order(info.nels,myorder+1);
  info.maxord = {*std::max_element(order.begin(), order.end())};
  info.nnodes = sum(order);

  // derive the global indexing of nodes:
  //info.nod(i,1) is the global number of j'th node in element i
  // std::vector<std::vector<int>> info.nod(info.nels, std::vector<int> (info.maxord-1,0));
  int n = 0;
  std::vector<int> tmp;
  for (int k=0; k<info.nels; k++){
    tmp.clear();
    for (int j=0; j<order[k]; j++){
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
    hi = h/(order[k]-1);
    for (int j=0; j<order[k]; j++){
      if (j==0)
        info.xnod[info.nod[k][j]] = xel[k] + delta;
      else
        info.xnod[info.nod[k][j]] = info.xnod[info.nod[k][j-1]] + hi;
    }
  }
  const int k = info.nels-1;
  h = b - xel[k] - (2*delta);
  hi = h/(order[k]-1);
  for (int j=0; j<order[k]; j++){
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

  info.order = order;
}

