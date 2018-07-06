#ifndef EigSolve
#define EigSolve
#include<vector>

using namespace std;
class matrix2D
{
public:
    void Init(int m, int n);
    vector<vector<double> > hess(int m, int n);
private:
    std::vector<std::vector<double> > mat;
};

#endif
