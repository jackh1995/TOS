#ifndef _maxcut_h_
#define _maxcut_h_
#include <vector>

struct MAXCUTinstance{

    double opt;
    int ell,n,m;
    std::vector<int> fvector;

};

double evaluateMAXCUT(char*, MAXCUTinstance*);
void loadMAXCUT(char*,char*, MAXCUTinstance*);

#endif