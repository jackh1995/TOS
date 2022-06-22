#include "traps.h"
#include <iostream>
#include <iomanip>
#include <cstring>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
using namespace std;

double trap(int n_ones, double fHigh, double fLow, int bb_size) {
    if (n_ones > bb_size)
        return 0;

    if (n_ones == bb_size) {
        return fHigh;
    } else {
        return fLow - n_ones * fLow / (bb_size-1);
    }
}

double evaluate_USal_NSize(int *x, USal_NSize_instance *inst) {
    size_t current_gene(0);
    int u;
    double gap;
    int res(0);

    for(vector<int>::const_iterator it=inst->bb_vector.begin(); it!=inst->bb_vector.end(); it++) {
        // cout << '(' << *it << ')';
        u = 0;
        gap = 1 - 1.0/double(*it);
        for (int i=0; i!=*it; ++i) {
            u += x[current_gene + i];
        }
        res += trap(u, 1, gap, *it);
        // cout << trap(u, 1, gap, *it)<< '/';
        current_gene += *it;
    }
    // cout << "\n\n";
    return res;
}

void load_USal_NSize(char *cnf_file_name, USal_NSize_instance *inst) {
	int temp;
	ifstream input;
	string line;

    input.open(cnf_file_name);
    getline (input, line);
	istringstream in(line);
	in >> inst->opt;
	in >> inst->n_bb;

	getline(input, line);
    istringstream ins(line);
    while (ins >> temp) {
        inst->bb_vector.push_back(temp);
    }
}
