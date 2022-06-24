#pragma once
#include <vector>
struct USal_NSize_instance{
    double opt;
    int n_bb;
    std::vector<int> bb_vector;
};

double trap(int n_ones, double fHigh, double fLow, int bb_size);
double evaluate_USal_NSize(int *x, USal_NSize_instance *inst);
void load_USal_NSize(char *cnf_file_name, USal_NSize_instance *inst);