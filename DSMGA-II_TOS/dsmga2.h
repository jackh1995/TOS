/***************************************************************************
 *   Copyright (C) 2015 Tian-Li Yu and Shih-Huan Hsu                       *
 *   tianliyu@ntu.edu.tw                                                   *
 ***************************************************************************/

#ifndef _DSMGA2_H_
#define _DSMGA2_H_

#include "chromosome.h"
#include "doublelinkedlistarray.h"
#include "fastcounting.h"
#include "statistics.h"
#include "trimatrix.h"
#include <list>
#include <unordered_map>
#include <set>

class DSMGA2 {

public:

    static unordered_map<unsigned long, unsigned> ch_count;

    // constructor
    DSMGA2(int n_ell, int n_nInitial, int n_maxGen, int n_maxFe, int fffff);
    
    // destructor
    ~DSMGA2();

    // run
    void oneRun();
    int doIt();
    bool shouldTerminate();
    bool foundOptima();
    
    // selection
    void selection();
    void tournamentSelection();
    
    // getters
    const Chromosome* const getPopulation(void) const;
    int getGeneration() const;
    int const getBestIndex(void) const;

private:

    static set<Chromosome> population_set;

    int ell;
    int nCurrent; 
    int selectionPressure;
    int maxGen;
    int maxFe;
    int repeat;
    int generation;
    int bestIndex;
    
    vector<int> *masks;
    int *selectionIndex;
    int *orderN; 
    int *orderELL;
    Chromosome *population;
    FastCounting *fastCounting;
    Statistics stFitness;
    unordered_map<unsigned long, double> pHash;
    bool EQ;

    // dsm
    TriMatrix<double> graph;
    TriMatrix<double> graph_size;

    void buildFastCounting();
    int countXOR(int, int) const;
    int countOne(int) const;
    double computeMI(double, double, double, double) const;
    void buildGraph();
    void buildGraph_sizecheck();
    
    // mixing
    void findMask(Chromosome &ch, vector<int> &mask, int startNode);
    void findMaskWithBound(Chromosome &ch, vector<int> &mask, int startNode, int bound);
    size_t findSize(Chromosome &, vector<int> &);
    void mixing();
    void restrictedMixing(Chromosome &);
    bool restrictedMixing(Chromosome &ch, vector<int> &mask);
    void backMixing(Chromosome &source, vector<int> &mask, Chromosome &des);
    void backMixingE(Chromosome &source, vector<int> &mask, Chromosome &des);
    bool isInP(const Chromosome &) const;
    
    // random traversal
    void genOrderN();
    void genOrderELL();

    // print
    void showStatistics();
    void countPop(void);
    void printPop(void);
    
    int *orderRM;
    void genOrderRM();
    #ifdef TRIMMING
    int rm_success_sizes_mean;
    vector<int> rm_success_sizes;
    #endif

    #ifdef ORDERING
    vector<int> orderILS;
    void print_orderILS(void) const;
    void reset_orderILS(void);
    void update_orderILS(size_t);
    #endif
    
    // [HAMMING]
    int get_ch_dist(Chromosome x, Chromosome y);
};

#endif