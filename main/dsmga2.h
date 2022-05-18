/*
 * dsmga2.h
 *
 *  Created on: May 2, 2011
 *      Author: tianliyu
 */

#ifndef _DSMGA2_H_
#define _DSMGA2_H_

#include <vector>
#include "chromosome.h"
#include "statistics.h"
#include "trimatrix.h"
#include "doublelinkedlistarray.h"
#include "fastcounting.h"

class DSMGA2 {
public:
    DSMGA2 (int n_ell, int n_nInitial, int n_maxGen, int n_maxFe, int fffff);

    ~DSMGA2 ();

    void selection ();
    /** tournament selection without replacement*/
    void tournamentSelection();

    void oneRun (bool output = true);
    int doIt (bool output = true);

    void buildGraph ();
    void mixing ();
    // 2016-11-26
    void findMask(Chromosome& ch, vector<int>& mask,int startNode);
    void findMask_size(Chromosome& ch, vector<int>& mask,int startNode,int bound);
    void buildGraph_sizecheck();
    // ****
    void restrictedMixing(Chromosome&);
    bool restrictedMixing(Chromosome& ch, vector<int>& mask);
    void backMixing(Chromosome& source, vector<int>& mask, Chromosome& des);
    void backMixingE(Chromosome& source, vector<int>& mask, Chromosome& des);

    bool shouldTerminate ();

    bool foundOptima ();

    int getGeneration () const {
        return generation;
    }

    bool isInP(const Chromosome& ) const;
    void genOrderN();
    void genOrderELL();

    void showStatistics ();

    bool isSteadyState ();

//protected:
public:

    int ell;                                  // chromosome length
    int nCurrent;                             // population size
    bool EQ;
    unordered_map<unsigned long, double> pHash; // to check if a chromosome is in the population


    vector<int> *masks;
    int *selectionIndex;
    int *orderN;                             // for random order
    int *orderELL;                             // for random order
    int selectionPressure;
    int maxGen;
    int maxFe;
    int repeat;
    int generation;
    int bestIndex;

    Chromosome* population;
    FastCounting* fastCounting;

    TriMatrix<double> graph;
   // 2016-11-26
    TriMatrix<double> graph_size;
    double previousFitnessMean;
    Statistics stFitness;

    // methods
    double computeMI(double, double, double, double) const;


    void findClique(int startNode, vector<int>& result);

    void buildFastCounting();
    int countXOR(int, int) const;
    int countOne(int) const;

    size_t findSize(Chromosome&, vector<int>&);
    size_t findSize(Chromosome&, vector<int>&, Chromosome&);
    
    bool converged();
    double lastMax, lastMean, lastMin;
    int convergeCount;

    int *orderRM;
    void genOrderRM();
    #ifdef TRIMMING
    int rm_success_sizes_mean;
    vector<int> rm_success_sizes;
    #endif

    #ifdef ORDERING
    vector<int> orderILS;
    void reset_orderILS(void);
    void update_orderILS(size_t);
    #endif

    #ifdef SIMILARITY_CHECK
    int get_ch_dist(Chromosome x, Chromosome y);
    #endif
};


#endif /* _DSMGA2_H_ */
