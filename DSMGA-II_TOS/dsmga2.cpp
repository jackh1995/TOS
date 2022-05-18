/***************************************************************************
 *   Copyright (C) 2015 Tian-Li Yu and Shih-Huan Hsu                       *
 *   tianliyu@ntu.edu.tw                                                   *
 ***************************************************************************/

#include "dsmga2.h"
#include "chromosome.h"
#include "fastcounting.h"
#include "statistics.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <vector>
#include <unordered_set>
#include <numeric>
using namespace std;

/* ------------------------------- CONSTRUCTOR ------------------------------ */

DSMGA2::DSMGA2(int _ell, int _nInitial, int _maxGan, int _maxFe, int fffff) :
    
    // stack memory allocation
    ell(_ell),
    nCurrent((_nInitial / 2) * 2),
    selectionPressure(2),
    maxGen(_maxGan),
    maxFe(_maxFe),
    repeat((ell > 50) ? ell / 50 : 1),
    generation(0),
    bestIndex(-1),
    
    // heap memory allocation
    masks(new vector<int>[ell]),
    selectionIndex(new int[nCurrent]),
    orderN(new int[nCurrent]),
    orderELL(new int[ell]),
    population(new Chromosome[nCurrent]),
    fastCounting(new FastCounting[ell]),

    orderRM(new int[nCurrent])

    #ifdef TRIMMING
    ,rm_success_sizes_mean(_ell/2)
    #endif
 
    #ifdef ORDERING
    ,orderILS(_ell/2) 
    #endif
    
    {

    /* ----------------------------- INITIALIZATION ----------------------------- */
    
    // global initialization
    Chromosome::function = (Chromosome::Function)fffff;
    Chromosome::nfe = 0,
    Chromosome::lsnfe = 0,
    Chromosome::rmnfe = 0,
    Chromosome::bmnfe = 0,
    Chromosome::bmeq = 0,
    Chromosome::bmgt = 0,
    Chromosome::bmfail = 0,
    Chromosome::initnfe = 0,
    Chromosome::hit = false,
    
    // dsm initialization
    graph.init(ell);
    graph_size.init(ell);
    
    // fastCounting initialization
    for (int i = 0; i < ell; i++)
        fastCounting[i].init(nCurrent);
    
    // population hashing (ch -> fitness)
    pHash.clear();
    for (int i = 0; i < nCurrent; ++i) {
        population[i].initR(ell);
        double f = population[i].getFitness(Chromosome::initnfe);
        pHash[population[i].getKey()] = f;
    }
    
    #ifdef ORDERING
    iota(orderILS.begin(), orderILS.end(), 0);
    #endif
    
    /* -------------------------- GREEDY HILL CLIMBING -------------------------- */

    // greedy hill climbing
    if (GHC) {
        for (int i = 0; i < nCurrent; i++)
            population[i].GHC();
    }
    
}

/* ------------------------------- DESTRUCTOR ------------------------------- */

DSMGA2::~DSMGA2() {
    delete[] masks;
    delete[] orderN;
    delete[] orderELL;
    delete[] selectionIndex;
    delete[] population;
    delete[] fastCounting;
}

/* ----------------------------------- RUN ---------------------------------- */

int DSMGA2::doIt () {
    
    generation = 0;
    
    while (!shouldTerminate()) {
        if (verbose == POPULATION || verbose == RM) {
            
            countPop();
            
            printPop();
        }
        oneRun();
    }
    
    return generation;
}

void DSMGA2::oneRun() {
    
    if (CACHE) {
        Chromosome::cache.clear();
    }
    
    mixing();
    
    #ifdef ORDERING
    reset_orderILS();
    #endif
    
    double max = -INF;
    stFitness.reset();
    
    for (int i = 0; i < nCurrent; ++i) {
        
        double fitness = population[i].getFitness();
        
        if (fitness > max) {
            max = fitness;
            bestIndex = i;
        }
        
        stFitness.record(fitness);
    }

    if (verbose != NO) {
        showStatistics();
    }
    
    ++generation;
}

bool DSMGA2::shouldTerminate() {
    
    bool termination = false;
    
    // reach maximum nfe
    if (maxFe != -1) {
        if (Chromosome::nfe > maxFe) {
            termination = true;
        }
    }

    // reach maximum generation
    if (maxGen != -1) {
        if (generation > maxGen) {
            termination = true;
        }
    }

    // reach optimial fitness
    if (population[0].getMaxFitness() - EPSILON <= stFitness.getMax()) {
        termination = true;
    }
    
    // population convergence
    if (stFitness.getMax() - EPSILON <= stFitness.getMean()) {
        termination = true;
    }
        
    return termination;
}

bool DSMGA2::foundOptima() {
    return (population[0].getMaxFitness() - EPSILON <= stFitness.getMax());
}

/* -------------------------------- SELECTION ------------------------------- */

void DSMGA2::selection() {
    tournamentSelection();
}

void DSMGA2::tournamentSelection() {
    
    int i, j;
    int randArray[selectionPressure * nCurrent];
    
    for (i = 0; i < selectionPressure; i++) {
        myRand.uniformArray(randArray + (i * nCurrent), nCurrent, 0, nCurrent - 1);
    }
    
    for (i = 0; i < nCurrent; i++) {

        int winner = 0;
        double winnerFitness = -INF;

        for (j = 0; j < selectionPressure; j++) {

            int challenger = randArray[selectionPressure * i + j];
            double challengerFitness = population[challenger].getFitness();

            if (challengerFitness > winnerFitness) {
                winner = challenger;
                winnerFitness = challengerFitness;
            }
        }

        selectionIndex[i] = winner;
    }
}

/* ----------------------------------- DSM ---------------------------------- */

void DSMGA2::buildFastCounting() {
    
    // build dsm from the after-selection population
    if (SELECTION) {
        for (int i = 0; i < nCurrent; i++) {
            for (int j = 0; j < ell; j++) {
                fastCounting[j].setVal(i, population[selectionIndex[i]].getVal(j));
            }
        }
    } 
    
    // build dsm without selection
    else {
        for (int i = 0; i < nCurrent; i++) {
            for (int j = 0; j < ell; j++) {
                fastCounting[j].setVal(i, population[i].getVal(j));
            }
        }
    }
}

// count the number of ones at gene x
int DSMGA2::countOne(int x) const {
    
    int n = 0;
    
    for (int i = 0; i < fastCounting[0].lengthLong; ++i) {
        unsigned long val = 0;
        val = fastCounting[x].gene[i];
        n += myBD.countOne(val);
    }

    return n;
}

// compute g[x] XOR g[y], where x and y are gene indices
int DSMGA2::countXOR(int x, int y) const {
    
    int n = 0;
    
    for (int i = 0; i < fastCounting[0].lengthLong; ++i) {
        unsigned long val = 0;
        val = fastCounting[x].gene[i];
        val ^= fastCounting[y].gene[i];
        n += myBD.countOne(val);
    }

    return n;
}

double DSMGA2::computeMI(double a00, double a01, double a10, double a11) const {

    double p0 = a00 + a01;
    double q0 = a00 + a10;
    double p1 = 1 - p0;
    double q1 = 1 - q0;
    double join = 0.0;

    if (a00 > EPSILON)
        join += a00 * log(a00);
    if (a01 > EPSILON)
        join += a01 * log(a01);
    if (a10 > EPSILON)
        join += a10 * log(a10);
    if (a11 > EPSILON)
        join += a11 * log(a11);

    double p = 0.0;

    if (p0 > EPSILON)
        p += p0 * log(p0);
    if (p1 > EPSILON)
        p += p1 * log(p1);

    double q = 0.0;

    if (q0 > EPSILON)
        q += q0 * log(q0);
    if (q1 > EPSILON)
        q += q1 * log(q1);

    return -p - q + join;
}

void DSMGA2::buildGraph() {
    
    int *one = new int[ell];
    
    for (int i = 0; i < ell; ++i) {
        one[i] = countOne(i);
    }
    
    for (int i = 0; i < ell; ++i) {
        for (int j = i + 1; j < ell; ++j) {
    
            int n00, n01, n10, n11;
            int nX = countXOR(i, j);
    
            n11 = (one[i] + one[j] - nX) / 2;
            n10 = one[i] - n11;
            n01 = one[j] - n11;
            n00 = nCurrent - n01 - n10 - n11;
    
            double p00 = (double)n00 / (double)nCurrent;
            double p01 = (double)n01 / (double)nCurrent;
            double p10 = (double)n10 / (double)nCurrent;
            double p11 = (double)n11 / (double)nCurrent;
            double p1_ = p10 + p11;
            double p0_ = p00 + p01;
            double p_0 = p00 + p10;
            double p_1 = p01 + p11;
            double linkage = computeMI(p00, p01, p10, p11);
            double linkage00 = 0.0, linkage01 = 0.0;
    
            if (p00 > EPSILON)
                linkage00 += p00 * log(p00 / p_0 / p0_);
            if (p11 > EPSILON)
                linkage00 += p11 * log(p11 / p_1 / p1_);
            if (p01 > EPSILON)
                linkage01 += p01 * log(p01 / p0_ / p_1);
            if (p10 > EPSILON)
                linkage01 += p10 * log(p10 / p1_ / p_0);
    
            if (Chromosome::nfe < 0) {
                pair<double, double> p(linkage, linkage);
                graph.write(i, j, p);
            } else {
                pair<double, double> p(linkage00, linkage01);
                graph.write(i, j, p);
            }
        }
    }
    
    delete[] one;
}

void DSMGA2::buildGraph_sizecheck() {

    int *one = new int[ell];

    for (int i = 0; i < ell; ++i) {
        one[i] = countOne(i);
    }

    for (int i = 0; i < ell; ++i) {
        for (int j = i + 1; j < ell; ++j) {

            int n00, n01, n10, n11;
            int nX = countXOR(i, j);

            n11 = (one[i] + one[j] - nX) / 2;
            n10 = one[i] - n11;
            n01 = one[j] - n11;
            n00 = nCurrent - n01 - n10 - n11;

            double p00 = (double)n00 / (double)nCurrent;
            double p01 = (double)n01 / (double)nCurrent;
            double p10 = (double)n10 / (double)nCurrent;
            double p11 = (double)n11 / (double)nCurrent;
            double p1_ = p10 + p11;
            double p0_ = p00 + p01;
            double p_0 = p00 + p10;
            double p_1 = p01 + p11;
            double linkage = computeMI(p00, p01, p10, p11);
            double linkage00 = 0.0, linkage01 = 0.0;

            if (p00 > EPSILON)
                linkage00 += p00 * log(p00 / p_0 / p0_);
            if (p11 > EPSILON)
                linkage00 += p11 * log(p11 / p_1 / p1_);
            if (p01 > EPSILON)
                linkage01 += p01 * log(p01 / p0_ / p_1);
            if (p10 > EPSILON)
                linkage01 += p10 * log(p10 / p1_ / p_0);

            pair<double, double> p(linkage, linkage);
            graph_size.write(i, j, p);
        }
    }

    delete[] one;
}

/* --------------------------------- MIXING --------------------------------- */

void DSMGA2::findMask(Chromosome &ch, vector<int> &result, int startNode) {
    
    result.clear();
    DLLA rest(ell);
    genOrderELL();
    
    for (int i = 0; i < ell; i++) {
        if (orderELL[i] == startNode)
            result.push_back(orderELL[i]);
        else
            rest.insert(orderELL[i]);
    }
    
    double *connection = new double[ell];
    
    for (DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++) {
    
        pair<double, double> p = graph(startNode, *iter);
        int i = ch.getVal(startNode);
        int j = ch.getVal(*iter);
        
        // p00 or p11
        if (i == j) 
            connection[*iter] = p.first;
    
        // p01 or p10
        else 
            connection[*iter] = p.second;
    }
    
    while (!rest.isEmpty()) {
    
        double max = -INF;
        int index = -1;
    
        for (DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++) {
            if (max < connection[*iter]) {
                max = connection[*iter];
                index = *iter;
            }
        }
    
        rest.erase(index);
        result.push_back(index);
    
        for (DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++) {
    
            pair<double, double> p = graph(index, *iter);
            int i = ch.getVal(index);
            int j = ch.getVal(*iter);
    
            // p00 or p11
            if (i == j) 
                connection[*iter] += p.first;
    
            // p01 or p10
            else 
                connection[*iter] += p.second;
        }
    }
    
    delete[] connection;
}

void DSMGA2::findMaskWithBound(Chromosome &ch, vector<int> &result, int startNode, int bound) {

    result.clear();
    
    DLLA rest(ell);
    
    for (int i = 0; i < ell; i++) {

        if (orderELL[i] == startNode)
            result.push_back(orderELL[i]);
        else
            rest.insert(orderELL[i]);
    }
    
    double *connection = new double[ell];
    
    for (DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++) {
        
        pair<double, double> p = graph_size(startNode, *iter);
        int i = ch.getVal(startNode);
        int j = ch.getVal(*iter);
    
        // p00 or p11
        if (i == j) 
            connection[*iter] = p.first;
    
        // p01 or p10
        else 
            connection[*iter] = p.second;
    }
    
    bound--;
    
    while (!rest.isEmpty() && bound > 0) {
    
        bound--;
    
        double max = -INF;
    
        int index = -1;
    
        for (DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++) {
            if (max < connection[*iter]) {
                max = connection[*iter];
                index = *iter;
            }
        }
    
        rest.erase(index);
    
        result.push_back(index);
    
        for (DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++) {
        
            pair<double, double> p = graph_size(index, *iter);
            int i = ch.getVal(index);
            int j = ch.getVal(*iter);
    
            // p00 or p11
            if (i == j) 
                connection[*iter] += p.first;
    
            // p01 or p10
            else 
                connection[*iter] += p.second;
        }
    }
    
    delete[] connection;
}

size_t DSMGA2::findSize(Chromosome &ch, vector<int> &mask) {
    
    DLLA candidate(nCurrent);
    
    for (int i = 0; i < nCurrent; ++i)
        candidate.insert(i);
    
    size_t size = 0;
    
    for (vector<int>::iterator it = mask.begin(); it != mask.end(); ++it) {
    
        int allele = ch.getVal(*it);
    
        for (DLLA::iterator it2 = candidate.begin(); it2 != candidate.end(); ++it2) {
    
            if (get_ch_dist(ch, population[*it2]) > ell/2 && USE_HAMMING) {
                candidate.erase(*it2);
            }
            else {
                if (population[*it2].getVal(*it) == allele)
                    candidate.erase(*it2);
            }
            
            // #ifdef SIMILARITY_CHECK
            // if (get_ch_dist(ch, population[*it2]) > ell/2) {
            //     candidate.erase(*it2);
            // }
            // #endif
            
            // if (population[*it2].getVal(*it) == allele)
            //     candidate.erase(*it2);

            if (candidate.isEmpty())
                break;
        }

        if (candidate.isEmpty())
            break;
        
        ++size;
    }
    
    return size;
}

void DSMGA2::mixing() {
    
    if (SELECTION) {
        selection();
    }
    
    buildFastCounting();
    
    buildGraph();
    
    buildGraph_sizecheck();

    // k epochs
    for (int k = 0; k < repeat; ++k) {
        
        genOrderRM();

        for (int i = 0; i < nCurrent/2; ++i) {
            
            restrictedMixing(population[orderRM[i]]);
            
            if (Chromosome::hit)
                break;
        }

        if (Chromosome::hit)
            break;
    
        #ifdef TRIMMING
        if (!rm_success_sizes.empty()){
        
        size_t sum(0);
        vector<int>::iterator it;
        
        for (it = rm_success_sizes.begin(); it != rm_success_sizes.end(); ++it){
            sum += (*it);
        }
        
        rm_success_sizes_mean = sum / rm_success_sizes.size() + 1;
        }
        #endif

        for (int i=nCurrent/2; i<nCurrent; ++i) {
         
            restrictedMixing(population[orderRM[i]]);
         
            if (Chromosome::hit) 
                break;
        }        

        #ifdef TRIMMING
            rm_success_sizes.clear();
        #endif
        
        if (Chromosome::hit)
            break;
    }
}

bool DSMGA2::restrictedMixing(Chromosome &ch, vector<int> &mask) {
    
    bool taken(false);
    size_t lastUB(0);
    Chromosome trial(ell);
    trial = ch;
    
    // the mask has passed the donor check
    for (size_t idx = 0; idx != mask.size(); ++idx) {
        
        size_t ub;
        #ifdef ORDERING
        ub = orderILS[idx]+1;
        #else
        ub = idx + 1;
        #endif

        #ifdef ORDERING
        if (ub > mask.size()) {
            continue;
        }
        #endif

        size_t size = 1;
        Chromosome trial(ell);
        trial = ch;
	    vector<int> takenMask;
        
        // flip the trial solution
        for (vector<int>::iterator it = mask.begin(); it != mask.end(); ++it) {
			
            takenMask.push_back(*it);
            trial.flip(*it);
            ++size;
            
            if (size > ub) {
                break;
            }
        }
        
        // discard the trial solution if we have seen it before
        if (isInP(trial)) {
            break;
        }
        
        // accept the trial solution if it is new and has better fitness
        if (trial.getFitness(Chromosome::rmnfe) >= ch.getFitness(Chromosome::rmnfe) - EPSILON) {
            pHash.erase(ch.getKey());
            pHash[trial.getKey()] = trial.getFitness();
            taken = true;
            ch = trial;
        }
        
        if (taken) {
        
            lastUB = ub;

            #ifdef ORDERING
            update_orderILS(lastUB);
            #endif

            #ifdef TRIMMING
            if (!ch.getTrim()) {
                rm_success_sizes.push_back(lastUB);
            }
            #endif
            
            break;
        }
    }
    
    // print rm status
    if (verbose == RM)
        if (mask.size()) {
            
            unordered_set<int> success;
            unordered_set<int> fail;
            
            for (size_t idx=0; idx!=mask.size(); ++idx) {
                if (idx < lastUB)
                    success.insert(mask[idx]);
                else
                    fail.insert(mask[idx]);
            }
            
            for (size_t idx=0; idx!=size_t(ell); ++idx) {
                if (success.count(idx)) {
                    cout << GREEN << idx%10 << RESET;
                }
                else if (fail.count(idx)) {
                    cout << RED << idx%10 << RESET;
                }
                else {
                    cout << idx%10;
                }
            } 

            cout << ' ' << lastUB << '/' << mask.size();
            
            if (!taken) {
                cout << endl;
            }
        }
    
    // rm succeeds
    if (lastUB != 0) {
        while (mask.size() > lastUB) {
            mask.pop_back();
        }
    }

    return taken;
}

void DSMGA2::restrictedMixing(Chromosome &ch) {
        
    int startNode = myRand.uniformInt(0, ell - 1);

    vector<int> mask;
    findMask(ch, mask, startNode);
    size_t size = findSize(ch, mask);

    vector<int> mask_size;
    findMaskWithBound(ch, mask_size, startNode, size);
    size_t size_original = findSize(ch, mask_size);
    
    if (size > size_original) {
        size = size_original;
    }
    
    if (size > (size_t)ell / 2) {
        size = ell / 2;
    }

    size_t size_max(0);

    #ifdef TRIMMING
    if (rm_success_sizes.empty()) {
        size_max = ell/2;
    } else {
        if (ch.getTrim()) {
            size_max = rm_success_sizes_mean;
        }
        else {
            size_max = ell/2;
        }
    }
    size_max = (size_max < (size_t)ell/2) ? size_max : ell/2;
    #else
    size_max = ell/2;
    #endif
    
    if (size > size_max) {
        size = size_max;
    }

    while (mask.size() > size) {
        mask.pop_back();
    }
    
    // restricted mixing
    bool taken = restrictedMixing(ch, mask);
    
    // Back mixing
    EQ = true;
    
    if (taken) {

        Chromosome::bmeq=0;
        Chromosome::bmgt=0;
        Chromosome::bmfail=0;

        genOrderN();

        for (int i = 0; i < nCurrent; ++i) {
        
            #ifdef SIMILARITY_CHECK
            if (get_ch_dist(ch, population[orderN[i]]) <= ell/2) {
                if (EQ) {
                    backMixingE(ch, mask, population[orderN[i]]);
                } else {
                    backMixing(ch, mask, population[orderN[i]]);
                }
            }
            #else
            if (EQ) {
                backMixingE(ch, mask, population[orderN[i]]);
            } else {
                backMixing(ch, mask, population[orderN[i]]);
            }
            #endif
        }
        
        // print rm status
        if (verbose == RM) {
            cout << ' ' << GREEN << Chromosome::bmgt << RESET << '/' \
            << YELLOW << Chromosome::bmeq << RESET << '/'  << \
            RED << Chromosome::bmfail << RESET << endl;
        }
    }
}

void DSMGA2::backMixing(Chromosome &source, vector<int> &mask, Chromosome &des) {

    Chromosome trial(ell);
    trial = des;

    for (vector<int>::iterator it = mask.begin(); it != mask.end(); ++it) {
        trial.setVal(*it, source.getVal(*it));
    }
        
    if (trial.getFitness(Chromosome::bmnfe) > des.getFitness(Chromosome::bmnfe)) {
        pHash.erase(des.getKey());
        pHash[trial.getKey()] = trial.getFitness();
        des = trial;
        ++Chromosome::bmgt;
        return;
    }

    ++Chromosome::bmfail;
}

void DSMGA2::backMixingE(Chromosome &source, vector<int> &mask, Chromosome &des) {
    Chromosome trial(ell);
    trial = des;

    for (vector<int>::iterator it = mask.begin(); it != mask.end(); ++it)
        trial.setVal(*it, source.getVal(*it));

    // strong bm
    if (trial.getFitness(Chromosome::bmnfe) > des.getFitness(Chromosome::bmnfe)) {
        pHash.erase(des.getKey());
        pHash[trial.getKey()] = trial.getFitness();
        EQ = false;
        des = trial;
        ++Chromosome::bmgt;
        return;
    }

    // weak bm
    if (trial.getFitness(Chromosome::bmnfe) >= des.getFitness(Chromosome::bmnfe) - EPSILON) {
        pHash.erase(des.getKey());
        pHash[trial.getKey()] = trial.getFitness();
        des = trial;
        ++Chromosome::bmeq;
        return;
    }

    ++Chromosome::bmfail;
}

bool DSMGA2::isInP(const Chromosome &ch) const {
    unordered_map<unsigned long, double>::const_iterator it = pHash.find(ch.getKey());
    return (it != pHash.end());
}

/* ---------------------------- RANDOM TRAVERSAL ---------------------------- */

void DSMGA2::genOrderN() {
    myRand.uniformArray(orderN, nCurrent, 0, nCurrent - 1);
}

void DSMGA2::genOrderELL() {
    myRand.uniformArray(orderELL, ell, 0, ell - 1);
}

/* ---------------------------------- PRINT --------------------------------- */

void DSMGA2::showStatistics() {
    
    printf("Gen:%d  Fitness:(Max/Mean/Min):%f/%f/%f \n",
           generation, 
           stFitness.getMax(), 
           stFitness.getMean(),
           stFitness.getMin()
    );
    
    fflush(NULL);
}

// count the number of a specific pattern
void DSMGA2::countPop(void) {
    
    ch_count.clear();
    
    for (int idx=0; idx!=nCurrent; ++idx) {
        
        // has seen
        if (ch_count.count(population[idx].getKey())) {
            ++ch_count[population[idx].getKey()];
        }

        // has not seen
        else {
            ch_count[population[idx].getKey()] = 1;
        }
    } 
    
    cout << endl;
}

// beautifully print the whole population
void DSMGA2::printPop(void) {
    
    population_set.clear();
    
    // insert every chromosome to the formatted population set
    for (int idx=0; idx!=nCurrent; ++idx) {
        population_set.insert(population[idx]);
    }

    set<Chromosome>::iterator it;

    // loop through every elememt in the formatted population set    
    for (it = population_set.begin(); it != population_set.end(); ++it) {
        cout << *it << ' ' << it->getCount() << endl;
    }
}

/* --------------------------------- GETTERS -------------------------------- */

const Chromosome* const DSMGA2::getPopulation(void) const {
    return population;
}

int DSMGA2::getGeneration() const {
    return generation;
}

int const DSMGA2::getBestIndex(void) const {
    return bestIndex;
}

void DSMGA2::genOrderRM() {
    
    // shuffle the order to trim
    myRand.uniformArray(orderRM, nCurrent, 0, nCurrent-1);
    
    #ifdef TRIMMING
    double probTrim(0.5);
    int numTrim(int(probTrim * nCurrent));
    
    // for every chromosome in the popoulation
    for (int idx=0; idx!=nCurrent; ++idx){
        
        // for the latter 50%, the max mask size is ell/2
        if (idx < numTrim){
            population[orderRM[idx]].setTrim(false);
        } 
        
        // trim the first 50% population
        else {
            population[orderRM[idx]].setTrim(true);
        }
    }
    #endif
}

#ifdef ORDERING
void DSMGA2::print_orderILS(void) const {
    
    for (int idx=0; idx!=ell/2; ++idx) {
        cout << orderILS[idx] << ' ';
    } 
    
    cout << endl;
}

void DSMGA2::reset_orderILS(void) {
    std::iota(orderILS.begin(), orderILS.end(), 0);
}

void DSMGA2::update_orderILS(size_t lastUB) {

    int last_UB_idx(0);
    
    for (int idx=0; idx!=ell/2; ++idx) {
        if (orderILS[idx] == (int(lastUB)-1)) {
            last_UB_idx = idx;
            break;
        }
    }

    if (last_UB_idx == 0) {
        return;
    }

    swap(orderILS[last_UB_idx-1], orderILS[last_UB_idx]);
}
#endif

int DSMGA2::get_ch_dist(Chromosome x, Chromosome y) {
    
    int result(0);
    
    // for (int i=0; i!=ell; ++i) {
    //     if (x.getVal(i) != y.getVal(i))
    //         ++result;
    // }
    
    // assert(result >=0 && result <= ell);

    for (int i=0; i<x.lengthLong; ++i) {
        result += myBD.countOne(x.gene[i] ^ y.gene[i]);
    }
    
    return result;
}