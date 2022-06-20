/* Copyright (C) 2015 by TEIL */

#include <cstdio>
#include <cstring>
#include "spin.h"
#include "chromosome.h"
#include "nk-wa.h"
#include "sat.h"
#include <iostream>
#include <bitset>
#include "dsmga2.h"
#include <cmath>

using namespace std;

#define TRAP_K 5

/* ------------------------------- CONSTRUCTOR ------------------------------ */

Chromosome::Chromosome () {

    length = 0;
    lengthLong = 0;
    gene = NULL;
    evaluated = false;

    #ifdef TRIMMING
    trim = false;
    #endif
}

Chromosome::Chromosome (int n_length) {
    
    gene = NULL;

    // initialization
    init(n_length);
}

/* ------------------------------ RULE OF THREE ----------------------------- */

// copy constructor
Chromosome::Chromosome(const Chromosome &c) {
    
    length = c.length;
    evaluated = c.evaluated;
    fitness = c.fitness;
    lengthLong = c.lengthLong;
    key = c.key;
    gene = new unsigned long[c.lengthLong];
    
    memcpy(gene, c.gene, sizeof(long) * lengthLong);
}

// copier
Chromosome& Chromosome::operator= (const Chromosome& c) {

    if (length != c.length) {
        length = c.length;
        init(length);
    }

    evaluated = c.evaluated;
    fitness = c.fitness;
    lengthLong = c.lengthLong;
    key = c.key;

    memcpy(gene, c.gene, sizeof(long) * lengthLong);

    return *this;
}

// destructor
Chromosome::~Chromosome () {
    if (gene != NULL) {
        delete []gene;
    }
}


/* ------------------------------- INITIALIZER ------------------------------ */

void Chromosome::init (int _length) {
    
    length = _length;
    lengthLong = quotientLong(length)+1;

    if (gene != NULL)
        delete []gene;

    gene = new unsigned long [lengthLong];
    gene[lengthLong-1] = 0;

    evaluated = false;
}

// zero initializaion
void Chromosome::init0 (int _length) {
    
    length = _length;
    lengthLong = quotientLong(length)+1;

    if (gene != NULL)
        delete []gene;

    gene = new unsigned long [lengthLong];

    for (int i=0; i<lengthLong; ++i)
        gene[i] = 0;

    key = 0;
    evaluated = false;
}

// random initialization
void Chromosome::initR (int _length) {
    
    length = _length;
    lengthLong = quotientLong(length)+1;

    if (gene != NULL)
        delete []gene;

    gene = new unsigned long [lengthLong];
    gene[lengthLong-1] = 0;

    key = 0;
    for (int i=0; i<length; ++i) {

        int val = myRand.flip();
        setValF(i, val);
        if (val == 1)
            key ^= zKey[i];
    }

    evaluated = false;
}

/* --------------------------------- GETTERS -------------------------------- */

// get the fitness of a chromosome
double Chromosome::getFitness (int& counter) {
    
    if (evaluated) {
        return fitness;
    } else {
      
        fitness = evaluate(counter);
      
        if (!hit && fitness > getMaxFitness()) {
            hit = true;
        }

        return fitness;
    }
}

double Chromosome::getFitness() const {
    if (evaluated) {
        return fitness;
    } else {
        throw runtime_error("Error: Chromosome not evaluated");
    }
}

int Chromosome::getLength () const {
    return length;
}

int Chromosome::getLength () {
    return length;
}

int Chromosome::getLengthLong () const { 
    return lengthLong;
}

int Chromosome::getLengthLong () {
    return lengthLong;
}

unsigned long * Chromosome::getGene() {
    return gene;
}    

unsigned long * Chromosome::getGene() const {
    return gene;
}

int Chromosome::getVal (int index) const {
        assert (index >= 0 && index < length);

        int q = quotientLong(index);
        int r = remainderLong(index);

        if ( (gene[q] & (1lu << r)) == 0 )
            return 0;
        else
            return 1;
    }

unsigned long Chromosome::getKey () const {
    return key;
}
// get the theoretical optimal fitness
double Chromosome::getMaxFitness () const {

    double maxF;

    switch (function) {
        case ONEMAX:
            maxF = length;
            break;
        case MKTRAP:
            maxF = length/TRAP_K;
            break;
        case FTRAP:
            maxF = length/6;
            break;
        case CYCTRAP:
            maxF =  length/(TRAP_K - 1);
            break;
        case SPINGLASS:
            maxF = mySpinGlassParams.opt;
            break;
        case NK:
            maxF = nkwa.maxF;
            break;
        case SAT:
            maxF = 0;
            break;
        case MAXCUT:
            maxF = myMAXCUT.opt;
            break;
        case USal_NSize:
            maxF = my_USal_NSize.opt;
            break;
        case USal_NSize_large:
            maxF = my_USal_NSize.opt;
            break;
        case linear_mktrap:
            maxF = (1 + length/TRAP_K) * length/TRAP_K / 2;
            break;
        case exponential_mktrap:
            maxF = pow2(length/TRAP_K) - 1;
            break;
        case four_five_six:
            // TRAP_K := 5
            maxF = length / 5;
            break;
        case three_four_five_six_seven:
            // TRAP_K := 5
            maxF = length / 5;
            break;
        case ftrap4:
            maxF = length/6;
            break;
        case ftrap6:
            maxF = length/6;
            break;
        default: 
            maxF = INF;
    }
    
    return maxF - EPSILON;
}

unsigned int Chromosome::getCount(void) const {
    return DSMGA2::ch_count[getKey()];
}

/* --------------------------------- SETTERS -------------------------------- */

void Chromosome::setVal (int index, int val) {

    assert (index >= 0 && index < length);

    if (getVal(index) == val) return;

    setValF(index, val);
    key ^= zKey[index];
}

void Chromosome::setValF (int index, int val) {

    assert (index >= 0 && index < length);
    //outputErrMsg ("Index overrange in Chromosome::operator[]");

    int q = quotientLong(index);
    int r = remainderLong(index);

    if (val == 1)
        gene[q] |= (1lu<<r);
    else
        gene[q] &= ~(1lu<<r);

    evaluated = false;
}

/* ---------------------------------- STATE --------------------------------- */

// whether the chromosome is evaluated or not
bool Chromosome::isEvaluated () const {
    return evaluated;
}

// whether the pattern has seen before
bool Chromosome::hasSeen() const {

    unordered_map<unsigned long, double>::iterator it = cache.find(key);
    
    if (it != cache.end()) {
        return true;
    } else {
        return false;
    }
}

/* ------------------------- EVOLUTIONARY OPERATORS ------------------------- */

// evaluate a chromosome
double Chromosome::evaluate(int& counter) {
    
    if (CACHE)
        if (hasSeen()) {
            evaluated = true;
            return cache[key];
        }
    
    if (!Chromosome::hit) {
        ++nfe;
        ++counter;
    }
        
    evaluated = true;
    double accum = 0.0;

    switch (function) {
        case ONEMAX:
            accum = oneMax();
            break;
        case MKTRAP:
            accum = mkTrap(1, 0.8);
            break;
        case CYCTRAP:
            accum = cycTrap(1, 0.8);
            break;
        case FTRAP:
            accum = fTrap();
            break;
        case SPINGLASS:
            accum = spinGlass();
            break;
        case NK:
            accum = nkFitness();
            break;
        case SAT:
            accum = satFitness();
            break;
        case MAXCUT:
            accum = maxcutFitness();
            break;
        case USal_NSize:
            accum = USal_NSize_fitness();
            break;
        case USal_NSize_large:
            accum = USal_NSize_fitness();
            break;
        case linear_mktrap:
            accum = linear_mktrap_fitness();
            break;
        case exponential_mktrap:
            accum = exponential_mktrap_fitness();
            break;
        case four_five_six:
            accum = four_five_six_fitness();
            break;
        case three_four_five_six_seven:
            accum = three_four_five_six_seven_fitness();
            break;
        case ftrap4:
            accum = ftrap4_fitness();
            break;
        case ftrap6:
            accum = ftrap6_fitness();
            break;
        default:
            accum = mkTrap(1, 0.8);
            break;
    }

    if (CACHE)
        cache[key]=accum;

    return accum;
}

void Chromosome::flip (int index) {
    assert (index >= 0 && index < length);

    int q = quotientLong(index);
    int r = remainderLong(index);

    gene[q] ^= (1lu<<r);
    key ^= zKey[index];

    evaluated = false;
}

// contribute to lsnfe
bool Chromosome::tryFlipping(int index) {
    
    // old fitness
    double oldF = getFitness(Chromosome::lsnfe);

    // flip the gene
    flip(index);

    // not improve
    if (getFitness(Chromosome::lsnfe) - EPSILON <= oldF) {

        flip(index);

        evaluated = true;
        fitness = oldF;
        return false;
    } 
    
    // improve
    else {
        return true;
    }
}

// greedy hill climbing
bool Chromosome::GHC() {

    int* order = new int [length];
    myRand.uniformArray(order, length, 0, length-1);

    bool flag = false;

    for (int i=0; i<length; ++i) {
        if (tryFlipping(order[i])) {
            flag = true;
        }
    }

    delete []order;
    return flag;
}

/* ----------------------------- TEST FUNCTIONS ----------------------------- */

double Chromosome::oneMax () const {

    double result = 0;

    for (int i = 0; i < length; ++i)
        result += getVal(i);

    return result;
}

double Chromosome::trap (int unitary, double fHigh, double fLow, int trapK) const {
    
    if (unitary > trapK)
        return 0;

    if (unitary == trapK) {
        return fHigh;
    } else {
        return fLow - unitary * fLow / (trapK-1);
    }
}

double Chromosome::mkTrap (double fHigh, double fLow) const {

    int i, j;
    int u;
    int TRAP_M = length / TRAP_K;

    if (length % TRAP_K != 0)
        outputErrMsg ("TRAP_K doesn't divide length");

    double result = 0;

    for (i = 0; i < TRAP_M; i++) {

        u = 0;

        for (j = 0; j < TRAP_K; j++) {
            u += getVal(i * TRAP_K + j);
        }

        result += trap (u, fHigh, fLow, TRAP_K);
    }

    return result;
}

double Chromosome::cycTrap(double fHigh, double fLow) const {

    int i, j;
    int u;
    int TRAP_M = length / (TRAP_K-1);

    if (length % (TRAP_K-1) != 0) {
        outputErrMsg ("TRAP_k doesn't divide length for Cyclic Setting");
    }
        
    double result = 0;

    for (i = 0; i < TRAP_M; i++) {
    
        u = 0;
        int idx = i * TRAP_K - i;
    
        for (j = 0; j < TRAP_K; j++) {
    
            int pos = idx + j;
    
            if (pos == length) {
                pos = 0;
            } else if (pos > length) {
                outputErrMsg ("CYCLIC BUG");
            }
                
            u += getVal(pos);
        }
        
        result += trap (u, fHigh, fLow, TRAP_K);
    }
    
    return result;
}

double Chromosome::fTrap() const {

    double result = 0.0;

    for (int i=0; i<length/6; ++i) {

        int u=0;

        for (int j=0; j<6; ++j)
            u += getVal(i*6+j);

        if (u==0)
            result += 1.0;
        else if (u==1)
            result += 0.0;
        else if (u==2)
            result += 0.4;
        else if (u==3)
            result += 0.8;
        else if (u==4)
            result += 0.4;
        else if (u==5)
            result += 0.0;
        else // u == 6
            result += 1.0;
    }

    return result;
}

double Chromosome::ftrap4_fitness() const {

    double result = 0.0;

    for (int i=0; i<length/6; ++i) {

        int u=0;

        for (int j=0; j<6; ++j)
            u += getVal(i*6+j);

        if (u==0)
            result += 1.0;
        else if (u==1)
            result += 0.0;
        else if (u==2)
            result += 0.4;
        else if (u==3)
            result += 0.8;
        else if (u==4)
            result += 0.4;
        else if (u==5)
            result += 0.0;
        else // u == 6
            result += 1.0;
    }

    double result2 = 0.0;

    for (int i=0; i<length/6; ++i) {

        int u=0;

        // optimal = 000111
        for (int j=0; j<3; ++j)
            u += 1 - getVal(i*6+j);
        for (int j=3; j<6; ++j)
            u += getVal(i*6+j);

        if (u==0)
            result2 += 1.0;
        else if (u==1)
            result2 += 0.0;
        else if (u==2)
            result2 += 0.4;
        else if (u==3)
            result2 += 0.8;
        else if (u==4)
            result2 += 0.4;
        else if (u==5)
            result2 += 0.0;
        else // u == 6
            result2 += 1.0;
    }
    
    if (result > result2) {
        return result;
    } else {
        return result2;
    }
}

double Chromosome::ftrap6_fitness() const {

    double result = 0.0;

    for (int i=0; i<length/6; ++i) {

        int u=0;

        for (int j=0; j<6; ++j)
            u += getVal(i*6+j);

        if (u==0)
            result += 1.0;
        else if (u==1)
            result += 0.0;
        else if (u==2)
            result += 0.4;
        else if (u==3)
            result += 0.8;
        else if (u==4)
            result += 0.4;
        else if (u==5)
            result += 0.0;
        else // u == 6
            result += 1.0;
    }

    double result2 = 0.0;

    for (int i=0; i<length/6; ++i) {

        int u=0;

        // optimal = 000111
        for (int j=0; j<3; ++j)
            u += 1 - getVal(i*6+j);
        for (int j=3; j<6; ++j)
            u += getVal(i*6+j);

        if (u==0)
            result2 += 1.0;
        else if (u==1)
            result2 += 0.0;
        else if (u==2)
            result2 += 0.4;
        else if (u==3)
            result2 += 0.8;
        else if (u==4)
            result2 += 0.4;
        else if (u==5)
            result2 += 0.0;
        else // u == 6
            result2 += 1.0;
    }

    double result3 = 0.0;

    for (int i=0; i<length/6; ++i) {

        int u=0;

        // optimal = 010101
        for (int j=0; j<6; ++j)
            if (j % 2) {
                u += 1 - getVal(i*6+j);
            } else {
                u += getVal(i*6+j);
            }

        if (u==0)
            result3 += 1.0;
        else if (u==1)
            result3 += 0.0;
        else if (u==2)
            result3 += 0.4;
        else if (u==3)
            result3 += 0.8;
        else if (u==4)
            result3 += 0.4;
        else if (u==5)
            result3 += 0.0;
        else // u == 6
            result3 += 1.0;
    }
    
    if (result >= result2 && result >= result3) {
        return result;
    } else if (result2 >= result && result2 >= result3) {
        return result2;
    } else { // result3 >= result && result3 >= result2
        return result3;
    }
}

double Chromosome::spinGlass () const {

    int *x = new int[length];
    double result;

    for (int i=0; i<length; i++)
        if (getVal(i) == 1)
            x[i] = 1;
        else
            x[i] = -1;

    result = evaluateSPIN(x, &mySpinGlassParams);

    delete[] x;

    return result;
}

double Chromosome::nkFitness() const {
    
    char *x = new char[length];

    for ( int i = 0; i < length; ++i) {
        x[i] = (char) getVal(i);
    }

    double result = evaluateNKProblem(x, &nkwa);
    
    delete[] x;
    
    return result;
}

double Chromosome::satFitness() const {
    
    int *x = new int[length];

    for ( int i = 0; i < length; ++i) {
        x[i] = getVal(i);
    }

    double result = evaluateSAT(x, &mySAT);
    
    delete []x;
    
    return result;
}

double Chromosome::maxcutFitness() const {

    char *x = new char[length];

    for (int i=0; i!=length; ++i){
        x[i] = getVal(i);
    }

    double result = evaluateMAXCUT(x, &myMAXCUT);

    delete[] x;

    return result;
}

double Chromosome::USal_NSize_fitness() const {

    int *x = new int[length];

    for (int i=0; i!=length; ++i){
        x[i] = getVal(i);
    }

    double result = evaluate_USal_NSize(x, &my_USal_NSize);

    delete[] x;

    return result;
}

double Chromosome::linear_mktrap_fitness() const {

    int i, j;
    int u;
    int TRAP_M = length / TRAP_K;

    if (length % TRAP_K != 0)
        outputErrMsg ("TRAP_K doesn't divide length");

    double result = 0;

    double fHigh = 1.0;
    double fLow = 0.8;

    for (i = 0; i < TRAP_M; i++) {

        u = 0;

        for (j = 0; j < TRAP_K; j++) {
            u += getVal(i * TRAP_K + j);
        }

        result += trap (u, fHigh, fLow, TRAP_K);
        fHigh += 1.0;
        fLow = fHigh * 4.0 / 5.0;
    }

    return result;
}

double Chromosome::exponential_mktrap_fitness() const {

    int i, j;
    int u;
    int TRAP_M = length / TRAP_K;

    if (length % TRAP_K != 0)
        outputErrMsg ("TRAP_K doesn't divide length");

    double result = 0;

    double fHigh = 1.0;
    double fLow = 0.8;

    for (i = 0; i < TRAP_M; i++) {

        u = 0;

        for (j = 0; j < TRAP_K; j++) {
            u += getVal(i * TRAP_K + j);
        }

        result += trap (u, fHigh, fLow, TRAP_K);
        fHigh *= 2.0;
        fLow = fHigh * 4.0 / 5.0;
    }

    return result;
}

double Chromosome::four_five_six_fitness() const {
    
    int i, j;
    int u;
    int milestone;

    if (length % 15 != 0)
        outputErrMsg ("length must be a multiple of 3 * 5");

    double result = 0.0;

    for (i = 0; i < length / 15; i++) {
        u = 0;
        for (j = 0; j < 4; j++) {
            u += getVal(i * 4 + j);
        }
        result += trap (u, 1.0, 0.75, 4);
    }

    milestone = length / 15 * 4;

    for (i = 0; i < length / 15; i++) {
        u = 0;
        for (j = 0; j < 5; j++) {
            u += getVal(milestone + i * 5 + j);
        }
        result += trap (u, 1.0, 0.8, 5);
    }

    milestone = length / 15 * 9;

    for (i = 0; i < length / 15; i++) {
        u = 0;
        for (j = 0; j < 6; j++) {
            u += getVal(milestone + i * 6 + j);
        }
        result += trap (u, 1.0, 1.0 - 1.0 / 6.0, 6);
    }

    return result;
}

double Chromosome::three_four_five_six_seven_fitness() const {
    
    int i, j;
    int u;
    int milestone;

    if (length % 25 != 0)
        outputErrMsg ("length must be a multiple of 5 * 5");

    double result = 0.0;

    for (i = 0; i < length / 25; i++) {
        u = 0;
        for (j = 0; j < 3; j++) {
            u += getVal(i * 3 + j);
        }
        result += trap (u, 1.0, 1.0 - 1.0 / 3.0, 3);
    }

    milestone = length / 25 * 3;
    for (i = 0; i < length / 25; i++) {
        u = 0;
        for (j = 0; j < 4; j++) {
            u += getVal(milestone + i * 4 + j);
        }
        result += trap (u, 1.0, 0.75, 4);
    }

    milestone = length / 25 * 7;
    for (i = 0; i < length / 25; i++) {
        u = 0;
        for (j = 0; j < 5; j++) {
            u += getVal(milestone + i * 5 + j);
        }
        result += trap (u, 1.0, 0.8, 5);
    }

    milestone = length / 25 * 12;
    for (i = 0; i < length / 25; i++) {
        u = 0;
        for (j = 0; j < 6; j++) {
            u += getVal(milestone + i * 6 + j);
        }
        result += trap (u, 1.0, 1.0 - 1.0 / 6.0, 6);
    }

    milestone = length / 25 * 18;
    for (i = 0; i < length / 25; i++) {
        u = 0;
        for (j = 0; j < 7; j++) {
            u += getVal(milestone + i * 7 + j);
        }
        result += trap (u, 1.0, 1.0 - 1.0 / 7.0, 7);
    }

    return result;
}

/* -------------------------- COMPARISON OPERATORS -------------------------- */

// compare if two chromosomes are equal
bool Chromosome::operator== (const Chromosome& c) const {
    
    if (length != c.length)
        return false;

    for (int i=0; i<lengthLong; i++)
        if (gene[i] != c.gene[i])
            return false;

    return true;
}

bool Chromosome::operator>(const Chromosome& c) const{
    if (getKey() == c.getKey()) {
        return false;
    }    
    if (getFitness() > c.getFitness()) {
        return true;
    } else if (getFitness() == c.getFitness() && getCount() >= c.getCount()) {
        return true;
    } else {
        return false;
    }
}

bool Chromosome::operator<(const Chromosome& c) const{
    if (getKey() == c.getKey()) {
        return false;
    }
    return *this > c;
    if (getFitness() < c.getFitness()) {
        return true;
    } else if (getFitness() == c.getFitness() && getCount() <= c.getCount()) {
        return true;
    } else {
        return false;
    }
}

#ifdef TRIMMING
bool Chromosome::getTrim() const {
    return trim;
}
void Chromosome::setTrim (bool state) {
    trim = state;
}
#endif

/* ------------------------------ MISCELLANEOUT ----------------------------- */

// print a chromosome
ostream& operator<<(ostream& os, const Chromosome& ch){
    for (int i=0; i!=ch.getLength(); ++i) {
        if (ch.getVal(i)) {
            cout << GREEN << 1 << RESET;
        } else {
            cout << RED << 0 << RESET;
        }
    } 
    
    printf(" %6.3f", ch.getFitness());
    
    return os;
}
