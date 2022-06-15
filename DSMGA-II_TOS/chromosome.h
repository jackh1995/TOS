/* Copyright (C) 2011 by TEIL */

#ifndef _CHROMOSOME_H
#define _CHROMOSOME_H

#include <unordered_map>
#include "global.h"
#include "nk-wa.h"

using namespace std;

class Chromosome {

public:

    static enum Function {
        ONEMAX=0,
        MKTRAP=1,
        FTRAP=2,
        CYCTRAP=3,
        NK=4,
        SPINGLASS=5,
        SAT=6,
        MAXCUT=7,
        USal_NSize=8,
        USal_NSize_large=9,
        linear_mktrap=10,
        exponential_mktrap=11,
    } function;

    // static variables
    static int nfe;
    static int lsnfe;
    static int initnfe;
    static int rmnfe;
    static int bmnfe;
    static int bmeq;
    static int bmgt;
    static int bmfail;
    static bool hit;
    static unordered_map<unsigned long, double> cache;
    
    // constructor
    Chromosome ();
    Chromosome (int n_ell);
    
    // rule of three
    Chromosome(const Chromosome &);
    Chromosome& operator= (const Chromosome & c);
    ~Chromosome ();
    
    // initializers
    void init (int _ell);
    void init0 (int _ell);
    void initR (int _ell);
    
    // getters
    int getVal(int index) const;
    unsigned long getKey () const;
    double getFitness (int& counter);
    double getFitness () const;
    int getLength ();
    int getLength () const;
    int getLengthLong ();
    int getLengthLong () const;
    unsigned long* getGene();
    unsigned long* getGene() const;
    double getMaxFitness () const;
    unsigned int getCount(void) const;
    
    // setters
    void setVal (int index, int val);
    void setValF (int index, int val);

    // state
    bool isEvaluated () const;
    bool hasSeen() const;

    // evolutionary operators
    double evaluate (int& counter);
    void flip (int index);
    bool tryFlipping (int index);
    bool GHC();

    // test functions
    double oneMax () const;
    double trap (int u, double high, double low, int trapK) const;
    double mkTrap (double high, double low) const;
    double cycTrap(double fHigh, double fLow) const;
    double fTrap () const;
    double spinGlass () const;
    double nkFitness() const;
    double satFitness() const;
    double maxcutFitness() const;
    double USal_NSize_fitness() const;
    double linear_mktrap_fitness() const;
    double exponential_mktrap_fitness() const;

    // comparison operators
    bool operator== (const Chromosome & c) const;
    bool operator>(const Chromosome&) const;
    bool operator<(const Chromosome&) const;
    
    // miscellaneous
    friend ostream& operator<<(ostream& os, const Chromosome& ch);

    #ifdef TRIMMING
    bool getTrim() const;
    void setTrim (bool state);
    #endif

    unsigned long *gene;
    int lengthLong;

private:
    
    int length;
    double fitness;
    bool evaluated;
    unsigned long key;

    #ifdef TRIMMING
    bool trim;
    #endif
};

#endif
