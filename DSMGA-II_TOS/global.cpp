/***************************************************************************
 *   Copyright (C) 2015 Tian-Li Yu and Shih-Huan Hsu                       *
 *   tianliyu@ntu.edu.tw                                                   *
 ***************************************************************************/

#include <cstdio>
#include <cstdlib>
#include <climits>
#include <cfloat>
#include "myrand.h"
#include "statistics.h"
#include "doublelinkedlistarray.h"
#include "zkey.h"
#include "chromosome.h"
#include "sat.h"
#include "dsmga2.h"

#ifdef SIMILARITY_CHECK
bool USE_HAMMING = true;
#endif

int maxMemory = 0;

bool GHC = true;
bool SELECTION = true;
bool CACHE = false;
bool SHOW_BISECTION = false;
bool SHOW_MP = true;

char outputFilename[100];
/* ------------------------ SECTION GLOBAL VARIABLES ------------------------ */
Chromosome::Function Chromosome::function;
int Chromosome::nfe;
int Chromosome::lsnfe;
int Chromosome::rmnfe;
int Chromosome::bmnfe;
int Chromosome::initnfe;
int Chromosome::bmeq;
int Chromosome::bmgt;
int Chromosome::bmfail;
bool Chromosome::hit;
unordered_map<unsigned long, unsigned> DSMGA2::ch_count;
set<Chromosome> DSMGA2::population_set;
Verbosity verbose;
/* -------------------------------------------------------------------------- */

unordered_map<unsigned long, double> Chromosome::cache;
ZKey zKey;
MyRand myRand;
BitwiseDistance myBD;
SPINinstance mySpinGlassParams;
NKWAProblem nkwa;
SATinstance mySAT;
MAXCUTinstance myMAXCUT;

void outputErrMsg(const char *errMsg) {
    printf("%s\n", errMsg);
    exit(1);
}

int pow2(int x) {
    return (1 << x);
}

