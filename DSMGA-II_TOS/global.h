/***************************************************************************
 *   Copyright (C) 2015 Tian-Li Yu and Shih-Huan Hsu                       *
 *   tianliyu@ntu.edu.tw                                                   *
 ***************************************************************************/

/* -------------------------------------------------------------------------- */

// #define TRIMMING
// #define MEDIAN
// #define MEAN
// #define MAX

/* -------------------------------------------------------------------------- */

// #define ORDERING

// #define EPOCH_RESTART
// #define GENERATION_RESTART

// #define COUNTING
// #define MOVING_ONE
// #define MOVING_FRONT

/* -------------------------------------------------------------------------- */

// #define SIMILARITY_CHECK

/* -------------------------------------------------------------------------- */


#ifndef _GLOBAL_H
#define _GLOBAL_H

#include <cassert>
#include <cmath>
#include "myrand.h"
#include "bitwisedistance.h"
#include "spin.h"
#include "nk-wa.h"
#include "doublelinkedlistarray.h"
#include "zkey.h"
#include "sat.h"
#include "maxcut.h"
#include "traps.h"
#include <unordered_map>
#include "color.h"
#include <string>
#define EPSILON (1e-8)
#define INF (1e10)

/* ----------------------------- global settings ---------------------------- */

extern bool USE_HAMMING;
extern bool GHC;
extern bool SELECTION;
extern bool CACHE;
extern bool SHOW_BISECTION;
extern bool SHOW_MP;

extern char outputFilename[100];
extern void gstop ();
extern void outputErrMsg (const char *errMsg);
extern int pow2 (int x);

extern ZKey zKey;
extern MyRand myRand;
extern BitwiseDistance myBD;
extern SPINinstance mySpinGlassParams;
extern SATinstance mySAT;
extern MAXCUTinstance myMAXCUT;
extern NKWAProblem nkwa;
extern USal_NSize_instance my_USal_NSize;

inline int quotientLong(int a) {
    return (a / (sizeof(unsigned long) * 8) );
}

inline int remainderLong(int a) {
    return (a & (sizeof(unsigned long) * 8 - 1));
}

inline double jointEntropy(double p00, double p01, double p10, double p11) {
    double result = 0.0;
    result -= p00 * log(p00);
    result -= p01 * log(p01);
    result -= p10 * log(p10);
    result -= p11 * log(p11);

    return result;
}

inline double mutualInformation(double p00, double p01, double p10, double p11) {
    double result = 0.0;

    double p0x = p00+p01;
    double p1x = p10+p11;
    double px0 = p00+p10;
    double px1 = p01+p11;

    result += (p00 < EPSILON) ? 0.0 : p00 * log (p00 / p0x / px0);
    result += (p01 < EPSILON) ? 0.0 : p01 * log (p01 / p0x / px1);
    result += (p10 < EPSILON) ? 0.0 : p10 * log (p10 / p1x / px0);
    result += (p11 < EPSILON) ? 0.0 : p11 * log (p11 / p1x / px1);

    return result;
}

inline double metric(double p00, double p01, double p10, double p11) {
    return mutualInformation(p00,p01,p10,p11)/jointEntropy(p00,p01,p10,p11);
}

inline double square(double a) {
    return a*a;
}

enum Problem { 
    onemax,
    mktrap,
    ftrap,
    cyctrap,
    nk,
    spin,
    maxsat,
    maxcut,
    USal_NSize,
    USal_NSize_large,
    linear_mktrap,
    exponential_mktrap, 
    four_five_six,
    three_four_five_six_seven,
    ftrap4,
    ftrap6};

enum Verbosity {
    NO,
    COMPACT,
    POPULATION,
    RM
};
extern Verbosity verbose;

static std::unordered_map<std::string, Problem> const Problem_table = { 
    {"onemax", Problem::onemax}, 
    {"mktrap", Problem::mktrap},
    {"ftrap", Problem::ftrap},
    {"cyctrap", Problem::cyctrap},
    {"spin", Problem::spin},
    {"nk", Problem::nk},
    {"maxsat", Problem::maxsat},
    {"maxcut", Problem::maxcut},
    {"USal_NSize", Problem::USal_NSize},
    {"USal_NSize_large", Problem::USal_NSize_large},
    {"linear_mktrap", Problem::linear_mktrap},
    {"exponential_mktrap", Problem::exponential_mktrap},
    {"four_five_six", Problem::four_five_six},
    {"three_four_five_six_seven", Problem::three_four_five_six_seven},
    {"ftrap4", Problem::ftrap4},
    {"ftrap6", Problem::ftrap6}};

#endif
