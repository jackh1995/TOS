/***************************************************************************
 *   Copyright (C) 2015 Tian-Li Yu and Shih-Huan Hsu                       *
 *   tianliyu@ntu.edu.tw                                                   *
 ***************************************************************************/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>

#include "statistics.h"
#include "dsmga2.h"
#include "global.h"

#define MAX_GEN 200

int step = 30;

using namespace std;

struct Record {
    int n;
    double nfe;
    double gen;
};

int main (int argc, char *argv[]) {

    if (argc != 5 && argc!=6 && argc !=7 && argc != 8) {
        printf ("sweep ell numConvergence function(0~3) [seed]\n");
        printf ("sweep ell numConvergence 4 [step #] [nk problem #] [seed]\n");
        printf ("sweep ell numConvergence 5 [spin problem #] [seed]\n");
        printf ("sweep ell numConvergence 6 [maxsat problem #] [seed]\n");
        printf ("sweep ell numConvergence 7 [maxcut problem #] [seed]\n");
        printf ("function: \n");
        printf ("     ONEMAX:  0\n");
        printf ("     MK    :  1\n");
        printf ("     FTRAP :  2\n");
        printf ("     CYC   :  3\n");
        printf ("     NK    :  4 (50, 100, 200, 400)      0-based\n");
        printf ("     SPIN  :  5 (36, 100, 196, 400, 784) 1-based\n");
        printf ("     MAXSAT:  6 (20, 50, 75, 100, 200)   1-based\n");
        printf ("     MAXCUT:  7 (50, 75, 100)            1-based\n");
        return -1;
    }

    int ell = atoi (argv[1]);
    int numConvergence = atoi (argv[2]);
    int fffff = atoi(argv[3]);
    int seed = -1;

    int problemNum = 0;
    int neighborNum = 0;
    int stepNum = 0;

    if (fffff == 0 || fffff == 1 || fffff == 2 || fffff == 3) {
        seed = atoi (argv[4]);
    }
    
    if (fffff == 4) {
        neighborNum = 4;
        stepNum = atoi (argv[4]);
        problemNum = atoi (argv[5]);
        seed = atoi (argv[6]);
    }

    if (fffff == 5 || fffff == 6 || fffff == 7) {
        problemNum = atoi (argv[4]);
        seed = atoi (argv[5]);
    }

    int nInitial = 10;

    Statistics st;

    Statistics stGen, stLS, stNFE;

    myRand.seed((unsigned long)seed);

    if (fffff == 5) {
	char filename[200];
        sprintf(filename, "./SPIN/%d/%d_%d",ell, ell, problemNum);
        if (SHOW_BISECTION) printf("Loading: %s\n", filename);
        loadSPIN(filename, &mySpinGlassParams);
    }

    if (fffff == 4) {
        char filename[200];
        sprintf(filename, "./NK_Instance/pnk%d_%d_%d_%d", ell, neighborNum, stepNum, problemNum);
        if (SHOW_BISECTION) printf("Loading: %s\n", filename);
        FILE *fp = fopen(filename, "r");
        loadNKWAProblem(fp, &nkwa);
        fclose(fp);
    }

    if (fffff == 6) {
        char filename[200];
        sprintf(filename, "./SAT/uf%d/uf%d-0%d.cnf",ell,ell,problemNum);
        if (SHOW_BISECTION) printf("Loading: %s\n", filename);
        loadSAT(filename, &mySAT);
    }

    if (fffff == 7) {
        char filename[200];
        char opt_filename[200];
        sprintf(filename, "./maxcut/w05_%d/w05_%d.%d", ell, ell, problemNum);
        sprintf(opt_filename, "./maxcut/g_w05_%d/g_w05_%d.%d", ell, ell, problemNum);
        if (SHOW_BISECTION) printf("Loading: %s\n", filename);
        if (SHOW_BISECTION) printf("Loading groundtruth: %s\n", opt_filename);
        loadMAXCUT(filename, opt_filename, &myMAXCUT);
    }

    bool foundOptima;
    Record rec[3];
    rec[0].n = nInitial;
    rec[1].n = nInitial+step;
    rec[2].n = nInitial+step+step;

    int popu;
    Record q1, q3;

    if (SHOW_BISECTION) printf("Bisection phase 1\n");

    for (int i=0; i<3; ++i) {
        popu = rec[i].n;

        if (SHOW_BISECTION) printf("[%d]: ", popu);
        if (SHOW_MP) {
            printf("[%d/%d] ", seed, popu);
            cout << flush;
        }

        foundOptima = true;

        stGen.reset();
        stNFE.reset();
        stLS.reset();

        for (int j=0; j<numConvergence; j++) {

            DSMGA2 ga(ell, popu, MAX_GEN, -1, fffff);
            ga.doIt();

            stGen.record(ga.getGeneration());
            stNFE.record(Chromosome::nfe);
            stLS.record(Chromosome::lsnfe);


            if (!ga.foundOptima()) {

                foundOptima = false;

                if (SHOW_BISECTION) {
                    printf("-");
                    fflush(NULL);
                }
                break;
            }

            if (SHOW_BISECTION) {
                printf("+");
                fflush(NULL);
            }
        }

        rec[i].gen = stGen.getMean();

        if (!foundOptima)
            rec[i].nfe = INF;
        else
            rec[i].nfe = stNFE.getMean();

        if (SHOW_BISECTION) printf(" : %f \n", rec[i].nfe);

    }

    while (rec[0].nfe < rec[1].nfe  && ((rec[2].n-rec[0].n)*20 > rec[1].n)) {

        rec[2] = rec[1];
        rec[1].n = (rec[0].n + rec[2].n) / 2;
        step /= 2;
        popu = rec[1].n;

        if (SHOW_BISECTION) printf("[%d]: ", popu);
        if (SHOW_MP) {
            printf("[%d/%d] ", seed, popu);
            cout << flush;
        }

        for (int j=0; j<numConvergence; j++) {

            DSMGA2 ga(ell, popu, MAX_GEN, -1, fffff);
            ga.doIt();

            stGen.record(ga.getGeneration());
            stNFE.record(Chromosome::nfe);
            stLS.record(Chromosome::lsnfe);


            if (!ga.foundOptima()) {

                foundOptima = false;

                if (SHOW_BISECTION) {
                    printf("-");
                    fflush(NULL);
                }
                break;
            }

            if (SHOW_BISECTION) {
                printf("+");
                fflush(NULL);
            }
        }


        rec[1].gen = stGen.getMean();

        if (!foundOptima)
            rec[1].nfe = INF;
        else
            rec[1].nfe = stNFE.getMean();

        if (SHOW_BISECTION) printf(" : %f \n", rec[1].nfe);

    }


    while ( (rec[1].nfe >= rec[0].nfe) || (rec[1].nfe >= rec[2].nfe)) {

        popu = rec[2].n + step;

        if (SHOW_BISECTION) printf("[%d]: ", popu);
        if (SHOW_MP) {
            printf("[%d/%d] ", seed, popu);
            cout << flush;
        }

        foundOptima = true;

        stGen.reset();
        stNFE.reset();
        stLS.reset();

        for (int j=0; j<numConvergence; j++) {

            DSMGA2 ga(ell, popu, MAX_GEN, -1, fffff);
            ga.doIt();

            stGen.record(ga.getGeneration());
            stNFE.record(Chromosome::nfe);
            stLS.record(Chromosome::lsnfe);


            if (!ga.foundOptima()) {

                foundOptima = false;

                if (SHOW_BISECTION) {
                    printf("-");
                    fflush(NULL);
                }
                break;
            }

            if (SHOW_BISECTION) {
                printf("+");
                fflush(NULL);
            }
        }


        rec[0] = rec[1];
        rec[1] = rec[2];
        rec[2].n = popu;
        rec[2].gen = stGen.getMean();

        if (!foundOptima)
            rec[2].nfe = INF;
        else
            rec[2].nfe = stNFE.getMean();

        if (SHOW_BISECTION) printf(" : %f \n", rec[2].nfe);

    }


    if (SHOW_BISECTION) printf("Bisection phase 2\n");

    while ( ((rec[2].n-rec[0].n)*20 > rec[1].n) && (rec[2].n>rec[1].n+1) && (rec[1].n>rec[0].n+1)) {

        q1.n = (rec[0].n + rec[1].n) / 2;

        if (SHOW_BISECTION) printf("[%d]: ", q1.n);

        foundOptima = true;

        for (int j=0; j<numConvergence; j++) {

            DSMGA2 ga(ell, q1.n, MAX_GEN, -1, fffff);
            ga.doIt();

            if (!ga.foundOptima()) {
                foundOptima = false;
                if (SHOW_BISECTION) {
                    printf("-");
                    fflush(NULL);
                }
                break;
            }
            if (SHOW_BISECTION) {
                printf("+");
                fflush(NULL);
            }
            if (j==0) {
                stGen.reset();
                stLS.reset();
                stNFE.reset();
            }
            stGen.record(ga.getGeneration());
            stNFE.record(Chromosome::nfe);
            stLS.record(Chromosome::lsnfe);
        }

        q1.gen = stGen.getMean();
        if (foundOptima)
            q1.nfe = stNFE.getMean();
        else
            q1.nfe = INF;


        if (SHOW_BISECTION) printf(" : %f \n", q1.nfe);


        q3.n = (rec[1].n + rec[2].n) / 2;

        if (SHOW_BISECTION) printf("[%d]: ", q3.n);

        foundOptima = true;

        for (int j=0; j<numConvergence; j++) {

            DSMGA2 ga(ell, q3.n, MAX_GEN, -1, fffff);
            ga.doIt();

            if (!ga.foundOptima()) {
                foundOptima = false;
                if (SHOW_BISECTION) {
                    printf("-");
                    fflush(NULL);
                }
                break;
            }
            if (SHOW_BISECTION) {
                printf("+");
                fflush(NULL);
            }
            if (j==0) {
                stGen.reset();
                stLS.reset();
                stNFE.reset();
            }
            stGen.record(ga.getGeneration());
            stNFE.record(Chromosome::nfe);
            stLS.record(Chromosome::lsnfe);
        }

        q3.gen = stGen.getMean();
        if (foundOptima) {
            q3.nfe = stNFE.getMean();
        } else
            q3.nfe = INF;

        if (SHOW_BISECTION) printf(" : %f \n", q3.nfe);

        if (rec[1].nfe < q1.nfe && rec[1].nfe < q3.nfe) {
            rec[0] = q1;
            rec[2] = q3;
        } else if (q1.nfe < rec[1].nfe && q1.nfe < q3.nfe) {
            rec[2] = rec[1];
            rec[1] = q1;
        } else { // q3nfe smallest
            rec[0] = rec[1];
            rec[1] = q3;
        }
    };

    if (fffff == 4)
        freeNKWAProblem(&nkwa);

    /* -------------------------------------------------------------------------- */
    /*                                 cout result                                */
    /* -------------------------------------------------------------------------- */
    
    
    cout << endl;
    cout << "--------------------------------------------------------------------------" << endl;
    cout << "seed " << seed << " done ";
    printf("pop=[%d] gen=[%f] nfe=[%f]\n", rec[1].n, rec[1].gen, rec[1].nfe);
    cout << "--------------------------------------------------------------------------" << endl;
    // printf("population: %d\n", rec[1].n);
    // printf("generation: %f\n", rec[1].gen);
    // printf("NFE: %f\n", rec[1].nfe);

    /* -------------------------------------------------------------------------- */
    /*                               Write to files                               */
    /* -------------------------------------------------------------------------- */
    
    std::ofstream outfile;
    string filename;
    
    string mode;
    if ((USE_ADAMEAN) && (USE_MOILS) && (USE_HAMMING)) {
        mode = "ttt";
    } else if ((!USE_ADAMEAN) && (USE_MOILS) && (USE_HAMMING)) {
        mode = "ftt";
    } else if ((USE_ADAMEAN) && (!USE_MOILS) && (USE_HAMMING)) {
        mode = "tft";
    } else if ((USE_ADAMEAN) && (USE_MOILS) && (!USE_HAMMING)) {
        mode = "ttf";
    } else if ((!USE_ADAMEAN) && (!USE_MOILS) && (USE_HAMMING)) {
        mode = "fft";
    } else if ((USE_ADAMEAN) && (!USE_MOILS) && (!USE_HAMMING)) {
        mode = "tff";
    } else if ((!USE_ADAMEAN) && (USE_MOILS) && (!USE_HAMMING)) {
        mode = "ftf";
    } else {
        mode = "fff";
    }
    
    mode = "-" + mode;

    switch(fffff) {
        case 0:
            filename = "./experiment/results/onemax" + std::to_string(ell) + mode + ".csv";
            break;
        case 1:
            filename = "./experiment/results/mktrap" + std::to_string(ell) + mode + ".csv";
            break;
        case 2:
            filename = "./experiment/results/ftrap" + std::to_string(ell) + mode + ".csv";
            break;
        case 3:
            filename = "./experiment/results/cyctrap" + std::to_string(ell) + mode + ".csv";
            break;
        case 4:
            filename = "./experiment/results/nks1-" + std::to_string(ell) + mode + ".csv";
            break;
        case 5:
            filename = "./experiment/results/spin" + std::to_string(ell) + mode + ".csv";
            break;
        case 6:
            filename = "./experiment/results/maxsat" + std::to_string(ell) + mode + ".csv";
            break;
        case 7:
            filename = "./experiment/results/maxcut" + std::to_string(ell) + mode + ".csv";
            break;
        default:
            cerr << "Problem not defined";
            return -1;
    }

    // cout << "Writing result to file " + filename << endl;
    outfile.open(filename, std::ios_base::app);
    outfile << rec[1].n << ',';
    outfile << rec[1].gen << ',';
    outfile << rec[1].nfe << '\n';
    outfile.close();

    return EXIT_SUCCESS;
}