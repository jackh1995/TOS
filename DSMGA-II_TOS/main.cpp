#include <math.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <time.h>
#include "statistics.h"
#include "dsmga2.h"
#include "global.h"
#include "chromosome.h"
using namespace std;

int main (int argc, char *argv[]) {
    
    if (argc != 9) {
        printf ("DSMGA2 ell nInitial function maxGen maxFe repeat display rand_seed\n");
        printf ("Functions: \n");
        printf ("    0. ONEMAX: Any ell\n");
        printf ("    1. MK    : Any ell\n");
        printf ("    2. FTRAP : Any ell\n");
        printf ("    3. CYC   : Any ell\n");
        printf ("    4. NK    : 25, 50, 100, 200, 400\n");
        printf ("    5. SPIN  : 36, 100, 196, 400, 784\n");
        printf ("    6. SAT   : 20, 50, 75, 100, 200\n");
        printf ("    7. MAXCUT: 50, 75, 100\n");
        printf ("    8. USal_NSize: 100, 200, 300, 400\n");
        printf ("    9. USal_NSize_large: 100, 200, 300, 400\n");
        printf ("    10. Linear MKTRAP\n");
        printf ("    11. Exponential MKTRAP\n");
        return -1;
    }

    int ell = atoi (argv[1]); // problem size
    int nInitial = atoi (argv[2]); // initial population size
    int fffff = atoi (argv[3]); // function
    int maxGen = atoi (argv[4]); // max generation
    int maxFe = atoi (argv[5]); // max fe
    int repeat = atoi (argv[6]); // how many time to repeat
    verbose = Verbosity(atoi (argv[7])); // display each generation or not
    int rand_seed = atoi (argv[8]);  // rand seed

    char* inst_num_env = getenv("DSMGA2_INSTANCE_NUMBER");
    int inst_num = 1;
    if(inst_num_env){ // overriden by environment variable
        inst_num = atoi(inst_num_env);
    }
    
    switch(fffff) {

        char filename[200];
        char opt_filename[200];
        FILE *fp;
        case onemax: case mktrap: case cyctrap: case ftrap: case linear_mktrap: case exponential_mktrap:
            break;
        case nk: 
            sprintf(filename, "./NK_Instance/pnk%d_%d_%d_%d", ell, 4, 1, inst_num);
            if (SHOW_BISECTION) printf("Loading: %s\n", filename);
            fp = fopen(filename, "r");
            loadNKWAProblem(fp, &nkwa);
            break;
        case spin:
            sprintf(filename, "./SPIN/%d/%d_%d",ell, ell, inst_num);
            if (SHOW_BISECTION) printf("Loading: %s\n", filename);
            loadSPIN(filename, &mySpinGlassParams); 
            break;
        case maxsat:
            sprintf(filename, "./SAT/uf%d/uf%d-0%d.cnf", ell, ell, inst_num);
            if (SHOW_BISECTION) printf("Loading: %s\n", filename);
            loadSAT(filename, &mySAT);
            break;
        case maxcut:
            sprintf(filename, "./maxcut/w05_%d/w05_%d.%d", ell, ell, inst_num);
            sprintf(opt_filename, "./maxcut/g_w05_%d/g_w05_%d.%d", ell, ell, inst_num);
            if (SHOW_BISECTION) printf("Loading: %s\n", filename);
            if (SHOW_BISECTION) printf("Loading groundtruth: %s\n", opt_filename);
            loadMAXCUT(filename, opt_filename, &myMAXCUT);
            break;
        case USal_NSize:
            sprintf(filename, "./traps/USal/NSize/%d_3_7_%d", ell, inst_num);
            if (SHOW_BISECTION) printf("Loading: %s\n", filename);
            load_USal_NSize(filename, &my_USal_NSize);
            break;
        case USal_NSize_large:
            sprintf(filename, "./traps/USal/NSize/%d_3_10_%d", ell, inst_num);
            if (SHOW_BISECTION) printf("Loading: %s\n", filename);
            load_USal_NSize(filename, &my_USal_NSize);
            break;
        default:
            throw runtime_error("Error: Problem not defined");
    }
    
    // set random seed
    if (rand_seed != -1)
        myRand.seed((unsigned long)rand_seed);
    
    // start evolving
    Statistics stGen, stFE, stINITFE, stLSFE, stRMFE, stBMFE;
    int usedGen;
    int failNum = 0;
    for (int i = 0; i < repeat; ++i) {
        DSMGA2 ga (ell, nInitial, maxGen, maxFe, fffff);
        usedGen = ga.doIt();
        if (!ga.foundOptima()) {
            failNum++;
            // if (repeat > 1)
            //     printf ("-");
            printf ("-");
        } else {
            stGen.record (usedGen);
            stFE.record (Chromosome::nfe);
            stINITFE.record(Chromosome::initnfe);
            stLSFE.record (Chromosome::lsnfe);
            stRMFE.record (Chromosome::rmnfe);
            stBMFE.record (Chromosome::bmnfe);
            // if (repeat > 1) {
            //     printf ("+");
            // }
            printf ("+");
        }
        
        // print the best seen chromosome if repeat == 1
        // if (repeat == 1) {
        //     cout << endl;
        //     cout << "SUMMARY: ";
        //     if (ga.foundOptima())
        //         cout << GREEN << "SUCCESS" << RESET << endl;
        //     else
        //         cout << RED << "FAIL" << RESET << endl;
        //     cout << ga.getPopulation()[ga.getBestIndex()];
        // }
        fflush (NULL);
    }
    
    // cout << fixed << setprecision(2) << endl;
    // cout << "GEN=[" << BOLDCYAN << stGen.getMean() << RESET << "] ";
    // cout << "NFE=[" << BOLDCYAN << stFE.getMean() << RESET << "] ";
    // cout << "INITNFE=[" << BOLDCYAN << stINITFE.getMean() << RESET << "] ";
    // cout << "LSNFE=[" << BOLDCYAN << stLSFE.getMean() << RESET << "] ";
    // cout << "RMNFE=[" << BOLDCYAN << stRMFE.getMean() << RESET << "] ";
    // cout << "BMNFE=[" << BOLDCYAN << stBMFE.getMean() << RESET << "] ";
    // cout << "SUCCESS=[" << BOLDCYAN << repeat-failNum << '/' << repeat << RESET << "]\n";

    cout << fixed << setprecision(2) << endl;
    cout << "GEN=[" << stGen.getMean() << "] ";
    cout << "NFE=[" << stFE.getMean() << "] ";
    cout << "INITNFE=[" << stINITFE.getMean() << "] ";
    cout << "LSNFE=[" << stLSFE.getMean() << "] ";
    cout << "RMNFE=[" << stRMFE.getMean() << "] ";
    cout << "BMNFE=[" << stBMFE.getMean() << "] ";
    cout << "SUCCESS=[" << repeat-failNum << '/' << repeat << "]\n";
    
    if (fffff == 4) freeNKWAProblem(&nkwa);
    return EXIT_SUCCESS;
}