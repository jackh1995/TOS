/***************************************************************************
 *   Copyright (C) 2015 Tian-Li Yu and Shih-Huan Hsu                       *
 *   tianliyu@ntu.edu.tw                                                   *
 ***************************************************************************/

#include <vector>
#include <vector>
#include <algorithm>
#include <iterator>

#include <iostream>
#include "chromosome.h"
#include "dsmga2.h"
#include "fastcounting.h"
#include "statistics.h"
#include <numeric>
#include <iomanip>
using namespace std;


DSMGA2::DSMGA2 (int n_ell, int n_nInitial, int n_maxGen, int n_maxFe, int fffff) {
    previousFitnessMean = -INF;
    ell = n_ell;
    nCurrent = (n_nInitial/2)*2;  // has to be even
    masks = new vector<int>[ell];
    orderRM = new int[nCurrent];
    Chromosome::function = (Chromosome::Function)fffff;
    Chromosome::nfe = 0;
    Chromosome::lsnfe = 0;
    Chromosome::hitnfe = 0;
    Chromosome::hit = false;

    selectionPressure = 2;
    maxGen = n_maxGen;
    maxFe = n_maxFe;

    graph.init(ell);
    graph_size.init(ell);
     
    bestIndex = -1;
    // masks = new vector<int>[ell];
    selectionIndex = new int[nCurrent];
    orderN = new int[nCurrent];
    orderELL = new int[ell];
    population = new Chromosome[nCurrent];
    fastCounting = new FastCounting[ell];

    for (int i = 0; i < ell; i++)
        fastCounting[i].init(nCurrent);

    #ifdef TRIMMING
    rm_success_sizes_mean = ell/2;
    #endif

    pHash.clear();
    for (int i=0; i<nCurrent; ++i) {
        population[i].initR(ell);
        double f = population[i].getFitness();
        pHash[population[i].getKey()] = f;
    }

    #ifdef ORDERING
    orderILS = vector<int>(ell/2);
    iota(orderILS.begin(), orderILS.end(), 0);
    #endif
    
    if (GHC) {
        for (int i=0; i < nCurrent; i++)
            population[i].GHC();
    }
}


DSMGA2::~DSMGA2 () {
    delete []masks;
    delete []orderN;
    delete []orderELL;
    delete []selectionIndex;
    delete []population;
    delete []fastCounting;
}



bool DSMGA2::isSteadyState () {

    if (stFitness.getNumber () <= 0)
        return false;

    if (previousFitnessMean < stFitness.getMean ()) {
        previousFitnessMean = stFitness.getMean () + EPSILON;
        return false;
    }

    return true;
}

bool DSMGA2::converged() {
    if (stFitness.getMax() == lastMax &&
        stFitness.getMean() == lastMean &&
        stFitness.getMin() == lastMin)
        convergeCount++;
    else
        convergeCount = 0;

    lastMax = stFitness.getMax();
    lastMean = stFitness.getMean();
    lastMin = stFitness.getMin();
    return (convergeCount > 300) ? true : false;
}

int DSMGA2::doIt (bool output) {
    generation = 0;
    while (!shouldTerminate ()) {
        oneRun (output);
       #ifdef DEBUG
       cin.get();
        #endif
    }
    return generation;
}


void DSMGA2::oneRun (bool output) {

    if (CACHE)
        Chromosome::cache.clear();

    mixing();

    #ifdef ORDERING
    reset_orderILS();
    #endif

    double max = -INF;
    stFitness.reset ();

    for (int i = 0; i < nCurrent; ++i) {
        double fitness = population[i].getFitness();
        if (fitness > max) {
            max = fitness;
            bestIndex = i;
        }
        stFitness.record (fitness);

    }

    if (output)
        showStatistics ();

    ++generation;
}


bool DSMGA2::shouldTerminate () {
    bool  termination = false;

    if (maxFe != -1) {
        if (Chromosome::nfe > maxFe)
            termination = true;
    }

    if (maxGen != -1) {
        if (generation > maxGen)
            termination = true;
    }

    if (population[0].getMaxFitness() <= stFitness.getMax() )
        termination = true;


    if (stFitness.getMax() - EPSILON <= stFitness.getMean() )
        termination = true;

    return termination;

}


bool DSMGA2::foundOptima () {
    return (stFitness.getMax() > population[0].getMaxFitness());
}


void DSMGA2::showStatistics () {

    printf ("Gen:%d  Fitness:(Max/Mean/Min):%f/%f/%f \n ",
            generation, stFitness.getMax (), stFitness.getMean (),
            stFitness.getMin ());
    fflush(NULL);
}



void DSMGA2::buildFastCounting() {

    if (SELECTION) {
        for (int i = 0; i < nCurrent; i++)
            for (int j = 0; j < ell; j++) {
                fastCounting[j].setVal(i, population[selectionIndex[i]].getVal(j));
            }

    } else {
        for (int i = 0; i < nCurrent; i++)
            for (int j = 0; j < ell; j++) {
                fastCounting[j].setVal(i, population[i].getVal(j));
            }
    }

}

int DSMGA2::countOne(int x) const {

    int n = 0;

    for (int i=0; i<fastCounting[0].lengthLong; ++i) {
        unsigned long val = 0;

        val = fastCounting[x].gene[i];

        n += myBD.countOne(val);
    }

    return n;
}


int DSMGA2::countXOR(int x, int y) const {

    int n = 0;

    for (int i=0; i<fastCounting[0].lengthLong; ++i) {
        unsigned long val = 0;


        val = fastCounting[x].gene[i];

        val ^= fastCounting[y].gene[i];

        n += myBD.countOne(val);
    }

    return n;
}

//2016-03-09
// Almost identical to DSMGA2::findClique
// except check 00 or 01 before adding connection
void DSMGA2::findMask(Chromosome& ch, vector<int>& result,int startNode){
    result.clear();

    
	DLLA rest(ell);
	genOrderELL();
	for( int i = 0; i < ell; i++){
		if(orderELL[i] == startNode)
			result.push_back(orderELL[i]);
		else
			rest.insert(orderELL[i]);
	}

	double *connection = new double[ell];

	for(DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++){
	    pair<double, double> p = graph(startNode, *iter);
		int i = ch.getVal(startNode);
		int j = ch.getVal(*iter);
		if(i == j)//p00 or p11
			connection[*iter] = p.first;
		else      //p01 or p10
			connection[*iter] = p.second;
	}
   
    while(!rest.isEmpty()){

	    double max = -INF;
		int index = -1;
		for(DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++){
		    if(max < connection[*iter]){
			    max = connection[*iter];
				index = *iter;
			}
		}

		rest.erase(index);
		result.push_back(index);

		for(DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++){
			pair<double, double> p = graph(index, *iter);
			int i = ch.getVal(index);
			int j = ch.getVal(*iter);
			if(i == j)//p00 or p11
				connection[*iter] += p.first;
			else      //p01 or p10
				connection[*iter] += p.second;
		}
	}

	delete []connection;
  
}

void DSMGA2::restrictedMixing(Chromosome& ch) {
    
    int startNode = myRand.uniformInt(0, ell - 1);    

    vector<int> mask;
	findMask(ch, mask,startNode);
    size_t size = findSize(ch, mask);
   
    vector<int> mask_size; 
    findMask_size(ch,mask_size,startNode,size);
    size_t size_original = findSize(ch,mask_size);

    if (size > size_original)
        size = size_original;
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

    while (mask.size() > size)
        mask.pop_back();
   
    bool taken = restrictedMixing(ch, mask);


    EQ = true;
    if (taken) {
    
        genOrderN();

        for (int i=0; i<nCurrent; ++i) {

            #ifdef SIMILARITY_CHECK
            if (get_ch_dist(ch, population[orderN[i]]) <= ell/2) {
                if (EQ) {
                    backMixingE(ch, mask, population[orderN[i]]);
                } else {
                    backMixing(ch, mask, population[orderN[i]]);
                }
            }
            #else
            if (EQ)
                backMixingE(ch, mask, population[orderN[i]]);
            else
                backMixing(ch, mask, population[orderN[i]]);
            #endif
        }
    }

}
void DSMGA2::findMask_size(Chromosome& ch, vector<int>& result,int startNode,int bound){
    result.clear();

    
	DLLA rest(ell);

	for( int i = 0; i < ell; i++){
		if(orderELL[i] == startNode)
			result.push_back(orderELL[i]);
		else
			rest.insert(orderELL[i]);
	}

	double *connection = new double[ell];

	for(DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++){
	    pair<double, double> p = graph_size(startNode, *iter);
		int i = ch.getVal(startNode);
		int j = ch.getVal(*iter);
		if(i == j)//p00 or p11
			connection[*iter] = p.first;
		else      //p01 or p10
			connection[*iter] = p.second;
	}
    bound--;
    while(!rest.isEmpty()&&bound>0){
        bound--;
	    double max = -INF;
		int index = -1;
		for(DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++){
		    if(max < connection[*iter]){
			    max = connection[*iter];
				index = *iter;
			}
		}

		rest.erase(index);
		result.push_back(index);

		for(DLLA::iterator iter = rest.begin(); iter != rest.end(); iter++){
			pair<double, double> p = graph_size(index, *iter);
			int i = ch.getVal(index);
			int j = ch.getVal(*iter);
			if(i == j)//p00 or p11
				connection[*iter] += p.first;
			else      //p01 or p10
				connection[*iter] += p.second;
		}
	}

    delete []connection;
  
}

void DSMGA2::backMixing(Chromosome& source, vector<int>& mask, Chromosome& des) {

    Chromosome trial(ell);
    trial = des;
    for (vector<int>::iterator it = mask.begin(); it != mask.end(); ++it)
        trial.setVal(*it, source.getVal(*it));

    if (trial.getFitness() > des.getFitness()) {
        pHash.erase(des.getKey());
        pHash[trial.getKey()] = trial.getFitness();
        des = trial;
          
        return;
    }

}

void DSMGA2::backMixingE(Chromosome& source, vector<int>& mask, Chromosome& des) {

    Chromosome trial(ell);
    trial = des;
    for (vector<int>::iterator it = mask.begin(); it != mask.end(); ++it)
        trial.setVal(*it, source.getVal(*it));

    if (trial.getFitness() > des.getFitness()) {
        pHash.erase(des.getKey());
        pHash[trial.getKey()] = trial.getFitness();

        EQ = false;
        des = trial;

return;
    }

    //2016-10-21
    if (trial.getFitness() >= des.getFitness() - EPSILON) {
        pHash.erase(des.getKey());
        pHash[trial.getKey()] = trial.getFitness();

        des = trial;
        return;
    }

}

bool DSMGA2::restrictedMixing(Chromosome& ch, vector<int>& mask) {

    // bool taken = false;
    // size_t lastUB = 0;

    // for (size_t ub = 1; ub <= mask.size(); ++ub) {

    //     size_t size = 1;
    //     Chromosome trial(ell);
    //     trial = ch;
	    
	// 	//2016-03-03
	//     vector<int> takenMask;

    //     for (vector<int>::iterator it = mask.begin(); it != mask.end(); ++it) {
            
    //         //2016-03-03
	// 		takenMask.push_back(*it);

    //         trial.flip(*it);

    //         ++size;
    //         if (size > ub) break;
    //     }

    //     //if (isInP(trial)) continue;
    //     //2016-10-21
    //     if (isInP(trial)) break;

    //     if (trial.getFitness() >= ch.getFitness() - EPSILON) {
    //         pHash.erase(ch.getKey());
    //         pHash[trial.getKey()] = trial.getFitness();

    //         taken = true;
    //         ch = trial;
    //     }

    //     if (taken) {
    //         lastUB = ub;
    //         break;
    //     }
    // }

    bool taken(false);
    size_t lastUB(0);
    Chromosome trial(ell);
    trial = ch;
    
    // the mask has passed the donor check
    for (size_t idx = 0; idx != mask.size(); ++idx) {
        
        size_t ub;
        #ifdef ORDERING
        ub = orderILS[idx]+1;
        if (ub > mask.size()) {
            continue;
        }
        #else
        ub = idx + 1;
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
        if (trial.getFitness() >= ch.getFitness() - EPSILON) {
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

    if (lastUB != 0) {
        while (mask.size() > lastUB)
            mask.pop_back();
    }

    return taken;

}

size_t DSMGA2::findSize(Chromosome& ch, vector<int>& mask) {

    DLLA candidate(nCurrent);
    for (int i=0; i<nCurrent; ++i)
        candidate.insert(i);

    size_t size = 0;
    for (vector<int>::iterator it = mask.begin(); it != mask.end(); ++it) {

        int allele = ch.getVal(*it);

        for (DLLA::iterator it2 = candidate.begin(); it2 != candidate.end(); ++it2) {
            #ifdef SIMILARITY_CHECK
            if (get_ch_dist(ch, population[*it2]) > ell/2) {
                candidate.erase(*it2);
            }
            #endif
            
            if (population[*it2].getVal(*it) == allele)
                candidate.erase(*it2);

            if (candidate.isEmpty())
                break;
        }

        if (candidate.isEmpty())
            break;

        ++size;
    }

    return size;


}

size_t DSMGA2::findSize(Chromosome& ch, vector<int>& mask, Chromosome& ch2) {

    size_t size = 0;
    for (vector<int>::iterator it = mask.begin(); it != mask.end(); ++it) {
        if (ch.getVal(*it) == ch2.getVal(*it)) break;
        ++size;
    }
    return size;
}

void DSMGA2::mixing() {

    if (SELECTION)
        selection();

    //* really learn model
    buildFastCounting();
    buildGraph();
    buildGraph_sizecheck();
    //for (int i=0; i<ell; ++i)
    //    findClique(i, masks[i]); // replaced by findMask in restrictedMixing

    int repeat = (ell>50)? ell/50: 1;

    for (int k = 0; k < repeat; ++k) {
        genOrderRM();
        #ifdef TRIMMING
        for (int i = 0; i < nCurrent/2; ++i) {
            restrictedMixing(population[orderRM[i]]);
            if (Chromosome::hit)
                break;
        }
        if (Chromosome::hit)
            break;
        if (!rm_success_sizes.empty()){
            size_t sum(0);
            vector<int>::iterator it;
            for (it = rm_success_sizes.begin(); it != rm_success_sizes.end(); ++it){
                sum += (*it);
            }
            rm_success_sizes_mean = sum / rm_success_sizes.size() + 1;
        }
        for (int i=nCurrent/2; i<nCurrent; ++i) {
            restrictedMixing(population[orderRM[i]]);
            if (Chromosome::hit) 
                break;
        }
        if (Chromosome::hit)
            break;
        rm_success_sizes.clear();
        #else
        for (int i=0; i<nCurrent; ++i) {
            restrictedMixing(population[orderRM[i]]);
            if (Chromosome::hit) break;
        }
        if (Chromosome::hit) break;
        #endif
    }
}

inline bool DSMGA2::isInP(const Chromosome& ch) const {

    unordered_map<unsigned long, double>::const_iterator it = pHash.find(ch.getKey());
    return (it != pHash.end());
}

inline void DSMGA2::genOrderN() {
    myRand.uniformArray(orderN, nCurrent, 0, nCurrent-1);
}

inline void DSMGA2::genOrderELL() {
    myRand.uniformArray(orderELL, ell, 0, ell-1);
}

void DSMGA2::buildGraph() {

    int *one = new int [ell];
    for (int i=0; i<ell; ++i) {
        one[i] = countOne(i);
    }

    for (int i=0; i<ell; ++i) {

        for (int j=i+1; j<ell; ++j) {

            int n00, n01, n10, n11;
            int nX =  countXOR(i, j);

            n11 = (one[i]+one[j]-nX)/2;
            n10 = one[i] - n11;
            n01 = one[j] - n11;
            n00 = nCurrent - n01 - n10 - n11;

            double p00 = (double)n00/(double)nCurrent;
            double p01 = (double)n01/(double)nCurrent;
            double p10 = (double)n10/(double)nCurrent;
            double p11 = (double)n11/(double)nCurrent;
            double p1_ = p10 + p11;
            double p0_ = p00 + p01;
            double p_0 = p00 + p10;
            double p_1 = p01 + p11;

            double linkage = computeMI(p00,p01,p10,p11);
            
            //2016-04-08_computeMI_entropy
            double linkage00 = 0.0, linkage01 = 0.0;
            if (p00 > EPSILON)
                linkage00 += p00*log(p00/p_0/p0_);
            if (p11 > EPSILON)
                linkage00 += p11*log(p11/p_1/p1_);
            if (p01 > EPSILON)
                linkage01 += p01*log(p01/p0_/p_1);
            if (p10 > EPSILON)
                linkage01 += p10*log(p10/p1_/p_0);
           
            if(Chromosome::nfe < 0){
                pair<double, double> p(linkage, linkage);
                graph.write(i, j, p);
            }
            else{
                pair<double, double> p(linkage00, linkage01);
                graph.write(i, j, p);
            }
        }
    }


    delete []one;

}
void DSMGA2::buildGraph_sizecheck() {

    int *one = new int [ell];
    for (int i=0; i<ell; ++i) {
        one[i] = countOne(i);
    }

    for (int i=0; i<ell; ++i) {

        for (int j=i+1; j<ell; ++j) {

            int n00, n01, n10, n11;
            int nX =  countXOR(i, j);

            n11 = (one[i]+one[j]-nX)/2;
            n10 = one[i] - n11;
            n01 = one[j] - n11;
            n00 = nCurrent - n01 - n10 - n11;

            double p00 = (double)n00/(double)nCurrent;
            double p01 = (double)n01/(double)nCurrent;
            double p10 = (double)n10/(double)nCurrent;
            double p11 = (double)n11/(double)nCurrent;
            double p1_ = p10 + p11;
            double p0_ = p00 + p01;
            double p_0 = p00 + p10;
            double p_1 = p01 + p11;

            double linkage = computeMI(p00,p01,p10,p11);
            
            //2016-04-08_computeMI_entropy
            double linkage00 = 0.0, linkage01 = 0.0;
            if (p00 > EPSILON)
                linkage00 += p00*log(p00/p_0/p0_);
            if (p11 > EPSILON)
                linkage00 += p11*log(p11/p_1/p1_);
            if (p01 > EPSILON)
                linkage01 += p01*log(p01/p0_/p_1);
            if (p10 > EPSILON)
                linkage01 += p10*log(p10/p1_/p_0);
        
	
            pair<double, double> p(linkage, linkage);
            graph_size.write(i, j, p);
			
        }
    }


    delete []one;

}


// from 1 to ell, pick by max edge
void DSMGA2::findClique(int startNode, vector<int>& result) {
    
   }

    
    double DSMGA2::computeMI(double a00, double a01, double a10, double a11) const {

    double p0 = a00 + a01;
    double q0 = a00 + a10;
    double p1 = 1-p0;
    double q1 = 1-q0;

    double join = 0.0;
    if (a00 > EPSILON)
        join += a00*log(a00);
    if (a01 > EPSILON)
        join += a01*log(a01);
    if (a10 > EPSILON)
        join += a10*log(a10);
    if (a11 > EPSILON)
        join += a11*log(a11);

    double p = 0.0;
    if (p0 > EPSILON)
        p += p0*log(p0);
    if (p1 > EPSILON)
        p += p1*log(p1);


    double q = 0.0;
    if (q0 > EPSILON)
        q += q0*log(q0);
    if (q1 > EPSILON)
        q += q1*log(q1);

    return -p-q+join;

}


void DSMGA2::selection () {
    tournamentSelection ();
}


// tournamentSelection without replacement
void DSMGA2::tournamentSelection () {
    int i, j;

    int randArray[selectionPressure * nCurrent];

    for (i = 0; i < selectionPressure; i++)
        myRand.uniformArray (randArray + (i * nCurrent), nCurrent, 0, nCurrent - 1);

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

void DSMGA2::genOrderRM() {
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

#ifdef SIMILARITY_CHECK
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
#endif