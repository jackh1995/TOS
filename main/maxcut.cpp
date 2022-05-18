#include <iostream>
#include <iomanip>
#include <cstring>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include "maxcut.h"

using namespace std;


double evaluateMAXCUT(char *x, MAXCUTinstance *inst) {
    int temp=1, counter=1;
    double result=0;

	for(vector<int>::iterator it = inst->fvector.begin(); it != inst->fvector.end(); it++) {

		if (counter%3 != 0) temp *= (x[*it-1]*2-1);

		else{
			result += ((*it)*(1-temp)/2);
			temp=1;
		}

        counter++;
    }
	return result;

}

void loadMAXCUT(char *cnf_file_name, char *cnf_opt_file_name, MAXCUTinstance *inst) {
	int temp;
	ifstream input;
	ifstream _input;
	string line;
	_input.open(cnf_opt_file_name);
    getline ( _input, line );
	istringstream _in( line );
	_in >> inst->opt;

	input.open(cnf_file_name);

    getline ( input, line );
	istringstream in( line );
	in >> inst->n;
	in >> inst->m;

	while(getline ( input, line )){
		if (line[0] == '\0') continue;
		istringstream ins( line );
        while ( ins >> temp ) inst->fvector.push_back(temp);
    }
}
