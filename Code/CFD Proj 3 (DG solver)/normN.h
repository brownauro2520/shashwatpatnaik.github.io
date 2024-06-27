#ifndef NORMN_H_INCLUDED
#define NORMN_H_INCLUDED
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

vector<double> normN(vector<vector<double>> N){
    vector<double> norm(N.size());

    for (int i =0; i < N.size(); i++){
        norm[i] = sqrt(pow(N[i][0],2)+pow(N[i][1],2));
    }

return norm;
}


#endif // NORMN_H_INCLUDED
