#ifndef RETURNCORD_H_INCLUDED
#define RETURNCORD_H_INCLUDED
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

vector<vector<double>> elmcord(vector<vector<double>> C, vector<double> nodes){
vector<vector<double>> cord;

for (int i = 0; i <nodes.size(); i++){
    double num = nodes[i]-1;
    vector<double> temp = C[num];
    cord.push_back(temp);
    temp.clear();
}

return cord;
}


#endif // RETURNCORD_H_INCLUDED
