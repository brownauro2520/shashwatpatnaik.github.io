#ifndef NORMALC_H_INCLUDED
#define NORMALC_H_INCLUDED
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

vector<double> idk(int side){
    vector<double> wtf;
    if(side == 1){
        wtf = {1,0};
    }
    else if(side == 2){
        wtf = {-1,1};
    }
    else{
        wtf = {0,-1};
    }
return wtf;
}


vector<vector<double>> normalC(int p, int q, int side, vector<double> quad, vector<vector<double>> cord) {
    vector<double> xref;
    vector<vector<double>> N;
    vector<double> xcor = idk(side);
    for (int i = 0; i < quad.size(); i+=1){
        xref = reflagrange(side,quad[i]);
        vector<vector<double>> J = jacobian(p, q, cord, xref);
        vector<double> dxeta = {J[0][0], J[1][0]};
        vector<double> dyneta = {J[0][1], J[1][1]};
        vector<double> temp1 = multiplyVS(xcor[0],dxeta);
        vector<double> temp2 = multiplyVS(xcor[1],dyneta);
        vector<double> Ntemp = add_v(temp1,temp2);
        N.push_back({Ntemp[1],-Ntemp[0]});
        Ntemp.clear(); temp1.clear(); temp2.clear();
    }
return N;
}


#endif // NORMALC_H_INCLUDED
