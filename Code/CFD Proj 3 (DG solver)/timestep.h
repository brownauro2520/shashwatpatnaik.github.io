#ifndef MASSINVMULT_H_INCLUDED
#define MASSINVMULT_H_INCLUDED
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>
using namespace std;
vector<double> massinvmult(vector<vector<double>> iM,vector<vector<double>> R, int p, int E_num){
int Np = (p+1)*(p+2)/2;
int Ndof = E_num*Np*4;
vector<double> F(Ndof,0.0);
for (int ie = 0; ie < E_num; ie++){
for (int i = 0; i<Np; i++){
    for (int j = 0; j<Np; j++) {
        F[ie*4*Np + i*4 + 0] = F[ie*4*Np + i*4 + 0] - (iM[ie*Np + i][j]*R[ie*Np + j][0]);
        F[ie*4*Np + i*4 + 1] = F[ie*4*Np + i*4 + 1] - (iM[ie*Np + i][j]*R[ie*Np + j][1]);
        F[ie*4*Np + i*4 + 2] = F[ie*4*Np + i*4 + 2] - (iM[ie*Np + i][j]*R[ie*Np + j][2]);
        F[ie*4*Np + i*4 + 3] = F[ie*4*Np + i*4 + 3] - (iM[ie*Np + i][j]*R[ie*Np + j][3]);
    }
}
}
return F;
}
#endif // MASSINVMULT_H_INCLUDED