#ifndef MASSMATRIX_H_INCLUDED
#define MASSMATRIX_H_INCLUDED
#include <iostream>
#include <vector>
#include <cmath>

vector<vector<double>> massmatrix(vector<double> basis, vector<vector<double>> xglobal, vector<double> quad, vector<double> wq, int q, int p){
    //basis - precomputed basis vector for each quad point
    //xglobal - coordinates of the triangle in real space
    //wq - vector of weights
    //quad - coordinates quadrature point 2D

    int prank = (p+1)*(p+2)/2;
    vector<vector<double>> M;
    vector<vector<double>> J;
    vector<vector<double>> sum(prank, vector<double>(prank, 0.0));

    vector<double> cords;
    double detJ;

    if (q == 1) {                                                                                       //if its a linear J constant
        for (int c = 0; c < quad.size(); c +=2) {                                                       //loop through point
            vector<double> grad(basis.begin()+((c/2)*prank),(basis.begin()+((c/2)*prank))+prank);      //compute basis for that point
            for (int i = 0; i < prank; i++) {                                                           //loop though basis
                for (int j = 0; j < prank; j++) {
                sum[i][j] += grad[i]*grad[j]*wq[c/2];
                }
            }
    }
    cords = {quad[0],quad[1]};
    J = jacobian(p,q,xglobal, cords);
    detJ = determinantJ(J);
    M = scalar_matrix_mult(detJ,sum);
    }
    else{                                                                                                   //if its a curve J not constant
    for (int c = 0; c < quad.size(); c+=2) {                                                                //loop through points
        cords = {quad[c], quad[c+1]};
        J = jacobian(p,q,xglobal, cords);                                                                  //computing J at that point for element e
        detJ = determinantJ(J);
        vector<double> grad(basis.begin()+((c/2)*prank),(basis.begin()+((c/2)*prank))+prank);               //compute basis for that point
        for (int i = 0; i < prank; i++){                                                                    //loop though basis
            for (int j = 0; j < prank; j++) {
            sum[i][j] += grad[i]*grad[j]*wq[c/2]*detJ;
            }
        }
    }
    M = sum;
    }
    int r = M.size();
    vector<double> M_V(r*r);
    vector<vector<double>> MI(r, vector<double>(r));
    for (int i =0; i<r; i++){
        for (int j =0; j<r; j++){
            M_V[i*r + j] = M[i][j];
        }
    }
    int gg = invmat(M_V,MI,r);
    if (gg != 0) {cout << "ur fucked";}
    sum.clear(); cords.clear(); J.clear();
return MI;
}
/* Testing mass matrix
vector<vector<double>> x = { {0,0}, {1,0}, {0,1} };
vector<double> xref1 = {0,0};
vector<vector<double>> Mtest = massmatrix(basis_li, x, get<1>(quad_li), get<2>(quad_li), q, p);
//vector<vector<double>> jtest = jacobian(p,ql,x,xref1);
//double Jdet = determinantJ(jtest);
//Mtest = scalar_matrix_mult(Jdet , Mtest);
cout << "mass" <<endl;
for (int i =0; i < Mtest.size(); i++){
    for (int j =0; j <Mtest[0].size(); j++){
        cout << Mtest[i][j] << " ";
    }
    cout<<endl;
}
cout << matrix_sum(Mtest);*/

#endif // MASSMATRIX_H_INCLUDED
