#ifndef JACOBIAN_H_INCLUDED
#define JACOBIAN_H_INCLUDED
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

vector<vector<double>> jacobian(int p, int q, vector<vector<double>> &x, vector<double> &xref) {

vector<vector<double>> J(2, vector<double>(2, 0.0));
vector<double> grad = gradientL(xref, q);
int prank = (q+1)*(q+2)/2;
vector<vector<double>> Jtemp(2, vector<double>(2));
vector<double> tempgrad;


for (int i = 0; i < prank; i++){
    tempgrad = {grad[i], grad[prank+i]};
    Jtemp = multiplyVV(x[i],tempgrad);
    J = addMatrices(J,Jtemp);
}

return J;
}

/* testing

vector<vector<double>> x = {{-1.000000000000000E+02,-1.000000000000000E+02}, {-7.142857142859999E+01,-1.000000000000000E+02}, {-7.822000646900000E+01, -7.998608171300000E+01}};
vector<double> xref1 = {0,0};
vector<vector<double>> J = jacobian(q, q, x, xref1);
double q1 = x[1][0]-x[0][0];
double q2 = x[2][0]-x[0][0];
double q3 = x[1][1]-x[0][1];
double q4 = x[2][1]-x[0][1];
vector<vector<double>> Jtest = { {q1 , q2} , {q3 , q4} };
cout << "func" <<endl;
for (int i = 0; i < J.size(); i++) {
        for (int j = 0; j < J[0].size(); j++) {
            cout << J[i][j] << " ";
        }
        cout << endl;
    }
cout << "actual" <<endl;
for (int i = 0; i < Jtest.size(); i++) {
        for (int j = 0; j < Jtest[0].size(); j++) {
            cout << Jtest[i][j] << " ";
        }
        cout <<
        endl;
    }
*/

#endif // JACOBIAN_H_INCLUDED
