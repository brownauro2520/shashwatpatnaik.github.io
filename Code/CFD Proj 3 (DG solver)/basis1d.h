#ifndef BASIS1D_H_INCLUDED
#define BASIS1D_H_INCLUDED
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

vector<double> basis1d(int p, vector<double> quad){
vector<double> basis;
for (int edge = 1; edge < 4; edge++){
    for (int i = 0; i < quad.size(); i++) {
        vector<double> xref = reflagrange(edge,quad[i]);
        vector<double> calc = shapeL(xref,p);
        basis.insert(basis.end(), calc.begin(), calc.end());
        calc.clear();

    }
}
return basis;
}

vector<double> basis1dcw(int p, vector<double> quad){
vector<double> basis;
for (int edge = 1; edge < 4; edge++){
    for (int i = 0; i < quad.size(); i++) {
        vector<double> xref = reflagrange(edge,1-quad[i]);
        vector<double> calc = shapeL(xref,p);
        basis.insert(basis.end(), calc.begin(), calc.end());
        calc.clear();

    }
}
return basis;
}

vector<double> basis2d(int p, vector<double> quad){
vector<double> calc;
vector<double> basis;
for (int i = 0; i < quad.size(); i+=2) {
    vector<double> xref = {quad[i], quad[i+1]};
    calc = shapeL(xref,p);
    basis.insert(basis.end(), calc.begin(), calc.end());
    calc.clear();
}
return basis;
}

vector<double> grad2d(int p, vector<double> quad){
vector<double> calc;
vector<double> grad;
for (int i = 0; i < quad.size(); i+=2) {
    vector<double> xref = {quad[i], quad[i+1]};
    calc = gradientL(xref,p);
    grad.insert(grad.end(), calc.begin(), calc.end());
    calc.clear();
}
return grad;
}




#endif // BASIS1D_H_INCLUDED
