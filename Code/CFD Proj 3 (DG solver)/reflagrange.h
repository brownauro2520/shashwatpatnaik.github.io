#ifndef REFLAGRANGE_H_INCLUDED
#define REFLAGRANGE_H_INCLUDED
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

vector<double> reflagrange(int side, double s){
vector<double> cord;
switch(side){
case 1:
    cord = {s, 0};
    break;
case 2:
    cord = {1-s, s};
    break;
case 3:
    cord = {0,1-s};
    break;
}

return cord;
}


#endif // REFLAGRANGE_H_INCLUDED
