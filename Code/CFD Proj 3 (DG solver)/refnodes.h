#ifndef REFNODES_H_INCLUDED
#define REFNODES_H_INCLUDED
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>
using namespace std;

vector<vector<double>> getLagrangeNodes(int p) {
    vector<vector<double>> nodes;
    if (p == 0) {
        nodes.push_back({1.0/3.0, 1.0/3.0});
    } else if (p == 1) {
        nodes.push_back({0.0, 0.0});
        nodes.push_back({1.0, 0.0});
        nodes.push_back({0.0, 1.0});
    } else if (p == 2){
		nodes.push_back({0.0, 0.0});
        nodes.push_back({0.5, 0.0});
        nodes.push_back({1.0, 0.0});
 		nodes.push_back({0.0, 0.5});
        nodes.push_back({0.5, 0.5});
        nodes.push_back({0.0, 1.0});
    } else if (p == 3){
		nodes.push_back({0.0, 0.0});
        nodes.push_back({1.0/3.0, 0.0});
        nodes.push_back({2.0/3.0, 0.0});
 		nodes.push_back({1.0, 0.0});
		nodes.push_back({0.0, 1.0/3.0});
        nodes.push_back({1.0/3.0, 1.0/3.0});
        nodes.push_back({2.0/3.0, 1.0/3.0});
        nodes.push_back({0.0, 2.0/3.0});
        nodes.push_back({1.0/3.0, 2.0/3.0});
        nodes.push_back({0.0, 1.0});
    }
    return nodes;
}


#endif // REFNODES_H_INCLUDED
