#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>
using namespace std;

vector<pair<double, double>> getLagrangeNodes(int p) {
    vector<pair<double, double>> nodes;
    if (p == 0) {
        nodes.push_back(make_pair(1.0/3.0, 1.0/3.0));
    } else if (p == 1) {
        nodes.push_back(make_pair(0.0, 0.0));
        nodes.push_back(make_pair(1.0, 0.0));
        nodes.push_back(make_pair(0.0, 1.0));
    } else if (p == 2){
		nodes.push_back(make_pair(0.0, 0.0));
        nodes.push_back(make_pair(0.5, 0.0));
        nodes.push_back(make_pair(1.0, 0.0));
 		nodes.push_back(make_pair(0.0, 0.5));
        nodes.push_back(make_pair(0.5, 0.5));
        nodes.push_back(make_pair(0.0, 1.0));       
    } else if (p == 3){
		nodes.push_back(make_pair(0.0, 0.0));
        nodes.push_back(make_pair(1.0/3.0, 0.0));
        nodes.push_back(make_pair(2.0/3.0, 0.0));
 		nodes.push_back(make_pair(1.0, 0.0));
		nodes.push_back(make_pair(0.0, 1.0/3.0));
        nodes.push_back(make_pair(1.0/3.0, 1.0/3.0));
        nodes.push_back(make_pair(2.0/3.0, 1.0/3.0));
        nodes.push_back(make_pair(0.0, 2.0/3.0));
        nodes.push_back(make_pair(1.0/3.0, 2.0/3.0));
        nodes.push_back(make_pair(0.0, 1.0));
    }
    return nodes;
}

vector<pair<double, double>> get1DNodesW(int order) {
    vector<pair<double, double>> nodes;
    if (order == 1) {
        nodes.push_back(make_pair(0.500000000000000, 1.000000000000000));
    } else if (order == 3) {
        nodes.push_back(make_pair(0.211324865405187, 0.500000000000000));
        nodes.push_back(make_pair(0.788675134594813, 0.500000000000000));
    } else if (order == 5){
		nodes.push_back(make_pair(0.112701665379258, 0.277777777777778));
        nodes.push_back(make_pair(0.500000000000000, 0.444444444444444));
        nodes.push_back(make_pair(0.887298334620742, 0.277777777777778));  
    } else if (order == 7){
		nodes.push_back(make_pair(0.069431844202974, 0.173927422568727));
        nodes.push_back(make_pair(0.330009478207572, 0.326072577431273));
        nodes.push_back(make_pair(0.669990521792428, 0.326072577431273));
 		nodes.push_back(make_pair(0.930568155797026, 0.173927422568727));
    } else if (order == 9){
		nodes.push_back(make_pair(0.046910077030668, 0.118463442528095));
        nodes.push_back(make_pair(0.230765344947158, 0.239314335249683));
        nodes.push_back(make_pair(0.500000000000000, 0.284444444444444));
 		nodes.push_back(make_pair(0.769234655052841, 0.239314335249683));
 		nodes.push_back(make_pair(0.953089922969332, 0.118463442528095));
    }
    return nodes;
}

int main() {		  
    int order_int_1D;
    int order_int_2D;
    int p = 3;
    int q = 2;
    // Get the number of Lagrange point and its coordinate in reference space
    int Np = (1+p)*(2+p)/2;
    vector<pair<double, double>> nodes = getLagrangeNodes(p);
    // Get the Integration order // if 1D order is an even number we do +1
    // 1D
    order_int_1D = 2*p + q;
    if (order_int_1D%2 == 0){
    order_int_1D = order_int_1D + 1;
	}
    // 2D
	order_int_2D = 2*p+1 + 2*(q-1);
    // Get the sq and weight for 1D quadrature
    vector<pair<double, double>> q_1Dnodes = get1DNodesW(order_int_1D);
    // Number of 1D quadrature point
    int n1 = q_1Dnodes.size();
  return 0;
}