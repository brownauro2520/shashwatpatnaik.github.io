#include <iostream>
#include <vector>
#include <cmath>
using namespace std;


std::vector<double> Scalar(const std::vector<double>& vector, double scalar) {
    std::vector<double> result(vector.size());

    for (int i = 0; i < vector.size(); i++) {
        result[i] = vector[i] * scalar;
    }

    return result;
}


vector<double> addVectors(const vector<double>& vector1, const vector<double>& vector2) {
    vector<double> result(vector1.size());
    for (int i = 0; i < vector1.size(); i++) {
        result[i] = vector1[i] + vector2[i];
    }
    return result;
}

vector<double> subtractVectors(const vector<double>& vector1, const vector<double>& vector2) {
    vector<double> result(vector1.size());
    for (int i = 0; i < vector1.size(); i++) {
        result[i] = vector1[i] - vector2[i];
    }
    return result;
}

pair<vector<double>, double> outflow(vector<double> &U, vector<double> &n) {

//initalizing outflow condition
double gamma = 1.4;
double gmi = gamma-1;
double p_t = 1/ gamma;
double p_b = 0.9*p_t;

//computing v_+. p_+ and c_+
double u = U[1]/U[0];
double v = U[2]/U[0];
vector<double> v_int = {u, v};
double v_int2 = pow(v_int[0],2)  + pow(v_int[1],2);
double p_plus = gmi*(U[3] - (0.5*U[0]*v_int2));
double c_plus = sqrt(gamma*p_plus/U[0]);

//calculating S_+, J_+ and un+
double S_plus = p_plus/pow(U[0],gamma);
double u_n_plus = dot_product(v_int,n);
double j_plus = u_n_plus+((2*c_plus)/(gmi));

//computing boundary nodes
double rho_b = pow((p_b/S_plus),(1/gamma));
double c_b = sqrt(gamma*p_b/rho_b);
double un_b = j_plus - (2*c_b/gmi);
vector<double> v_b = addVectors(subtractVectors(v_int,Scalar(n,dot_product(v_int,n))),Scalar(n,un_b));
double v_b2 = pow(v_b[0],2)  + pow(v_b[1],2);
double rhoE_b = p_b/gmi + 0.5*rho_b*v_b2;

//computing flux
vector<double> F(4,0);
F[0] = (rho_b*v_b[0]*n[0]) + (rho_b*v_b[1]*n[1]);
F[1] = ((rho_b*pow(v_b[0],2) + p_b)*n[0]) + (rho_b*v_b[1]*v_b[0]*n[1]);
F[2] = ((rho_b*v_b[1]*v_b[0]*n[0])) + ((rho_b*pow(v_b[1],2) + p_b)*n[1]);
F[3] = (v_b[0]*(rhoE_b+p_b)*n[0]) + (v_b[1]*(rhoE_b+p_b)*n[1]);

//evaluating wave speed
double smag = sqrt(v_b2)+c_b;

return make_pair(F, smag);
}
