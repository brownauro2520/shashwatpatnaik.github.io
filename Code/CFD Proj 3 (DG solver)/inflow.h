#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

double dot_product(const vector<double> &v1, const vector<double> &v2) {
    double result = 0;
    if (v1.size() != v2.size()) {
        cout << "Invalid size" << endl;
        return 0;
    }
    for (int i = 0; i < v1.size(); i++) {
        result += v1[i] * v2[i];
    }
    return result;
}

pair<vector<double>, double> inflow(vector<double> &U, vector<double> &n) {
//inflow condition
double gamma = 1.4;
double c = 1;
double rho = 1;
double alpha = 1*M_PI/180;
double gmi = gamma-1;
double gmratio = gamma/gmi;


//calculating v_+, p_+, and c_+
double u = U[1]/U[0];
double v = U[2]/U[0];
vector<double> v_int = {u, v};
double v_int2 = pow(v_int[0],2)  + pow(v_int[1],2);
double p_plus = gmi*(U[3] - (0.5*U[0]*v_int2));
double c_plus = sqrt(gamma*p_plus/U[0]);

//calculating u_n_+
double u_n_plus = dot_product(v_int,n);
//calculating j+
double j_plus = u_n_plus+((2*c_plus)/(gmi));

//initializing Tt*R,Pt and dn
double p_t = pow(c,2)*rho / gamma;
double T_t_R = p_t/(rho);

//computing N_in and dn
vector<double> normal_in = {cos(alpha), sin(alpha)};
double dn = dot_product(normal_in,n);

//calculating Mb
double Mb;
double a_Mb = (gamma * T_t_R * pow(dn, 2)) - (gmi / 2) * pow(j_plus, 2);
double b_Mb = (4 * gamma * T_t_R * dn / gmi);
double c_Mb = (4 * gamma * T_t_R) / (pow(gmi, 2)) - pow(j_plus, 2);

double discriminant = b_Mb * b_Mb - 4 * a_Mb * c_Mb;

if (discriminant >= 0) {
    double x1 = (-b_Mb + sqrt(discriminant)) / (2 * a_Mb);
    double x2 = (-b_Mb - sqrt(discriminant)) / (2 * a_Mb);
    if (x1 < 0 && x2 >= 0){Mb = x2;}
    else if (x2 < 0 && x1 >= 0){Mb = x1;}
    else if (x2 >= 0 and x1 >= 0){
            if (x2<=x1){Mb = x2;}
            else {Mb = x1;}
            }
    else {
            if (x2>=x1){Mb = x2;}
            else {Mb = x1;}
         }
} else {
    std::cout << "The equation has no real solutions for Mb." << std::endl;
}


//evaluating boundary states
double T_b_R = T_t_R / (1+(0.5*gmi*pow(Mb,2)));
double p_b = p_t*pow((T_b_R/T_t_R),(gmratio));
double rho_b = p_b/(T_b_R);
double c_b = sqrt(gamma*p_b/rho_b);
vector<double> v_b = {Mb*c_b*normal_in[0],Mb*c_b*normal_in[1]};
double v_b2 = pow(v_b[0],2)  + pow(v_b[1],2);
double rhoE_b = p_b/gmi  +  0.5*rho_b*v_b2;

//computing Flux
vector<double> F(4,0);
F[0] = (rho_b*v_b[0]*n[0]) + (rho_b*v_b[1]*n[1]);
F[1] = ((rho_b*pow(v_b[0],2) + p_b)*n[0]) + (rho_b*v_b[1]*v_b[0]*n[1]);
F[2] = ((rho_b*v_b[1]*v_b[0]*n[0])) + ((rho_b*pow(v_b[1],2) + p_b)*n[1]);
F[3] = (v_b[0]*(rhoE_b+p_b)*n[0]) + (v_b[1]*(rhoE_b+p_b)*n[1]);



//evaluating wave speed
double smag = sqrt(v_b2)+c_b;

return make_pair(F, smag);
}
