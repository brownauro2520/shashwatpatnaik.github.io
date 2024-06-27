#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>
#include "reflagrange.h"
#include "quad1d.h"
#include "quad2d.h"
#include "shape.h"
#include "matinv.h"
#include "refnodes.h"
#include "vectoralgebra.h"
#include "jacobian.h"
#include "basis1d.h"
#include "massmatrix.h"
#include "normalC.h"
#include "readsavefiles.h"
#include "returncord.h"
#include "normN.h"
#include "timestep.h"
#include "residual.h"


using namespace std;


int main(int argc, char* argv[]){



float M = 0.25; //subsonic = 0.25
string mesh = "0"; //coarse = 0, medium = 1, fine = 2
string flow = "S"; //always subsonic
double CFL = stod(argv[2]);
//p 0 - 0.45 CFL


string sol = string("U0") + string("P") + argv[1] + string(".txt");
string res = string("Rr0") + string("P") + argv[1] + string(".txt");


vector<double> U;
U = readvector(sol); //number should be iter -1
cout << "done U reading " << U.size()<<endl;

//--------------------------------------------------------------------------------------------------------------------
//define free stream, p, q
double alpha = 1* M_PI / 180;
double gamma = 1.4;
double gmi = gamma -1;
vector<double> u_inf = {1, M*cos(alpha), M*sin(alpha), (1/(gmi*gamma))+(pow(M,2)/2) };
const int p = stod(argv[1]);
const int q = 2; //order for curve elements
const int ql = 1; //order for linear elements
const int point = (p+1)*(p+2)/2;

//defining order of intergation
const int linear = (2*p + 1);
const int curve_i = ((2*p + 1) + (2*(q-1)));
const int curve_e = ((2*p + 1) + (q-1));








//--------------------------------------------------------------------------------------------------------------------
//read the mesh files

vector<vector<double>> E2N = readmatrix("E2N.txt");
vector<vector<double>> E2NC = readmatrix("E2NC.txt");
vector<vector<double>> C = readmatrix("C.txt");
vector<double> car = readvector("cornot.txt");
vector<vector<double>> I2E = readmatrix("I2E.txt");
vector<vector<double>> B2E = readmatrix("B2E.txt");
vector<vector<double>> In = readmatrix("In.txt");
vector<vector<double>> Bn = readmatrix("Bn.txt");
//vector<vector<double>> C = readmatrix("C"+mesh+".txt");
vector<double> Il = readvector("Il.txt");
vector<double> Bl = readvector("Bl.txt");
vector<double> A = readvector("A.txt");
int N = E2N.size();

//----------------------------------------------------------------------------------
//compute basis function @quadrature point
//obtain reference point nodes

vector<vector<double>> xrefq1 = getLagrangeNodes(1); // linear elements
vector<vector<double>> xrefq2 = getLagrangeNodes(2); //curve elements
vector<vector<double>> xlagrange = getLagrangeNodes(p); //depending on p


//q = 1 element - interior
// n, coor, wq1
tuple<double, vector<double>, vector<double>> quad_li = quad2d(linear);
vector<double> basis_li = basis2d(p, get<1>(quad_li));
vector<double> grad_li = grad2d(p, get<1>(quad_li));

//q = 1 element - edge
tuple<double, vector<double>, vector<double>> quad_le = quad1d(linear);
vector<double> basis_le = basis1d(p, get<1>(quad_le));
vector<double> basis_lecw = basis1dcw(p, get<1>(quad_le));

//curved element q >1 - interior
tuple<double, vector<double>, vector<double>> quad_ci = quad2d(curve_i);
vector<double> basis_ci = basis2d(p, get<1>(quad_ci));
vector<double> grad_ci = grad2d(p, get<1>(quad_ci));

//curved element q >1 - edge
tuple<double, vector<double>, vector<double>> quad_ce = quad1d(curve_e);
vector<double> basis_ce = basis1d(p, get<1>(quad_ce));
vector<double> basis_cecw = basis1dcw(p, get<1>(quad_ce));

//computing normal for curve

vector<vector<double>> BnC;
int pointer = 0;
for (int i = 0; i < B2E.size(); i++) {
    if (B2E[i][2] == -3){
    double s = B2E[i][1];
    double elm = B2E[i][0];
    vector<double> node = E2NC[pointer];
    vector<vector<double>> cord = elmcord(C, node);
    vector<vector<double>> tempnorm = normalC(p, q, s, get<1>(quad_ce), cord);
    BnC = combinematrix(tempnorm, BnC);
    pointer +=1;
    }
}

vector<double> J_edge = normN(BnC);
for (int i = 0; i < J_edge.size(); i++) {
    BnC[i][0] = BnC[i][0]/J_edge[i];
    BnC[i][1] = BnC[i][1]/J_edge[i];
}

//--------------------------------------------------------------------------------------------------------------------
//compute mass matrix @reference

vector<vector<double>> Mtest;
vector<vector<double>> iM;

int new1 = 0;
for (int elm = 0; elm < N; elm++){
     //if it is a linear element
     if (car[elm]==0){
         //getting global coordinates for that element
         vector<double> node = E2N[elm];
         vector<vector<double>> cord = elmcord(C, node);
         Mtest = massmatrix(basis_li, cord, get<1>(quad_li), get<2>(quad_li), ql, p);
     }
     //if it is a curve element
     else {
         //getting global coordinates for that element
         vector<double> node = E2NC[new1];
         vector<vector<double>> cord = elmcord(C, node);
         Mtest = massmatrix(basis_ci, cord, get<1>(quad_ci), get<2>(quad_ci), q, p);
         new1 += 1;
     }
     iM = combinematrix(Mtest, iM);
     Mtest.clear();
 }

cout << "done mass" <<endl;
//--------------------------------------------------------------------------------------------------------------------
//initializing state @lagrange nodes
/*
vector<double> U;
vector<double> basis_L;

for (int num = 0; num < xlagrange.size(); num++){       //computing basis at ref lagrange
    vector<double> calc = basis2d(p,xlagrange[num]);
    basis_L.insert(basis_L.end(), calc.begin(), calc.end());
    calc.clear();
}

for (int i = 0; i < N; i++){
    for(int y = 0; y< point; y++){
        vector<double> ustate(4,0);
        for (int j =(y*point); j < point+point*y; j++){
            for (int k = 0; k < 4; k++){
                ustate[k] += basis_L[j]*u_inf[k];
            }
        }
        U.insert(U.end(), ustate.begin(), ustate.end());
        ustate.clear();
    }
}
// vector<double> U = readvector("U0P2HLLE.txt");
cout << "done U " << U.size()<<endl;
*/



//--------------------------------------------------------------------------------------------------------------------
//creating a empty vector to store L1 norm
vector<double> store_res;
 double total_res = 100;
// //setting tolerance
double tolerance = pow(10,-5);
int iter = 1;
//time stepping
while (total_res > tolerance) {
// Cal R0

pair<vector<vector<double>>, vector<double>> R0 = residual(p, N, U, E2N, E2NC, I2E, B2E, C, In, Bn, BnC, get<1>(quad_li),get<1>(quad_ci),
Il, Bl, u_inf, car, basis_li, basis_ci, grad_li, get<2>(quad_li), grad_ci, get<2>(quad_ci), get<2>(quad_le), get<2>(quad_ce),
basis_le, basis_lecw, basis_ce, A, CFL, J_edge);
// for (int i =0; i < U.size(); i++){
// for (int j =0; j < 4; j++){
//     cout << R0.first[0][1] <<'\n';
// }
// cout << endl;
// }
// F0
vector<double> F0 = massinvmult(iM,R0.first, p, N);

// // 1
vector<double> U1(F0.size());
for (int i =0; i < F0.size(); i++){
        //cout<<i<< " "; cout << R0.second[i]<<endl;
     U1[i] = U[i]+R0.second[floor(i/(4))]*F0[i];
}

// // Cal R1
pair<vector<vector<double>>, vector<double>> R1 = residual(p, N, U1, E2N, E2NC, I2E, B2E, C, In, Bn, BnC, get<1>(quad_li),get<1>(quad_ci),
Il, Bl, u_inf, car, basis_li, basis_ci, grad_li, get<2>(quad_li), grad_ci, get<2>(quad_ci), get<2>(quad_le), get<2>(quad_ce),
basis_le, basis_lecw, basis_ce, A, CFL, J_edge);

// // F1
vector<double> F1 = massinvmult(iM,R1.first, p, N);


// // Final
vector<double> Ufinal(U.size());
for (int i =0; i < U.size(); i++){
     Ufinal[i] = U[i]+(0.5*R0.second[floor(i/(4))])*(F0[i] + F1[i]);
 }
 U = Ufinal;
// //computing total residual
 total_res = sum_V(R0.first);
 cout << "Residual " << total_res << " iter " << iter << " " << p <<endl;
iter +=1;
// //storing residual
 store_res.push_back(abs(total_res));

if (iter%1000 == 0 ){
   savefilev(U, sol);
savefilev(store_res, res);
}
// //error is less than tolerance
if (abs(total_res) < tolerance){break;}
 }
savefilev(U, sol);
savefilev(store_res, res);
//saving state and residual to text


return 0;
}
