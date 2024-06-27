#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include "wall.h"
#include "roe.h"
#include "extrafunctions.h"
#include "residual.h"

using namespace std;

int main(int argc, char* argv[]) {

    //------------------------------------------------------------------------------------------------------------------
    //pre processing
    double M = 0.5; //subsonic = 0.25
    double CFL = 1;

    string mesh = argv[1]; //coarse = 0, medium = 1, fine = 2

    //--------------------------------------------------------------------------------------------------------------------
    //define free stream,

    double alpha = 8 * M_PI / 180;
    double gamma = 1.4;
    double gmi = gamma - 1;
    vector<double> u_inf = { 1, M * cos(alpha), M * sin(alpha), (1 / (gmi * gamma)) + (pow(M,2) / 2) };

    //--------------------------------------------------------------------------------------------------------------------
    //read the mesh files
    vector<vector<double>> E2N = readmatrix("matrices/fine/E2N.txt");
    vector<vector<double>> I2E = readmatrix("matrices/fine/I2E.txt");
    vector<vector<double>> B2E = readmatrix("matrices/fine/B2E.txt");
    vector<vector<double>> In = readmatrix("matrices/fine/In.txt");
    vector<vector<double>> Bn = readmatrix("matrices/fine/Bn.txt");
    vector<double> Il = readvector("matrices/fine/Il.txt");
    vector<double> Bl = readvector("matrices/fine/Bl.txt");
    vector<double> A = readvector("matrices/fine/A.txt");
    int N = E2N.size();


    //--------------------------------------------------------------------------------------------------------------------
    //initializing state @lagrange nodes
    vector<double> U;
    vector<double> R(N * 4);

    U = readvector("sols/inj/Uinj.txt");
    vector<double> dT(N);
    residual(R, dT, U, I2E, B2E, In, Bn, Il, Bl, A, u_inf, N, CFL);
    savefilev(R, "sols/res/RHh.txt");
    return 0;
}