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
    cout << "Preporocessing" << endl;

    double M = 0.5; //subsonic = 0.25
    double CFL = 1;

    string mesh = argv[1]; //coarse = 0, medium = 1, fine = 2
    string space = argv[3];

    //--------------------------------------------------------------------------------------------------------------------
    //define free stream,
    cout << "Obtaining free stream for mesh" << mesh  << endl;

    double alpha = 8 * M_PI / 180;
    double gamma = 1.4;
    double gmi = gamma - 1;
    vector<double> u_inf = { 1, M * cos(alpha), M * sin(alpha), (1 / (gmi * gamma)) + (pow(M,2) / 2) };

    //--------------------------------------------------------------------------------------------------------------------
    //read the mesh files
    cout << "Reading files" << endl;
    vector<vector<double>> E2N = readmatrix("matrices/"+space+"/E2N.txt");
    vector<vector<double>> I2E = readmatrix("matrices/"+space+"/I2E.txt");
    vector<vector<double>> B2E = readmatrix("matrices/"+space+"/B2E.txt");
    vector<vector<double>> In = readmatrix("matrices/"+space+"/In.txt");
    vector<vector<double>> Bn = readmatrix("matrices/"+space+"/Bn.txt");
    vector<double> Il = readvector("matrices/"+space+"/Il.txt");
    vector<double> Bl = readvector("matrices/"+space+"/Bl.txt");
    vector<double> A = readvector("matrices/"+space+"/A.txt");
    int N = E2N.size();


    //--------------------------------------------------------------------------------------------------------------------
    //initializing state @lagrange nodes

    cout << "Initiliazing States" << endl;

    vector<double> U(N*4);
    if (string(argv[2]) == "inf"){
        for (int i = 0; i < N*4; i+=4) {
            U[i] = u_inf[0]; U[i+1] = u_inf[1]; U[i+2] = u_inf[2]; U[i+3] = u_inf[3];
        }
    } else if (string(argv[2]) == "load"){
        
        if (space == "coarse"){
            U = readvector("sols/U"+mesh+".txt");
            
        } else if (space == "fine"){
            U = readvector("sols/inj/Uinj.txt");
        }
    }

    //--------------------------------------------------------------------------------------------------------------------
    //time stepping

    cout << "Time Stepping \n";

    //allocating dT vector memory


    //creating a empty vector to store L1 norm
    vector<double> store_res;

    //setting tolerance
    double tolerance = pow(10, -5);

    //allocating f0, f1, u_half
    vector<double> f0(N * 4);
    vector<double> f1(N * 4);
    vector<double> u_half(N * 4);

    // in the loop ------
    double total_res;
    vector<double> dT(N);
    vector<double> dT_half(N);
    vector<double> R(N * 4);
    vector<double> R_half(N * 4);


    for (int iter = 0; iter < 1000000; iter++) {

        //initializing Residual and dT vector to zero
        fill(dT.begin(), dT.end(), 0.0);
        fill(dT_half.begin(), dT_half.end(), 0.0);
        fill(R.begin(), R.end(), 0.0);
        fill(R_half.begin(), R_half.end(), 0.0);
        total_res = 0;

        //first step, computing f0
        residual(R, dT, U, I2E, B2E, In, Bn, Il, Bl, A, u_inf, N, CFL);

        for (int i = 0; i < N * 4; i++) {
            f0[i] = R[i] / (-A[floor(i / 4)]);
        }


        //half-step
        // computing u_half
        for (int i = 0; i < N * 4;i++) {
            u_half[i] = U[i] + (f0[i] * dT[floor(i / 4)]);
        }


        //creating f1
        residual(R_half, dT_half, u_half, I2E, B2E, In, Bn, Il, Bl, A, u_inf, N, CFL);



        for (int i = 0; i < N * 4; i++) {
            f1[i] = R_half[i] / (-A[floor(i / 4)]);
        }

        //full step, marching forward
        for (int i = 0; i < N * 4; i++) {
            U[i] += (f0[i] + f1[i]) * 0.5 * dT[floor(i / 4)];
        }

        //computing total residual
        for (int i = 0; i < N * 4;i++) {
            total_res += abs(R[i]);
        }

        cout << "Residual " << total_res << " iter " << iter << endl;

        //storing residual
        store_res.push_back(abs(total_res));

        //error is less than tolerance
        if (abs(total_res) < tolerance) {
            if (space == "fine"){
                savefilev(U, "sols/U.txt");
            } else if (space == "coarse"){
                savefilev(U, "sols/U"+mesh+".txt");
            }
            
            break;
        }

    }
    return 0;
}