#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include "wall.h"
#include "extrafunctions.h"
#include "residual.h"
#include <adolc/adolc.h>
#include <Eigen/Sparse>
using namespace std;
using namespace Eigen;
using Eigen::MatrixXd;

int main(int argc, char* argv[]) {
    double M = 0.5; //subsonic = 0.25
    double CFL = 1;

    double alpha = 8 * M_PI / 180;
    double gamma = 1.4;
    double gmi = gamma - 1;
    vector<adouble> u_inf = { 1, M * cos(alpha), M * sin(alpha), (1 / (gmi * gamma)) + (pow(M,2) / 2) };

    //--------------------------------------------------------------------------------------------------------------------
    //read the mesh files
    cout << "Reading Files" << endl;
    vector<vector<double>> E2N = readmatrix("matrices/fine/E2N.txt");
    vector<vector<double>> I2E = readmatrix("matrices/fine/I2E.txt");
    vector<vector<double>> B2E = readmatrix("matrices/fine/B2E.txt");
    vector<vector<double>> In = readmatrix("matrices/fine/In.txt");
    vector<vector<double>> Bn = readmatrix("matrices/fine/Bn.txt");
    vector<double> Il = readvector("matrices/fine/Il.txt");
    vector<double> Bl = readvector("matrices/fine/Bl.txt");
    vector<double> A = readvector("matrices/fine/A.txt");
    int N = E2N.size();

    vector<double> dT(N);
    vector<adouble> U(N*4);
    vector<adouble> R(N * 4);

    int n = N*4; // input dimension is n

    vector<double> UU(N*4);
    vector<double> RR(N*4);
    UU = readvector("sols/U.txt");

    trace_on(1); // start tracking computation with ADOL-C
    for (int i = 0; i < n; i++)
        U[i] <<= UU[i]; // set the values of the input variables

    residual(R, dT, U, I2E, B2E, In, Bn, Il, Bl, A, u_inf, N, CFL);
    
    for(int i = 0; i < n; i++)
        R[i] >>= RR[i]; // Use >>= to tell ADOL-C that y[i] are the output variables

    trace_off();

    double x[n]; // independent vector x
    for (int i = 0; i < n; i++)
        x[i] = UU[i]; // set the pointers to the input values
        
    //////////////////////////////////////////
    /* coordinate format for Jacobian */
    unsigned int *rind  = NULL;        /* row indices    */
    unsigned int *cind  = NULL;        /* column indices */
    double       *values = NULL;       /* values         */
    int nnz;
    int options[4];

    options[0] = 0;          /* sparsity pattern by index domains (default) */ 
    options[1] = 0;          /*                         safe mode (default) */ 
    options[2] = 0;          /*              not required if options[0] = 0 */ 
    options[3] = 0;          /*                column compression (default) */ 

    sparse_jac(1, n, n, 0, x, &nnz, &rind, &cind, &values, options);
    cout << "Jacobian Solved" << endl;

    //////////////////////////////////////////
    vector<double> b_vec;
    b_vec = readvector("sols/outjacs/DLDU.txt");
    cout << "Output Jacobian Formed" << endl;
    
    // Create the sparse matrix object from the vectors
    Eigen::SparseMatrix<double> J(n, n);
    for (int i = 0; i < nnz; i++) {
        J.coeffRef(rind[i], cind[i]) = values[i];
    }
    J.makeCompressed();

    // Define the right-hand side vector
    Eigen::VectorXd b(b_vec.size());
    b << Eigen::Map<Eigen::VectorXd>(b_vec.data(), b_vec.size());

    // Solve the system of equations using SparseLU solver
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(J);
    solver.factorize(J);
    Eigen::VectorXd adj = solver.solve(b);

    // Print the solution vector
    //std::cout << "Solution vector: " << std::endl << adj << std::endl;
    
    free(rind); rind=NULL;
    free(cind); cind=NULL;
    free(values); values=NULL;
    std::ofstream file("sols/adjoints/adj.txt");
    IOFormat fmt(StreamPrecision, DontAlignCols, ", ", "\n");
    file << adj.format(fmt);
    file.close();
    cout << "Saved Adjoint" << endl;
    return 0;
}