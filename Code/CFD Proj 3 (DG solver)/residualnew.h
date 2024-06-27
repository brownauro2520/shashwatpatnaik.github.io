#ifndef RESIDUAL_H_INCLUDED
#define RESIDUAL_H_INCLUDED
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>
#include <tuple>
#include "inviscid.h"
#include "roe.h"
#include "HLLE.h"
using namespace std;

pair<vector<vector<double>>, vector<double>> residual(
int p,
int E_num,
vector<vector<double>> U,
vector<vector<double>> E2N,
vector<vector<double>> E2NC,
vector<vector<double>> I2E,
vector<vector<double>> B2E,
vector<vector<double>> C,
vector<vector<double>> In,
vector<vector<double>> Bn,
vector<vector<double>> BnC,
vector<double> q2DcoorL,
vector<double> q2DcoorC,
vector<double> Il,
vector<double> Bl,
vector<double> u_inf,
vector<double> CURVE,
vector<vector<double>> PhiN,
vector<vector<double>> PhiC,
vector<vector<double>> GPhiNx,
vector<vector<double>> GPhiNy,
vector<double> wqN,
vector<vector<double>> GPhiCx,
vector<vector<double>> GPhiCy,
vector<double> wqC,
vector<double> wq1N,
vector<double> wq1C,
vector<double> edgephix,
vector<double> edgephiy,
vector<double> edgephiBC,
vector<double> A,
double CFL,
vector<double> J_edge_matrix)
{
	// std::vector<std::vector<double>> U(state.size()/4, std::vector<double>(4));
	// for (int i = 0; i < state.size()/4 ; i++) {
    // 	for (int j = 0; j < 4; j++) {
    // 	 U[i][j] = state[(i*4) + j];
	// }
	//     }
	// Find Np (Number of Lagrange polynomial)
	int Np = (p+1)*(p+2)/2;
	// Cal DOF (NUmber of Unknown, Same for all element)
	int Ndof = E_num*Np;
	// Create dT matrix
	vector<double> dT(Ndof,0.0);
	// Create zero Residual Matrix NDOF x 4 (Number of Unknown x 4 States)
	vector<vector<double>> R;
	for (int i = 0; i < Ndof; i++) {
		R.push_back({0.0,0.0,0.0,0.0});
		}
	// Get the size of the vector
	int NqN = PhiN.size();
	// Get the size of the vector
	int NqC = PhiC.size();

	// // Change basis vector to basis matrix
	// std::vector<std::vector<double>> PhiN(NqN, std::vector<double>(Np));
	// for (int i = 0; i < NqN ; i++) {
    // 	for (int j = 0; j < Np; j++) {
    // 	 PhiN[i][j] = PhiNV[(i*Np) + j];
	// }
	//     }
	// // Change basis vector to basis matrix
	// std::vector<std::vector<double>> PhiC(NqC, std::vector<double>(Np));
	// for (int i = 0; i < NqC ; i++) {
    // 	for (int j = 0; j < Np; j++) {
    // 	 PhiC[i][j] = PhiCV[(i*Np) + j];
	// }
	//     }
	// Change basis vector to basis matrix
	// std::vector<std::vector<double>> GPhiNx(NqN, std::vector<double>(Np));
	// for (int i = 0; i < NqN ; i++) {
    //     int j_0 = 0;
    // 	for (int j = 0; j < Np; j++) {
    // 	 GPhiNx[i][j] = GPhiNV[(i*Np*2) + j_0];
    //      j_0 = j_0 + 2;
	// }
	//     }
	// Change basis vector to basis matrix
	// std::vector<std::vector<double>> GPhiNy(NqN, std::vector<double>(Np));
	// for (int i = 0; i < NqN ; i++) {
    //     int j_0 = 1;
    // 	for (int j = 0; j < Np; j++) {
    // 	 GPhiNy[i][j] = GPhiNV[(i*Np*2) + j_0];
    //      j_0 = j_0 + 2;
	// }
	//     }
	// Change basis vector to basis matrix
	// std::vector<std::vector<double>> GPhiCx(NqC, std::vector<double>(Np));
    // for (int i = 0; i < NqC ; i++) {
    //     int j_0 = 0;
    // 	for (int j = 0; j < Np; j++) {
    // 	 GPhiCx[i][j] = GPhiCV[(i*Np*2) + j_0];
    //      j_0 = j_0 + 2;
	// }
	//     }
	// // Change basis vector to basis matrix
	// std::vector<std::vector<double>> GPhiCy(NqC, std::vector<double>(Np));
	// for (int i = 0; i < NqC ; i++) {
    //     int j_0 = 1;
    // 	for (int j = 0; j < Np; j++) {
    // 	 GPhiCy[i][j] = GPhiCV[(i*Np*2) + j_0];
    //      j_0 = j_0 + 2;
	// }
	//     }

    std::vector<std::vector<double>> q2DcordL(NqN, std::vector<double>(2));
	for (int i = 0; i < NqN ; i++) {
    	for (int j = 0; j < 2; j++) {
    	 q2DcordL[i][j] = q2DcoorL[i*2 + j];
	}
	    }

    std::vector<std::vector<double>> q2DcordC(NqC, std::vector<double>(2));
	for (int i = 0; i < NqC ; i++) {
    	for (int j = 0; j < 2; j++) {
    	 q2DcordC[i][j] = q2DcoorC[i*2 + j];
	}
	    }

	// %%%%%%%%%%%%%%%%%%%% Contributions from element interiors %%%%%%%%%%%%%%%%%%%%
    int new1 = 0;
	// Loop through each element
	for (int i_E = 0; i_E < E_num; i_E++) {

    int curve_check = CURVE[i_E];
	int Nq;
    int q;
	std::vector<std::vector<double>> Phi;
	std::vector<std::vector<double>> GPhix;
	std::vector<std::vector<double>> GPhiy;
	std::vector<double> wq;
    std::vector<std::vector<double>> q2Dcoor;
    std::vector<vector<double>> J;
    std::vector<vector<double>> INV_J(2, vector<double>(2));
    vector<double> cords;
    vector<vector<double>> cord;
    double det_J;
	if (curve_check == 1){
	Nq = NqC;
	Phi = PhiC;
	GPhix = GPhiCx;
	GPhiy = GPhiCy;
	wq = wqC;
    q = 2;
    // Coordinate of the reference point should be here
    q2Dcoor = q2DcordC;
	}
    else {
    Nq = NqN;
	Phi = PhiN;
	GPhix = GPhiNx;
	GPhiy = GPhiNy;
	wq = wqN;
    q = 1;
    // Coordinate of the reference point should be here
    q2Dcoor = q2DcordL;
	}

	// Initialize Ur in each quadrature point
	vector<vector<double>> Ur;
	for (int i = 0; i < Nq; i++) {
		Ur.push_back({0.0,0.0,0.0,0.0});
	}
	// Calculate Ur in each quadrature point
	for (int i = 0; i < Nq ; i++) {
    	for (int j = 0; j < Np; j++) {
    	 Ur[i][0] = Ur[i][0] +  Phi[i][j]*U[j + i_E*Np][0];
    	 Ur[i][1] = Ur[i][1] +  Phi[i][j]*U[j + i_E*Np][1];
    	 Ur[i][2] = Ur[i][2] +  Phi[i][j]*U[j + i_E*Np][2];
    	 Ur[i][3] = Ur[i][3] +  Phi[i][j]*U[j + i_E*Np][3];
	}
	    }

	// Calculate Flux vector @ quadrature point
	vector<vector<double>> Fx;
	vector<vector<double>> Fy;
	for (int i = 0; i < Nq ; i++) {
	double E = Ur[i][3]/Ur[i][0];
	double u = Ur[i][1]/Ur[i][0];
	double v = Ur[i][2]/Ur[i][0];
	double P = (1.4-1.0)*(Ur[i][3] - 0.5*Ur[i][0]*(pow(u,2) + pow(v,2)));
	double H = E + P/Ur[i][0];
	Fx.push_back({Ur[i][1], Ur[i][0]*pow(u,2) + P, Ur[i][0]*u*v, Ur[i][0]*u*H});
	Fy.push_back({Ur[i][2], Ur[i][0]*u*v, Ur[i][0]*pow(v,2) + P, Ur[i][0]*v*H});
		}

	// cout << Fx[0][0] << " " << Fx[0][1] << " " << Fx[0][2] << " " << Fx[0][3] << '\n';
	// cout << Fy[0][0] << " " << Fy[0][1] << " " << Fy[0][2] << " " << Fy[0][3] << '\n';

	// Calculate Residual here
	// Needed: Jacobian, Weight2D, GPhi
	for (int i = 0; i < Nq ; i++) {

        // Use ie to find the global coordinate
        // Global Lagrange Coordinate
        //if it is a linear element
        if (curve_check==0){

            //getting global coordinates for that element
            vector<double> node = E2N[i_E];
            cord = elmcord(C, node);

        }
        //if it is a curve element
        else {
            //getting global coordinates for that element
            vector<double> node = E2NC[new1];
            cord = elmcord(C, node);
            if (i == Nq-1){new1 += 1;}
        }
        // Get reference coordinate @ each quadrature point
        double xi = q2Dcoor[i][0];
        double eta = q2Dcoor[i][1];
        cords = {xi,eta};
        // Call Jacobian function, we wll get 2x2 matrix for each elem, each q point
        J = jacobian(p,q,cord, cords);
        // Do det and Inverse
        // det_J = (J[0][0]*J[1][1])-(J[1][0]*J[0][1]);
        det_J = determinantJ(J);
        INV_J[0][0] = J[1][1]/det_J;
        INV_J[0][1] = -J[0][1]/det_J;
        INV_J[1][0] = -J[1][0]/det_J;
        INV_J[1][1] = J[0][0]/det_J;
		// cout << "INV_J " << endl;
		// for (int i=0; i<2; i++){
		// 	cout << INV_J[i][0] << " " << INV_J[i][1] << endl;
		// }
    	for (int j = 0; j < Np; j++) {
        // We need to call Jacobian to obtain Jacobian matrix and the det of the jacobian
        R[j + i_E*Np][0] = R[j + i_E*Np][0] - ((GPhix[i][j]*INV_J[0][0] + GPhiy[i][j]*INV_J[1][0])*wq[i]*Fx[i][0]*det_J +(GPhix[i][j]*INV_J[0][1] + GPhiy[i][j]*INV_J[1][1])*wq[i]*Fy[i][0]*det_J);
        R[j + i_E*Np][1] = R[j + i_E*Np][1] - ((GPhix[i][j]*INV_J[0][0] + GPhiy[i][j]*INV_J[1][0])*wq[i]*Fx[i][1]*det_J +(GPhix[i][j]*INV_J[0][1] + GPhiy[i][j]*INV_J[1][1])*wq[i]*Fy[i][1]*det_J);
        R[j + i_E*Np][2] = R[j + i_E*Np][2] - ((GPhix[i][j]*INV_J[0][0] + GPhiy[i][j]*INV_J[1][0])*wq[i]*Fx[i][2]*det_J +(GPhix[i][j]*INV_J[0][1] + GPhiy[i][j]*INV_J[1][1])*wq[i]*Fy[i][2]*det_J);
        R[j + i_E*Np][3] = R[j + i_E*Np][3] - ((GPhix[i][j]*INV_J[0][0] + GPhiy[i][j]*INV_J[1][0])*wq[i]*Fx[i][3]*det_J +(GPhix[i][j]*INV_J[0][1] + GPhiy[i][j]*INV_J[1][1])*wq[i]*Fy[i][3]*det_J);
    	//  R[j + i_E*Np][0] = R[j + i_E*Np][0] - J*J*(GPhix[i][j]*wq[i]*Fx[i][0] + GPhiy[i][j]*wq[i]*Fy[i][0]);
    	//  R[j + i_E*Np][1] = R[j + i_E*Np][1] - J*J*(GPhix[i][j]*wq[i]*Fx[i][1] + GPhiy[i][j]*wq[i]*Fy[i][1]);
    	//  R[j + i_E*Np][2] = R[j + i_E*Np][2] - J*J*(GPhix[i][j]*wq[i]*Fx[i][2] + GPhiy[i][j]*wq[i]*Fy[i][2]);
    	//  R[j + i_E*Np][3] = R[j + i_E*Np][3] - J*J*(GPhix[i][j]*wq[i]*Fx[i][3] + GPhiy[i][j]*wq[i]*Fy[i][3]);
	}
	
	    }
    }
// 	for (int i =0; i < U.size(); i++){
// for (int j =0; j < 4; j++){
//     cout << R[i][j] <<'\n';
// }
// cout << endl;
// }
// cout << "Interior" << endl;
    // cout << "edge" << endl;
	int Nq1 = wq1N.size();
	// cout << "going in" <<endl;
	vector<vector<double>> edgephi;
	for (int i = 0; i < 3*Nq1*Np; i++) {
		// cout << edgephix[i] << " " << edgephiy[i] << endl;
		edgephi.push_back({edgephix[i],edgephiy[i]});
		}
	int elemL;
	int elemR;
	int edgeL;
	int edgeR;
	int IL;
	int IR;
	int PhiL_index;
	int PhiR_index;
	double PhiL;
	double PhiR;
	vector<double> smag(Nq1,0.0);
	// Contributions from interior edges
	// Looping each interior edge using I2E
	for (int ie = 0; ie < I2E.size(); ie++) {
	// Get elementL, elementR from I2E
	elemL = I2E[ie][0];
	elemR = I2E[ie][2];
	// Get edgeL, edgeR from I2E
	edgeL = I2E[ie][1];
	edgeR = I2E[ie][3];
	std::vector<std::vector<double>> matrixPhiL(Nq1, std::vector<double>(Np));
	std::vector<std::vector<double>> matrixPhiR(Nq1, std::vector<double>(Np));
	vector<vector<double>> uL;
	for (int i = 0; i < Nq1; i++) {
		uL.push_back({0.0,0.0,0.0,0.0});
		}
	vector<vector<double>> uR;
	for (int i = 0; i < Nq1; i++) {
		uR.push_back({0.0,0.0,0.0,0.0});
		}
		// cout << "edge h1" << endl;
	// Loop over each q point
	for (int q1 = 0; q1 < Nq1; q1++) {
	// Looping over each Unknown
	for (int k = 0; k < Np; k++) {

            // cout << "edge h2" << endl;
	// Get PhiL index
	// PhiL_index = (edgeL-1)*Np*Nq1 + q1*Np + k;
	// PhiL_index = (edgeL-1)*Np + q1*Np*3 + k;
	PhiL_index = (edgeL-1)*Np*Nq1 + q1*Np + k;
	PhiR_index = (edgeR-1)*Np*Nq1 + q1*Np + k;
	// cout << "PhiL_index: " << PhiL_index<< endl;
	// Get PhiR index
	// PhiR_index = (edgeR-1)*Np + q1*Np*3 + k;
	// cout << "PhiR_index: " << PhiR_index<< endl;
	// Get Phi L and R for each edge
	PhiL = edgephi[PhiL_index][0];
	PhiR = edgephi[PhiR_index][1];
	// cout << edgeL << " " << edgeR << endl;
	// cout << "PhiL: " << PhiL<< " ";
	// cout << "PhiR: " << PhiR<< endl;
	// Store as a matrix
	matrixPhiL[q1][k] = PhiL;
	matrixPhiR[q1][k] = PhiR;
    // cout << "edge1" << endl;
	// Calculate uL
	uL[q1][0] = uL[q1][0] +  matrixPhiL[q1][k]*U[(elemL - 1)*Np + k][0];
	//cout << matrixPhiL[q1][k] << " " << U[(elemL - 1)*Np + k][0] << " " << matrixPhiL[q1][k]*U[(elemL - 1)*Np + k][0] << endl;
	uL[q1][1] = uL[q1][1] +  matrixPhiL[q1][k]*U[(elemL - 1)*Np + k][1];
    uL[q1][2] = uL[q1][2] +  matrixPhiL[q1][k]*U[(elemL - 1)*Np + k][2];
    uL[q1][3] = uL[q1][3] +  matrixPhiL[q1][k]*U[(elemL - 1)*Np + k][3];

   // Calculate uR
	uR[q1][0] = uR[q1][0] +  matrixPhiR[q1][k]*U[(elemR - 1)*Np + k][0];
	uR[q1][1] = uR[q1][1] +  matrixPhiR[q1][k]*U[(elemR - 1)*Np + k][1];
    uR[q1][2] = uR[q1][2] +  matrixPhiR[q1][k]*U[(elemR - 1)*Np + k][2];
    uR[q1][3] = uR[q1][3] +  matrixPhiR[q1][k]*U[(elemR - 1)*Np + k][3];
   }
	}
	//Call the flux function to get Fhat (Flux at edge quad point)
	vector<double> n = In[ie];
    double J = Il[ie];
	for (int q1_1 = 0; q1_1 < Nq1; q1_1++) {
	vector<double> UL = uL[q1_1];
	// for (int i =0; i < UL.size(); i++){
	// 		cout << "Ul " << UL[i]<< " ";}
	// 	cout << endl;
	
    vector<double> UR = uR[q1_1];
		// for (int i =0; i < UR.size(); i++){
		// 	cout <<"Ur " << UR[i]<< " ";}
		// cout << endl;
	// Flux at each quad point
	// pair<vector<double>, double> FSI = HLLE(UL, UR, n, 1.4);
	// pair<vector<double>, double> FSI = HLLE(UL, UR, n, 1.4);
	pair<vector<double>, double> FSI = roe(UL, UR, n, 1.4);
	// cout << ie <<endl;
	// for (int i =0; i < FSI.first.size(); i++){cout << FSI.first[i]<< " ";}
	// cout <<endl;
	// Store smag as a vector
	smag[q1_1] = abs(FSI.second);
	// cout << "edge2" << endl;
	for (int k_1 = 0; k_1 < Np; k_1++) {
	// Cal the residual
	R[(elemL - 1)*Np + k_1][0] = R[(elemL - 1)*Np + k_1][0] + J*matrixPhiL[q1_1][k_1]*wq1N[q1_1]*FSI.first[0];
	// cout << J<< " " << matrixPhiL[q1_1][k_1] << " " << wq1N[q1_1] << " " << FSI.first[0]<< endl;
	R[(elemL - 1)*Np + k_1][1] = R[(elemL - 1)*Np + k_1][1] + J*matrixPhiL[q1_1][k_1]*wq1N[q1_1]*FSI.first[1];
	R[(elemL - 1)*Np + k_1][2] = R[(elemL - 1)*Np + k_1][2] + J*matrixPhiL[q1_1][k_1]*wq1N[q1_1]*FSI.first[2];
	R[(elemL - 1)*Np + k_1][3] = R[(elemL - 1)*Np + k_1][3] + J*matrixPhiL[q1_1][k_1]*wq1N[q1_1]*FSI.first[3];

	R[(elemR - 1)*Np + k_1][0] = R[(elemR - 1)*Np + k_1][0] - J*matrixPhiR[q1_1][k_1]*wq1N[q1_1]*FSI.first[0];
	R[(elemR - 1)*Np + k_1][1] = R[(elemR - 1)*Np + k_1][1] - J*matrixPhiR[q1_1][k_1]*wq1N[q1_1]*FSI.first[1];
	R[(elemR - 1)*Np + k_1][2] = R[(elemR - 1)*Np + k_1][2] - J*matrixPhiR[q1_1][k_1]*wq1N[q1_1]*FSI.first[2];
	R[(elemR - 1)*Np + k_1][3] = R[(elemR - 1)*Np + k_1][3] - J*matrixPhiR[q1_1][k_1]*wq1N[q1_1]*FSI.first[3];
	}
	}
	// cout << "Boundary" << endl;
// 	for (int i =0; i < matrixPhiR.size(); i++){
// for (int j =0; j < matrixPhiR[0].size(); j++){
//     cout << matrixPhiR[i][j] <<' ';
// }
// cout << endl;
// }
	// cout << "edge3" << endl;
	// Find maximum smag @ each edge
   auto max_smag = std::max_element(smag.begin(), smag.end());
   double max = *max_smag;
   double dl = Il[ie];
   double localt = max*dl;
   // Looping over each Unknown
	for (int k = 0; k < Np; k++) {
   dT[(elemL - 1)*Np + k] += localt;
   dT[(elemR - 1)*Np + k] += localt;
	}
}

// cout << "Boundary" << endl;
// 	for (int i =0; i < U.size(); i++){
// for (int j =0; j < 4; j++){
//     cout << R[i][j] <<'\n';
// }
// cout << endl;
// }
// cout << "edge4" << endl;
// Setup matrix for boundary
	// Get new Nq1
	int Nq1L = wq1N.size();
	int Nq1C = wq1C.size();
	int curve = -1;
	double Phi_B;
	std::vector<double> edgephiB;
    std::vector<double> wq1;
    double J_edge;
	// Contributions from Boundary edges
	// Looping each interior edge using B2E
	for (int be = 0; be < B2E.size(); be++) {
	// Check if it's curve
	int check = B2E[be][2];
	if (check == -1){ // Farfield
	Nq1 = Nq1L;
	edgephiB = edgephix;
	wq1 = wq1N;
    J_edge = Bl[be];
// 			for (int j =0; j < edgephiB.size(); j++){
//     cout << edgephiB[j] <<'\n';
// }
	}
   else {
   curve = curve +1;
	Nq1 = Nq1C;
	edgephiB = edgephiBC;
	wq1 = wq1C;
	}
// 				for (int j =0; j < edgephiB.size(); j++){
//     cout << edgephiB[j] <<'\n';
// }
	// cout << "edge5" << endl;
	vector<double> smag(Nq1,0.0);
	std::vector<std::vector<double>> matrixPhiB(Nq1, std::vector<double>(Np));
	vector<vector<double>> uB;
	for (int i = 0; i < Nq1; i++) {
		uB.push_back({0.0,0.0,0.0,0.0});
		}
	// Get element from B2E
	int elemB = B2E[be][0];
	// Get edge from B2E
	int edgeB = B2E[be][1];
	// Looping over each Unknown
	for (int k = 0; k < Np; k++) {
	// Loop over each q point
	for (int q1 = 0; q1 < Nq1; q1++) {
	// Get Phi index
	// int Phi_indexB = (edgeB-1)*Np*Nq1 + q1*Np + k;
	int Phi_indexB = (edgeB-1)*Np*Nq1 + q1*Np + k;
	// int Phi_indexB = (edgeB-1)*Np + q1*Np*3 + k;
	// cout << Phi_indexB << '\n';
	// cout << edgephiB[0] << '\n';
	Phi_B = edgephiB[Phi_indexB];
	// Phi_B =  edgephiB[0];
	// cout << Phi_B << '\n';
	// Store as a matrix
	matrixPhiB[q1][k] = Phi_B;
	// cout << Phi_B << endl;
	// Calculate u
	uB[q1][0] = uB[q1][0] +  matrixPhiB[q1][k]*U[(elemB - 1)*Np + k][0];
	// cout << matrixPhiB[q1][k] << " " << U[(elemB - 1)*Np + k][0] << " " << matrixPhiB[q1][k]*U[(elemB - 1)*Np + k][0] << endl;
	uB[q1][1] = uB[q1][1] +  matrixPhiB[q1][k]*U[(elemB - 1)*Np + k][1];
    uB[q1][2] = uB[q1][2] +  matrixPhiB[q1][k]*U[(elemB - 1)*Np + k][2];
    uB[q1][3] = uB[q1][3] +  matrixPhiB[q1][k]*U[(elemB - 1)*Np + k][3];
   }
	}
	// cout << "edge6" << endl;
	// Call the flux function to get Fhat (Flux at edge quad point)
	for (int q1_1 = 0; q1_1 < Nq1; q1_1++) {
	pair<vector<double>, double> FSI;
	vector<double> UB = uB[q1_1];
		// 	for (int i =0; i < UB.size(); i++){
		// 	cout << UB[i]<< " ";}
		// cout << endl;
	if (check == -1){ // Means farfield Linear
	// The normal at each quad point
	vector<double> n = Bn[be];
	// Flux at each quad point
	// FSI = HLLE(UB, u_inf, n, 1.4);
	FSI = roe(UB, u_inf, n, 1.4);
	// cout << "here"<< endl;
	// cout << FSI.first[0] << " " << FSI.first[1] << " " << FSI.first[2] << " " << FSI.first[3]<< endl;
	}
	else { // Curve
	// Get normal vector at each 1D quadrature point I'm gonna call it BnC
	vector<double> n = BnC[curve*Nq1 + q1_1];
	// Flux at each quad point
	//FSI = HLLE(UB, u_inf, n, 1.4);
	FSI = inviscid(UB, n);
    // Set J_edge
    J_edge = J_edge_matrix[curve*Nq1 + q1_1];
	}
	// cout << "edge8" << endl;
	// Store smag as a vector
	smag[q1_1] = abs(FSI.second);
	for (int k_1 = 0; k_1 < Np; k_1++) {
	// Cal the residual
	R[(elemB - 1)*Np + k_1][0] = R[(elemB - 1)*Np + k_1][0] + J_edge*matrixPhiB[q1_1][k_1]*wq1[q1_1]*FSI.first[0];
	// cout << J_edge*matrixPhiB[q1_1][k_1]*wq1[q1_1]*FSI.first[0] << '\n';
	// cout << matrixPhiB[q1_1][k_1] << " "<< FSI.first[0] << '\n';
	R[(elemB - 1)*Np + k_1][1] = R[(elemB - 1)*Np + k_1][1] + J_edge*matrixPhiB[q1_1][k_1]*wq1[q1_1]*FSI.first[1];
	R[(elemB - 1)*Np + k_1][2] = R[(elemB - 1)*Np + k_1][2] + J_edge*matrixPhiB[q1_1][k_1]*wq1[q1_1]*FSI.first[2];
	R[(elemB - 1)*Np + k_1][3] = R[(elemB - 1)*Np + k_1][3] + J_edge*matrixPhiB[q1_1][k_1]*wq1[q1_1]*FSI.first[3];
	// cout << R[(elemB - 1)*Np + k_1][0] << " " <<R[(elemB - 1)*Np + k_1][1] << " "<< R[(elemB - 1)*Np + k_1][2] << " "<< R[(elemB - 1)*Np + k_1][3] << endl;
	}
	// cout << "edge9 " << elemB<< endl;
	}
	// Find maximum smag @ each edge
   auto max_smag = std::max_element(smag.begin(), smag.end());
   double max = *max_smag;
   double dl = Bl[be];
   double localt = max*dl;
   // Looping over each Unknown
	for (int k = 0; k < Np; k++) {
   dT[(elemB - 1)*Np + k] += localt;
	}
	// cout << "edge10 " <<  endl;
	}
// 		for (int i =0; i < U.size(); i++){
// for (int j =0; j < 4; j++){
//     cout << R[i][j] <<' ';
// }
// cout << endl;
// }
// 
	// cout << "edge11 " <<  endl;
  //computing time stepping
int index;
for (int i=0; i < Ndof; i++){
    index = floor(i/Np);
    dT[i] = (2*CFL*A[index])/dT[i];
}
// for (int j =0; j < dT.size(); j++){
//     cout << dT[j] <<'\n';
// }
// cout << "edge12 " <<  endl;
  return make_pair(R, dT);
}
#endif // RESIDUAL_H_INCLUDED
