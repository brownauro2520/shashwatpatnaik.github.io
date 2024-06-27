#ifndef SHAPE_H_INCLUDED
#define SHAPE_H_INCLUDED
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

/******************************************************************/
//   FUNCTION Definition: shapeL

vector<double> shapeL(vector<double> &xref, int p)
{
  /* Returns Lagrange shape functions on a triangle, for a given order
     0<=p<=3 and reference coordinates in xref.  pphi is a pointer to
     the phi vector, which gets reallocated. */
  int prank;
  double x, y;

  prank = (p+1)*(p+2)/2;
  double* phi_ptr = new double[prank];
  vector<double> phi(phi_ptr, phi_ptr + prank);
  x = xref[0];
  y = xref[1];

  switch (p) {

  case 0:
    phi[0] = 1.0;
    break;

  case 1:
    phi[0] = 1-x-y;
    phi[1] =   x  ;
    phi[2] =     y;
    break;

  case 2:
    phi[0] = 1.0-3.0*x-3.0*y+2.0*x*x+4.0*x*y+2.0*y*y;
    phi[2] = -x+2.0*x*x;
    phi[5] = -y+2.0*y*y;
    phi[4] = 4.0*x*y;
    phi[3] = 4.0*y-4.0*x*y-4.0*y*y;
    phi[1] = 4.0*x-4.0*x*x-4.0*x*y;
    break;

  case 3:
    phi[0] = 1.0-11.0/2.0*x-11.0/2.0*y+9.0*x*x+18.0*x*y+9.0*y*y-9.0/2.0*x*x*x-27.0/2.0*x*x*y-27.0/2.0*x*y*y-9.0/2.0*y*y*y;
    phi[3] = x-9.0/2.0*x*x+9.0/2.0*x*x*x;
    phi[9] = y-9.0/2.0*y*y+9.0/2.0*y*y*y;
    phi[6] = -9.0/2.0*x*y+27.0/2.0*x*x*y;
    phi[8] = -9.0/2.0*x*y+27.0/2.0*x*y*y;
    phi[7] = -9.0/2.0*y+9.0/2.0*x*y+18.0*y*y-27.0/2.0*x*y*y-27.0/2.0*y*y*y;
    phi[4] = 9.0*y-45.0/2.0*x*y-45.0/2.0*y*y+27.0/2.0*x*x*y+27.0*x*y*y+27.0/2.0*y*y*y;
    phi[1] = 9.0*x-45.0/2.0*x*x-45.0/2.0*x*y+27.0/2.0*x*x*x+27.0*x*x*y+27.0/2.0*x*y*y;
    phi[2] = -9.0/2.0*x+18.0*x*x+9.0/2.0*x*y-27.0/2.0*x*x*x-27.0/2.0*x*x*y;
    phi[5] = 27.0*x*y-27.0*x*x*y-27.0*x*y*y;
    break;

  default:
    printf("Unrecognized p in shape.\n");
    break;

  }
  return phi;
}

/******************************************************************/
//   FUNCTION Definition: gradientL

vector<double> gradientL(vector<double> &xref, int p)
{
  /* Returns gradients (x and y derivatives) of the Lagrange shape
     functions for a given order p and reference coordinates.  pghi is
     a pointer to the vector that gets reallocated.  In pgphi, all x
     derivatives are stored first, then the y derivatives. */

  int prank;
  double x, y;



  prank = (p+1)*(p+2)/2;
  double* gphi_ptr = new double[2*prank];
  vector<double> gphi(gphi_ptr, gphi_ptr + 2*prank);
  x = xref[0];
  y = xref[1];

  switch (p){

  case 0:
    gphi[0] =  0.0;
    gphi[1] =  0.0;
    break;

  case 1:
    gphi[0] =  -1.0;
    gphi[1] =  1.0;
    gphi[2] =  0.0;
    gphi[3] =  -1.0;
    gphi[4] =  0.0;
    gphi[5] =  1.0;

    break;

  case 2:
    gphi[0] =  -3.0+4.0*x+4.0*y;
    gphi[1] =  4.0-8.0*x-4.0*y;
    gphi[2] =  -1.0+4.0*x;
    gphi[3] =  -4.0*y;
    gphi[4] =  4.0*y;
    gphi[5] =  0;
    gphi[6] =  -3.0+4.0*x+4.0*y;
    gphi[7] =  -4.0*x;
    gphi[8] =  0;
    gphi[9] =  4.0-4.0*x-8.0*y;
    gphi[10] =  4.0*x;
    gphi[11] =  -1.0+4.0*y;

    break;

  case 3:
    gphi[0] =  -11.0/2.0+18.0*x+18.0*y-27.0/2.0*x*x-27.0*x*y-27.0/2.0*y*y;
    gphi[3] =  1.0-9.0*x+27.0/2.0*x*x;
    gphi[9] =  0.0;
    gphi[6] =  -9.0/2.0*y+27.0*x*y;
    gphi[8] =  -9.0/2.0*y+27.0/2.0*y*y;
    gphi[7] =  9.0/2.0*y-27.0/2.0*y*y;
    gphi[4] =  -45.0/2.0*y+27.0*x*y+27.0*y*y;
    gphi[1] =  9.0-45.0*x-45.0/2.0*y+81.0/2.0*x*x+54.0*x*y+27.0/2.0*y*y;
    gphi[2] =  -9.0/2.0+36.0*x+9.0/2.0*y-81.0/2.0*x*x-27.0*x*y;
    gphi[5] =  27.0*y-54.0*x*y-27.0*y*y;
    gphi[10] =  -11.0/2.0+18.0*x+18.0*y-27.0/2.0*x*x-27.0*x*y-27.0/2.0*y*y;
    gphi[13] =  0.0;
    gphi[19] =  1.0-9.0*y+27.0/2.0*y*y;
    gphi[16] =  -9.0/2.0*x+27.0/2.0*x*x;
    gphi[18] =  -9.0/2.0*x+27.0*x*y;
    gphi[17] =  -9.0/2.0+9.0/2.0*x+36.0*y-27.0*x*y-81.0/2.0*y*y;
    gphi[14] =  9.0-45.0/2.0*x-45.0*y+27.0/2.0*x*x+54.0*x*y+81.0/2.0*y*y;
    gphi[11] =  -45.0/2.0*x+27.0*x*x+27.0*x*y;
    gphi[12] =  9.0/2.0*x-27.0/2.0*x*x;
    gphi[15] =  27.0*x-27.0*x*x-54.0*x*y;

    break;
  default:
    cout << "you done bro" << endl;
    break;
  }
  return gphi;
} // gradientL


#endif // SHAPE_H_INCLUDED
