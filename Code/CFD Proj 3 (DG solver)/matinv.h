#ifndef MATINV_H_INCLUDED
#define MATINV_H_INCLUDED
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

#define MEPS 1e-15

/******************************************************************/
// FUNCTION Definition: computePLU
static int computePLU(vector<double> &A, vector<int> &P, int r) {
  /* Computes the permuted LU factorization of A, which is rxr.  The
     LU decomposition overwrites A. */

  int k, j, i, kmax, tk;
  int kr, jr, ik;
  double l, Akk, Amax, tmp, t0, t1, t2, t3;

  /* initialize P vector */
  for (k = 0; k < r; k++) P[k] = k;

  for (k = 0; k < (r - 1); k++) { /* loop through columns */
    kr = k * r;

    /* find maximum (absolute value) entry below pivot */
    Amax = fabs(A[kr + k]);
    kmax = k;
    for (i = (k + 1); i < r; i++) {
      if (fabs(A[i * r + k]) > Amax) {
        Amax = fabs(A[i * r + k]);
        kmax = i;
      }
    }

    /* switch rows k and kmax if necessary */
    if (kmax != k) {
      for (i = 0; i < r; i++) {
        tmp = A[kr + i];
        A[kr + i] = A[kmax * r + i];
        A[kmax * r + i] = tmp;
      }
      tk = P[k];
      P[k] = P[kmax];
      P[kmax] = tk;
    }

    Akk = A[kr + k];
    if (fabs(Akk) < MEPS) {
      cout << "Warning: Matrix near singular\n";
      return -1;
    }

    for (j = (k + 1); j < r; j++) {
      jr = j * r;
      A[jr + k] = A[jr + k] / Akk;
    }

    i = k + 1;
    while ((i + 3) < r) {
      ik = kr + i;
      t0 = A[ik + 0];
      t1 = A[ik + 1];
      t2 = A[ik + 2];
      t3 = A[ik + 3];
      for (j = (k + 1); j < r; j++) {
        jr = j * r + i;
        l = A[j * r + k];
        A[jr + 0] -= l * t0;
        A[jr + 1] -= l * t1;
        A[jr + 2] -= l * t2;
        A[jr + 3] -= l * t3;
      }
      i += 4;
    }
    while (i < r) {
      ik = kr + i;
      t0 = A[ik];
      for (j = (k + 1); j < r; j++) {
        jr = j * r;
        A[jr + i] -= A[jr + k] * t0;
      }
      i += 1;
    }
  }

  return 0;
} /* end computePLU */


static int solvePLU(vector<double>& LU, vector<int>& P, vector<double>& b, vector<double>& u, int r)
{
    /* Solves P^-1 (LU) u = b, where LU is rxr. */

    int k, kr, j;
    vector<double> y(r, 0.0);

    /* Solve Ly=Pb */
    for (k = 0; k < r; k++){
        kr = k * r;
        double temp = b[P[k]];
        for (j = 0; j < k; j++)
            temp -= LU[kr + j] * y[j];
        y[k] = temp;
    }

    /* Solve Uu=y */
    for (k = (r - 1); k >= 0; k--){
        kr = k * r;
        double temp = y[k];
        for (j = (k + 1); j < r; j++)
            temp -= LU[kr + j] * u[j];
        u[k] = temp / LU[kr + k];
    }

    return 0;
}

int invmat(vector<double> &A, vector<vector<double>> &iA, int r) {
  /* Sets iA = A^-1, where the matrices are rxr */
  int i, j, ierr;
  vector<int> P(r);
  vector<double> I(r), ia(r);

  ierr = computePLU(A, P, r);
  if (ierr != 0) {return ierr;}


  for (i = 0; i < r; i++) {
    for (j = 0; j < r; j++){
      I[j] = 0;}
    I[i] = 1.0;

    ierr = solvePLU(A, P, I, ia, r);

    if (ierr != 0) {return ierr;}

    for (j = 0; j < r; j++){
      iA[j][i] = ia[j];}
  }


  return 0;
}


#endif // MATINV_H_INCLUDED