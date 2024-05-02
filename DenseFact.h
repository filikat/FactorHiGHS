#ifndef DENSE_FACT_H
#define DENSE_FACT_H

#include <time.h>
#include "Blas_declaration.h"

double GetTime() {
  struct timespec now;
  clock_gettime(CLOCK_REALTIME, &now);
  return now.tv_sec + now.tv_nsec * 1e-9;
}

#define max(i, j) ((i) >= (j) ? (i) : (j))
#define min(i, j) ((i) >= (j) ? (j) : (i))

// block size
const int nb = 256;

// variables for BLAS calls
double d_one = 1.0;
double d_zero = 0.0;
double d_m_one = -1.0;
int i_one = 1;
char LL = 'L';
char NN = 'N';
char RR = 'R';
char TT = 'T';
char UU = 'U';

#endif