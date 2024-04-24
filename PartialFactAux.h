#ifndef PARTIAL_FACT_H
#define PARTIAL_FACT_H

#include <time.h>

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

// BLAS declaration
void dcopy(int* n, double* dx, int* incx, double* dy, int* incy);
void dscal(int* n, double* da, double* dx, int* incx);
double ddot(int* n, double* dx, int* incx, double* dy, int* incy);
void dgemv(char* trans, int* m, int* n, double* alpha, double* A, int* lda,
           double* x, int* incx, double* beta, double* y, int* incy);
void dtrsm(char* side, char* uplo, char* trans, char* diag, int* m, int* n,
           double* alpha, double* a, int* lda, double* b, int* ldb);
void dsyrk(char* uplo, char* trans, int* n, int* k, double* alpha, double* a,
           int* lda, double* beta, double* c, int* ldc);
void dgemm(char* transa, char* transb, int* m, int* n, int* k, double* alpha,
           double* A, int* lda, double* B, int* ldb, double* beta, double* C,
           int* ldc);

#endif