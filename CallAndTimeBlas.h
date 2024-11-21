#ifndef CALL_AND_TIME_BLAS_H
#define CALL_AND_TIME_BLAS_H

#include "DataCollector.h"

// level 1
void callAndTime_daxpy(int n, double da, const double* dx, int incx, double* dy,
                       int incy, DataCollector& DC);
void callAndTime_dcopy(int n, const double* dx, int incx, double* dy, int incy,
                       DataCollector& DC);
void callAndTime_dscal(int n, const double da, double* dx, int incx,
                       DataCollector& DC);
void callAndTime_dswap(int n, double* dx, int incx, double* dy, int incy,
                       DataCollector& DC);

// level 2
void callAndTime_dgemv(char trans, int m, int n, double alpha, const double* A,
                       int lda, const double* x, int incx, double beta,
                       double* y, int incy, DataCollector& DC);
void callAndTime_dtpsv(char uplo, char trans, char diag, int n,
                       const double* ap, double* x, int incx,
                       DataCollector& DC);
void callAndTime_dtrsv(char uplo, char trans, char diag, int n, const double* A,
                       int lda, double* x, int incx, DataCollector& DC);
void callAndTime_dger(int m, int n, double alpha, const double* x, int incx,
                      const double* y, int incy, double* A, int lda,
                      DataCollector& DC);

// level 3
void callAndTime_dgemm(char transa, char transb, int m, int n, int k,
                       double alpha, const double* A, int lda, const double* B,
                       int ldb, double beta, double* C, int ldc,
                       DataCollector& DC);
void callAndTime_dsyrk(char uplo, char trans, int n, int k, double alpha,
                       const double* a, int lda, double beta, double* c,
                       int ldc, DataCollector& DC);
void callAndTime_dtrsm(char side, char uplo, char trans, char diag, int m,
                       int n, double alpha, const double* a, int lda, double* b,
                       int ldb, DataCollector& DC);

// kernel
int callAndTime_denseFactK(char uplo, int n, double* A, int lda,
                           int* pivot_sign, double thresh, double* regul,
                           int* swaps, double* pivot_2x2, DataCollector& DC,
                           int sn, int bl);

#endif