#ifndef CALL_AND_TIME_BLAS_H
#define CALL_AND_TIME_BLAS_H

#include "DataCollector.h"

void callAndTime_dsyrk(char uplo, char trans, int n, int k, double alpha,
                       const double* a, int lda, double beta, double* c,
                       int ldc, DataCollector& DC);
void callAndTime_dgemm(char transa, char transb, int m, int n, int k,
                       double alpha, const double* A, int lda, const double* B,
                       int ldb, double beta, double* C, int ldc,
                       DataCollector& DC);
void callAndTime_dtpsv(char uplo, char trans, char diag, int n,
                       const double* ap, double* x, int incx,
                       DataCollector& DC);
void callAndTime_dtrsv(char uplo, char trans, char diag, int n, const double* A,
                       int lda, double* x, int incx, DataCollector& DC);
void callAndTime_dtrsm(char side, char uplo, char trans, char diag, int m,
                       int n, double alpha, const double* a, int lda, double* b,
                       int ldb, DataCollector& DC);
void callAndTime_dgemv(char trans, int m, int n, double alpha, const double* A,
                       int lda, const double* x, int incx, double beta,
                       double* y, int incy, DataCollector& DC);
void callAndTime_dcopy(int n, const double* dx, int incx, double* dy, int incy,
                       DataCollector& DC);
void callAndTime_daxpy(int n, double da, const double* dx, int incx, double* dy,
                       int incy, DataCollector& DC);
void callAndTime_dscal(int n, const double da, double* dx, int incx,
                       DataCollector& DC);

int callAndTime_denseFactK(char uplo, int n, double* A, int lda,
                           const int* pivot_sign, double thresh, double* regul,
                           DataCollector& DC);

#endif