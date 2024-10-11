#ifndef CALL_AND_TIME_BLAS_H
#define CALL_AND_TIME_BLAS_H

void callAndTime_dsyrk(char uplo, char trans, int n, int k, double alpha,
                       const double* a, int lda, double beta, double* c,
                       int ldc, double* times);
void callAndTime_dgemm(char transa, char transb, int m, int n, int k,
                       double alpha, const double* A, int lda, const double* B,
                       int ldb, double beta, double* C, int ldc, double* times);
void callAndTime_dtrsm(char side, char uplo, char trans, char diag, int m,
                       int n, double alpha, const double* a, int lda, double* b,
                       int ldb, double* times);

void callAndTime_dcopy(int n, const double* dx, int incx, double* dy, int incy,
                       double* times);
void callAndTime_daxpy(int n, double da, const double* dx, int incx, double* dy,
                       int incy, double* times);
void callAndTime_dscal(int n, const double da, double* dx, int incx,
                       double* times);

int callAndTime_fduf(char uplo, int n, double* A, int lda, double thresh,
                     double* regul, double* times);
int callAndTime_fiuf(char uplo, int n, double* A, int lda,
                     const int* pivot_sign, double thresh, double* regul,
                     int* n_reg_piv, double* times);

#endif