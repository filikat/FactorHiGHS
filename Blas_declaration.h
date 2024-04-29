#ifndef BLAS_DECLARATION_H
#define BLAS_DECLARATION_H

// declaration of BLAS functions

#ifdef __cplusplus
extern "C" {
#endif

// level 1
void daxpy(int* n, double* alpha, double* dx, int* incx, double* dy, int* incy);
void dcopy(int* n, double* dx, int* incx, double* dy, int* incy);
double ddot(int* n, double* dx, int* incx, double* dy, int* incy);
void dscal(int* n, double* da, double* dx, int* incx);

// level 2
void dgemv(char* trans, int* m, int* n, double* alpha, const double* A, int* lda,
           double* x, int* incx, double* beta, double* y, int* incy);
void dtrsv(char* uplo, char* trans, char* diag, int* n, const double* A, int* lda,
           double* x, int* incx);

// level 3
void dgemm(char* transa, char* transb, int* m, int* n, int* k, double* alpha,
           double* A, int* lda, double* B, int* ldb, double* beta, double* C,
           int* ldc);
void dsyrk(char* uplo, char* trans, int* n, int* k, double* alpha, double* a,
           int* lda, double* beta, double* c, int* ldc);
void dtrsm(char* side, char* uplo, char* trans, char* diag, int* m, int* n,
           double* alpha, double* a, int* lda, double* b, int* ldb);

#ifdef __cplusplus
}
#endif

#endif