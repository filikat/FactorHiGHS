#ifndef BLAS_DECLARATION_H
#define BLAS_DECLARATION_H

// declaration of BLAS functions

#ifdef __cplusplus
extern "C" {
#endif

// level 1
void daxpy_(const int* n, const double* alpha, const double* dx,
            const int* incx, double* dy, const int* incy);
void dcopy_(const int* n, const double* dx, const int* incx, double* dy,
            const int* incy);
double ddot_(const int* n, const double* dx, const int* incx, const double* dy,
             const int* incy);
void dscal_(const int* n, const double* da, double* dx, const int* incx);
void dswap_(const int* n, double* dx, const int* incx, double* dy,
            const int* incy);

// level 2
void dgemv_(const char* trans, const int* m, const int* n, const double* alpha,
            const double* A, const int* lda, const double* x, const int* incx,
            const double* beta, double* y, const int* incy);
void dtpsv_(const char* uplo, const char* trans, const char* diag, const int* n,
            const double* ap, double* x, const int* incx);
void dtrsv_(const char* uplo, const char* trans, const char* diag, const int* n,
            const double* A, const int* lda, double* x, const int* incx);
void dger_(const int* m, const int* n, const double* alpha, const double* x,
           const int* incx, const double* y, const int* incy, double* A,
           const int* lda);

// level 3
void dgemm_(const char* transa, const char* transb, const int* m, const int* n,
            const int* k, const double* alpha, const double* A, const int* lda,
            const double* B, const int* ldb, const double* beta, double* C,
            const int* ldc);
void dsyrk_(const char* uplo, const char* trans, const int* n, const int* k,
            const double* alpha, const double* a, const int* lda,
            const double* beta, double* c, const int* ldc);
void dtrsm_(const char* side, const char* uplo, const char* trans,
            const char* diag, const int* m, const int* n, const double* alpha,
            const double* a, const int* lda, double* b, const int* ldb);

#ifdef __cplusplus
}
#endif

#endif