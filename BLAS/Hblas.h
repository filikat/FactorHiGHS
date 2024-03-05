#ifndef HBLAS_H
#define HBLAS_H

#include <iostream>

void H_dsyrk(char* uplo, char* trans, int* n, int* k, double* alpha, double* a,
             int* lda, double* beta, double* c, int* ldc);

void H_dtrsm(char* side, char* uplo, char* trans, char* diag, int* m, int* n,
             double* alpha, double* a, int* lda, double* b, int* ldb);

#endif