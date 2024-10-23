#ifndef DENSE_FACT_H
#define DENSE_FACT_H

#include "DataCollector.h"

// dense factorization kernel
int denseFactK(char uplo, int n, double* A, int lda, const int* pivot_sign,
               double thresh, double* regul, DataCollector& DC);

// dense partial factorization, with blocks
int denseFactF(int n, int k, int nb, double* A, int lda, double* B, int ldb,
               const int* pivot_sign, double thresh, double* regul,
               DataCollector& DC);

// dense partial factorization, in blocked-hybrid format
int denseFactHP(int n, int k, int nb, double* A, double* B,
                const int* pivot_sign, double thresh, double* regul,
                DataCollector& DC);

// dense partial factorization, in blocked-hybrid format with hybrid Schur
// complement
int denseFactHH(int n, int k, int nb, double* A, double* B,
                const int* pivot_sign, double thresh, double* regul,
                DataCollector& DC);

// function to convert A from lower packed, to lower-blocked-hybrid format
int denseFactL2H(double* A, int nrow, int ncol, int nb, DataCollector& DC);

#endif