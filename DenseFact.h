#ifndef DENSE_FACT_H
#define DENSE_FACT_H

#include "DataCollector.h"

// dense factorization kernel
int denseFactK(char uplo, int n, double* A, int lda, const int* pivot_sign,
               double thresh, double* regul, DataCollector& DC, int sn, int bl);

// dense partial factorization, in full format
int denseFactF(int n, int k, int nb, double* A, int lda, double* B, int ldb,
               const int* pivot_sign, double thresh, double* regul,
               DataCollector& DC, int sn);

// dense partial factorization, in packed format with full diagonal blocks
int denseFactFP(int n, int k, int nb, double* A, double* B,
                const int* pivot_sign, double thresh, double* regul,
                DataCollector& DC, int sn);

// dense partial factorization, in "hybrid formats"
int denseFactH(char format, int n, int k, int nb, double* A, double* B,
               const int* pivot_sign, double thresh, double* regul,
               DataCollector& DC, int sn);

// function to convert A from lower packed, to lower-blocked-hybrid format
int denseFactP2H(double* A, int nrow, int ncol, int nb, DataCollector& DC);

#endif