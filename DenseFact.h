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
int denseFactFH(char format, int n, int k, int nb, double* A, double* B,
                const int* pivot_sign, double thresh, double* regul,
                DataCollector& DC, int sn);
int denseFactFH_2(char format, int n, int k, int nb, double* A, double* B,
                  const int* pivot_sign, double thresh, double* regul,
                  DataCollector& DC, int sn);

// function to convert A from lower packed, to lower-blocked-hybrid format
int denseFactFP2FH(double* A, int nrow, int ncol, int nb, DataCollector& DC);

// objects to call blas in parallel
// they are needed, otherwise the lambda is too large to fit into a task

class GemmCaller {
  char transa_;
  char transb_;
  int m_;
  int n_;
  int k_;
  double alpha_;
  const double* A_;
  int lda_;
  const double* B_;
  int ldb_;
  double beta_;
  double* C_;
  int ldc_;
  DataCollector& DC_;

 public:
  GemmCaller(char transa, char transb, int m, int n, int k, double alpha,
             const double* A, int lda, const double* B, int ldb, double beta,
             double* C, int ldc, DataCollector& DC);
  void run();
};

class TrsmCaller {
  char side_;
  char uplo_;
  char trans_;
  char diag_;
  int m_;
  int n_;
  double alpha_;
  const double* A_;
  int lda_;
  double* B_;
  int ldb_;
  DataCollector& DC_;

 public:
  TrsmCaller(char side, char uplo, char trans, char diag, int m, int n,
             double alpha, const double* A, int lda, double* B, int ldb,
             DataCollector& DC);
  void run();
};

class GemmCaller_2 {
  int n_;
  int k_;
  int nb_;
  double* A_;
  double* T_;
  double* R_;
  int* diag_start_;
  DataCollector& DC_;

 public:
  GemmCaller_2(int n, int k, int nb, double* A, double* T, double* R,
               int* diag_start, DataCollector& DC);
  void run(int jj, int j, int v);
};

#endif