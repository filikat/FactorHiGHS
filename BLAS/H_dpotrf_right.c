#include "Hblas.h"

void dsyrk(char* uplo, char* trans, int* n, int* k, double* alpha, double* a,
           int* lda, double* beta, double* c, int* ldc);

void dtrsm(char* side, char* uplo, char* trans, char* diag, int* m, int* n,
           double* alpha, double* a, int* lda, double* b, int* ldb);

void dgemm(char* transa, char* transb, int* m, int* n, int* k, double* alpha,
           double* A, int* lda, double* B, int* ldb, double* beta, double* C,
           int* ldc);

void H_dpotf2(char* uplo, int* n, double* restrict A, int* ldA, int* info);

void H_dpotrf_right(char* uplo, int* n, double* restrict A, int* ldA, int* info,
              int* k) {
  // ===========================================================================
  // HiGHS BLAS function
  // DPOTRF: Double POsitive definite TRiangular Factorization
  // This version calls level 3 BLAS functions
  // Right looking version
  //
  // Perform one of:
  //  A = U^T * U   [1]
  //  A = L^T * L   [2]
  //
  // Arguments:
  // - uplo   : perform [1] if uplo = 'U','u';
  //            perform [2] if uplo = 'L','l'.
  // - n      : dimension of matrix A.
  // - A      : array of size (lda * n). To be accessed by columns.
  // - lda    : leading dimension of A.
  // - info   : information about outcome.
  // - k      : number of columns to factorize;
  //            if k < n, a partial factorization is computed.
  //
  //
  // Filippo Zanetti, 2024
  // ===========================================================================

  // ===========================================================================
  // Check input
  // ===========================================================================
  if (!uplo || !n || !A || !ldA || !info || !k) {
    printf("Invalid pointer\n");
    return;
  }

  int upper = (*uplo == 'U' || *uplo == 'u');
  int na = *n;
  int lda = *ldA;
  *info = 0;

  if (!upper && !(*uplo == 'L' || *uplo == 'l')) {
    printf("Invalid parameter uplo\n");
    *info = -1;
    return;
  }
  if (*n < 0) {
    printf("Invalid parameter n\n");
    *info = -2;
    return;
  }
  if (lda < max(1, na)) {
    printf("Invalid parameter lda\n");
    *info = -4;
    return;
  }
  if (*k < 0) {
    printf("Invalid parameter k\n");
    *info = -6;
    return;
  }

  // Quick return
  if (*n == 0) return;

  // block size
  const int nb = 64;

  // ===========================================================================
  // Blocked case
  // ===========================================================================
  char Transa, Transb, Side, Unit;
  int N, M, K;
  double MinusOne = -1.0;
  double One = 1.0;

  if (upper) {
    // Compute A = U^T * U using BLAS level 3

    Transa = 'T';
    Transb = 'N';
    Side = 'L';
    Unit = 'N';

    // j is the starting col of the block of columns
    for (int j = 0; j < *k; j += nb) {
      // jb is the size of the block
      int jb = min(nb, *k - j);

      N = jb;
      K = j;
      M = na - j - jb;

      // factorize diagonal block
      H_dpotf2(uplo, &N, &A[j + lda * j], ldA, info);
      if (*info != 0) {
        *info += j - 1;
        return;
      }

      if (j + jb < na) {
        // solve block of columns with diagonal block
        dtrsm(&Side, uplo, &Transa, &Unit, &N, &M, &One, &A[j + lda * j], ldA,
              &A[j + (j + jb) * lda], ldA);

        // update Schur complement
        dsyrk(uplo, &Transa, &M, &N, &MinusOne, &A[j + lda * (j + jb)], ldA,
              &One, &A[j + jb + lda * (j + jb)], ldA);
      }
    }
  } else {
    // Compute A = L * L^T using BLAS level 3

    Transa = 'N';
    Transb = 'T';
    Side = 'R';
    Unit = 'N';

    // j is the starting col of the block of columns
    for (int j = 0; j < *k; j += nb) {
      // jb is the size of the block
      int jb = min(nb, *k - j);

      N = jb;
      K = j;
      M = na - j - jb;

      // factorize diagonal block
      H_dpotf2(uplo, &N, &A[j + lda * j], ldA, info);
      if (*info != 0) {
        *info += j - 1;
        return;
      }

      if (j + jb < na) {
        // solve block of columns with diagonal block
        dtrsm(&Side, uplo, &Transb, &Unit, &M, &N, &One, &A[j + lda * j], ldA,
              &A[j + jb + lda * j], ldA);

        // update Schur complement
        dsyrk(uplo, &Transa, &M, &N, &MinusOne, &A[j + jb + lda * j], ldA, &One,
              &A[j + jb + lda * (j + jb)], ldA);
      }
    }
  }
}