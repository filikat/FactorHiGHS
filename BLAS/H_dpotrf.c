#include "Hblas.h"

void dsyrk(char* uplo, char* trans, int* n, int* k, double* alpha, double* a,
           int* lda, double* beta, double* c, int* ldc);

void dtrsm(char* side, char* uplo, char* trans, char* diag, int* m, int* n,
           double* alpha, double* a, int* lda, double* b, int* ldb);

void dgemm(char* transa, char* transb, int* m, int* n, int* k, double* alpha,
           double* A, int* lda, double* B, int* ldb, double* beta, double* C,
           int* ldc);

void H_dpotf2(char* uplo, int* n, double* restrict A, int* ldA, int* info);

void H_dpotrf(char* uplo, int* n, double* restrict A, int* ldA,
              double* restrict B, int* ldB, int* info, int* k) {
  // ===========================================================================
  // HiGHS BLAS function
  // DPOTRF: Double POsitive definite TRiangular Factorization
  // This version calls level 3 BLAS functions
  // Partial factorization version
  //
  // Perform partial factorization of matrix M with left-looking approach.
  // In lower format:
  // A is used to access the first k columns, i.e., M(0:n-1,0:k-1).
  // B is used to access the remaining lower triangle, i.e., M(k:n-1,k:n-1).
  // In upper format:
  // A is used to access the first k rows, i.e., M(0:k-1,0:n-1).
  // B is used to access the remaining upper triangle, i.e., M(k:n-1,k:n-1).
  //
  // Arguments:
  // - uplo   : Upper triangular form if uplo = 'U','u';
  //            Lower triangular form if uplo = 'L','l'.
  // - n      : Dimension of matrix M.
  // - A      : Array of size (lda * k). To be accessed by columns.
  //            On input, it contains the first k columns/rows of M.
  //            On output, it contains the trapezoidal factor of the first k
  //            columns/rows of M.
  // - lda    : Leading dimension of A.
  //            It must be at least n, for lower format, and at least k for
  //            upper format.
  // - B      : Array of size (ldb * (n-k)). To be accessed by columns.
  //            On input, it contains the remaining (n-k) columns/rows of M.
  //            On output, it contains the Schur complement.
  //            It can be null if k >= n.
  // - ldb    : Leading dimension of B.
  //            It must be at least (n-k), if k < n.
  //            It can be null if k >= n.
  // - info   : Information about outcome.
  //            info = 0, factorization successful
  //            info = i > 0, minor of order i not positive definite.
  //            info = -i < 0, i-th input argument illegal.
  // - k      : Number of columns/rows to factorize.
  //            If k < n, a partial factorization is computed.
  //            If k >= n, a full factorization is computed and B is not used.
  //
  //
  // Filippo Zanetti, 2024
  // ===========================================================================

  // ===========================================================================
  // Check input
  // ===========================================================================
  if (!uplo || !n || !A || !ldA || !info || !k || ((!B || !ldB) && *k < *n)) {
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
  if ((!upper && lda < max(1, na)) || (upper && lda < max(1, *k))) {
    printf("Invalid parameter lda\n");
    *info = -4;
    return;
  }
  if (ldB && *ldB < max(1, na - *k)) {
    printf("Invalid parameter ldb\n");
    *info = -6;
  }
  if (*k < 0) {
    printf("Invalid parameter k\n");
    *info = -8;
    return;
  }

  // Quick return
  if (*n == 0) return;

  const int nb = 64;

  // ===========================================================================
  // Main operations
  // ===========================================================================
  char Transt = 'T';
  char Transn = 'N';
  char Side = upper ? 'L' : 'R';
  char Unit = 'N';
  int N, M, K;
  double MinusOne = -1.0;
  double One = 1.0;

  if (upper) {
    // Compute A = U^T * U using BLAS level 3

    // j is the starting col of the block of columns
    for (int j = 0; j < *k; j += nb) {
      // jb is the size of the block
      int jb = min(nb, *k - j);

      // sizes for blas calls
      N = jb;
      K = j;
      M = na - j - jb;

      // update diagonal block
      dsyrk(uplo, &Transt, &N, &K, &MinusOne, &A[lda * j], ldA, &One,
            &A[j + lda * j], ldA);

      // factorize diagonal block
      H_dpotf2(uplo, &N, &A[j + lda * j], ldA, info);
      if (*info != 0) {
        *info += j - 1;
        return;
      }

      if (j + jb < na) {
        // update block of rows
        dgemm(&Transt, &Transn, &N, &M, &K, &MinusOne, &A[lda * j], ldA,
              &A[lda * (j + jb)], ldA, &One, &A[j + (j + jb) * lda], ldA);

        // solve block of columns with diagonal block
        dtrsm(&Side, uplo, &Transt, &Unit, &N, &M, &One, &A[j + lda * j], ldA,
              &A[j + (j + jb) * lda], ldA);
      }
    }
  } else {
    // Compute A = L * L^T using BLAS level 3

    // j is the starting col of the block of columns
    for (int j = 0; j < *k; j += nb) {
      // jb is the size of the block
      int jb = min(nb, *k - j);

      // sizes for blas calls
      N = jb;
      K = j;
      M = na - j - jb;

      // update diagonal block
      dsyrk(uplo, &Transn, &N, &K, &MinusOne, &A[j], ldA, &One, &A[j + lda * j],
            ldA);

      // factorize diagonal block
      H_dpotf2(uplo, &N, &A[j + lda * j], ldA, info);
      if (*info != 0) {
        *info += j - 1;
        return;
      }

      if (j + jb < na) {
        // update block of columns
        dgemm(&Transn, &Transt, &M, &N, &K, &MinusOne, &A[j + jb], ldA, &A[j],
              ldA, &One, &A[j + jb + lda * j], ldA);

        // solve block of columns with diagonal block
        dtrsm(&Side, uplo, &Transt, &Unit, &M, &N, &One, &A[j + lda * j], ldA,
              &A[j + jb + lda * j], ldA);
      }
    }
  }

  // update Schur complement if partial factorization is required
  if (*k < na) {
    N = na - *k;
    if (upper) {
      dsyrk(uplo, &Transt, &N, k, &MinusOne, &A[*k * lda], ldA, &One, B, ldB);
    } else {
      dsyrk(uplo, &Transn, &N, k, &MinusOne, &A[*k], ldA, &One, B, ldB);
    }
  }
}