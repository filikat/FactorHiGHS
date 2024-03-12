#include "Hblas.h"

double ddot(int* n, double* dx, int* incx, double* dy, int* incy);
void dgemv(char* trans, int* m, int* n, double* alpha, double* A, int* lda,
           double* x, int* incx, double* beta, double* y, int* incy);
void dscal(int* n, double* da, double* dx, int* incx);

void H_dpotf2(char* uplo, int* n, double* restrict A, int* ldA, int* info) {
  // ===========================================================================
  // HiGHS BLAS function
  // DPOTF2: Double POsitive definite Triangular Factorization
  // This version calls level 2 BLAS functions
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
  //
  //
  // Filippo Zanetti, 2024
  // ===========================================================================

  // ===========================================================================
  // Check input
  // ===========================================================================
  if (!uplo || !n || !A || !ldA || !info) {
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

  // Quick return
  if (*n == 0) return;

  // ===========================================================================
  // General case
  // ===========================================================================

  if (upper) {
    // Compute A = U^T * U using BLAS level 2

    int N, M;
    char Trans = 'T';
    double One = 1.0;
    double MinusOne = -1.0;
    int IncOne = 1;
    double coeff;

    for (int j = 0; j < na; ++j) {
      N = j;
      M = na - j - 1;

      // compute diagonal element
      double Ajj =
          A[j + lda * j] - ddot(&N, &A[lda * j], &IncOne, &A[lda * j], &IncOne);
      if (Ajj <= 0.0 || isnan(Ajj)) {
        A[j + lda * j] = Ajj;
        *info = j;
        return;
      }
      Ajj = sqrt(Ajj);
      A[j + lda * j] = Ajj;
      coeff = 1.0 / Ajj;

      // compute row j
      if (j < na - 1) {
        dgemv(&Trans, &N, &M, &MinusOne, &A[lda * (j + 1)], ldA, &A[lda * j],
              &IncOne, &One, &A[j + (j + 1) * lda], ldA);
        dscal(&M, &coeff, &A[j + (j + 1) * lda], ldA);
      }
    }

  } else {
    // Compute A = L * L^T using BLAS level 2

    int N, M;
    char Trans = 'N';
    double One = 1.0;
    double MinusOne = -1.0;
    int IncOne = 1;
    double coeff;

    for (int j = 0; j < na; ++j) {
      N = j;
      M = na - j - 1;

      // compute diagonal element
      double Ajj = A[j + lda * j] - ddot(&N, &A[j], ldA, &A[j], ldA);
      if (Ajj <= 0.0 || isnan(Ajj)) {
        A[j + lda * j] = Ajj;
        *info = j;
        return;
      }
      Ajj = sqrt(Ajj);
      A[j + lda * j] = Ajj;
      coeff = 1.0 / Ajj;

      // compute column j
      if (j < na - 1) {
        dgemv(&Trans, &M, &N, &MinusOne, &A[j + 1], ldA, &A[j], ldA, &One,
              &A[j + 1 + j * lda], &IncOne);
        dscal(&M, &coeff, &A[j + 1 + j * lda], &IncOne);
      }
    }
  }
}