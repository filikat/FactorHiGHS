#include "Hblas.h"

void H_dgemm(char* transa, char* transb, int* m, int* n, int* k, double* alpha,
             double* restrict A, int* ldA, double* restrict B, int* ldB,
             double* beta, double* restrict C, int* ldC) {
  // ===========================================================================
  // HiGHS BLAS function
  // DGEMM: Double GEneral Matrix Matrix product
  //
  // Perform one of:
  //  C = beta * C + alpha * A * B       [1]
  //  C = beta * C + alpha * A^T * B     [2]
  //  C = beta * C + alpha * A * B^T     [3]
  //  C = beta * C + alpha * A^T * B^T   [4]
  //
  // Arguments:
  // - transa : perform operation with A ('N','n') or A^T ('T','t','C','c').
  // - transb : perform operation with B ('N','n') or B^T ('T','t','C','c').
  // - m      : number of rows of C.
  // - n      : number of columns of C.
  // - k      : number of columns of A, for [1],[3].
  //            number of rows of A, for [2],[4].
  //            number of columns of B, for [3],[4].
  //            number of rows of B, for [1],[2].
  // - alpha  : scalar.
  // - A      : array of size (lda * k), for operation [1],[3];
  //            array of size (lda * m), for operation [2],[4].
  //            To be accessed by columns.
  // - lda    : leading dimension of A.
  // - B      : array of size (lda * n), for operation [1],[2];
  //            array of size (lda * k), for operation [3],[4].
  //            To be accessed by columns.
  // - ldb    : leading dimension of B.
  // - beta   : scalar.
  // - C      : array of size (ldc * n), to be accessed by columns.
  // - ldc    : leading dimension of C.
  //
  //
  // Filippo Zanetti, 2024
  // ===========================================================================

  // ===========================================================================
  // Check input
  // ===========================================================================
  if (!transa || !transb || !m || !n || !k || !alpha || !A || !ldA || !B ||
      !ldB || !beta || !C || !ldC) {
    printf("Invalid pointer\n");
    return;
  }

  int notransa = (*transa == 'n' || *transa == 'N');
  int notransb = (*transb == 'n' || *transb == 'N');
  int nrowa = notransa ? *m : *k;
  int nrowb = notransb ? *k : *n;
  int nrowc = *m;
  int ncolc = *n;
  int lda = *ldA;
  int ldb = *ldB;
  int ldc = *ldC;

  if (!notransa && !(*transa == 't' || *transa == 'T') &&
      !(*transa == 'c' || *transa == 'C')) {
    printf("Invalid parameter transa\n");
    return;
  }
  if (!notransb && !(*transb == 't' || *transb == 'T') &&
      !(*transb == 'c' || *transb == 'C')) {
    printf("Invalid parameter transb\n");
    return;
  }
  if (*m < 0) {
    printf("Invalid parameter m\n");
    return;
  }
  if (*n < 0) {
    printf("Invalid parameter n\n");
    return;
  }
  if (*k < 0) {
    printf("Invalid parameter k\n");
    return;
  }
  if (lda < max(1, nrowa)) {
    printf("Invalid parameter lda\n");
    return;
  }
  if (ldb < max(1, nrowb)) {
    printf("Invalid parameter ldb\n");
    return;
  }
  if (ldc < max(1, nrowc)) {
    printf("Invalid parameter ldc\n");
    return;
  }

  // Quick return
  if (*m == 0 || *n == 0 || ((*alpha == 0.0 || *k == 0) && *beta == 1.0))
    return;

  // ===========================================================================
  // Case alpha = 0
  // ===========================================================================
  if (*alpha == 0.0) {
    if (*beta == 0.0)
      for (int j = 0; j < nrowc * ncolc; ++j) C[j] = 0.0;
    else
      for (int j = 0; j < nrowc * ncolc; ++j) C[j] *= *beta;
    return;
  }

  // ===========================================================================
  // General case
  // ===========================================================================
  double temp;

  if (notransb) {
    if (notransa) {
      // operation [1]
      for (int j = 0; j < ncolc; ++j) {
        // Prepare C based on beta
        if (*beta == 0.0) {
          for (int i = 0; i < nrowc; ++i) C[i + ldc * j] = 0.0;
        } else if (*beta != 1.0) {
          for (int i = 0; i < nrowc; ++i) C[i + ldc * j] *= *beta;
        }

        // Compute values
        for (int l = 0; l < *k; ++l) {
          temp = *alpha * B[l + ldb * j];
          for (int i = 0; i < nrowc; ++i)
            C[i + ldc * j] += temp * A[i + lda * l];
        }
      }
    } else {
      // operation [2]
      for (int j = 0; j < ncolc; ++j) {
        for (int i = 0; i < nrowc; ++i) {
          // Compute values
          temp = 0.0;
          for (int l = 0; l < *k; ++l) temp += A[l + lda * i] * B[l + ldb * j];

          // Store in C
          if (*beta == 0.0)
            C[i + ldc * j] = *alpha * temp;
          else
            C[i + ldc * j] = *alpha * temp + *beta * C[i + ldc * j];
        }
      }
    }
  } else {
    if (notransa) {
      // operation [3]
      for (int j = 0; j < ncolc; ++j) {
        // Prepare C based on beta
        if (*beta == 0.0) {
          for (int i = 0; i < nrowc; ++i) C[i + ldc * j] = 0.0;
        } else if (*beta != 1.0) {
          for (int i = 0; i < nrowc; ++i) C[i + ldc * j] *= *beta;
        }

        // Compute values
        for (int l = 0; l < *k; ++l) {
          temp = *alpha * B[j + ldb * l];
          for (int i = 0; i < nrowc; ++i)
            C[i + ldc * j] += temp * A[i + lda * l];
        }
      }
    } else {
      // operation [4]
      for (int j = 0; j < ncolc; ++j) {
        for (int i = 0; i < nrowc; ++i) {
          // Compute values
          temp = 0.0;
          for (int l = 0; l < *k; ++l) temp += A[l + lda * i] * B[j + ldb * l];

          // Store in C
          if (*beta == 0.0)
            C[i + ldc * j] = *alpha * temp;
          else
            C[i + ldc * j] = *alpha * temp + *beta * C[i + ldc * j];
        }
      }
    }
  }
}