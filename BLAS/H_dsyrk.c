#include "Hblas.h"

void H_dsyrk(char* uplo, char* trans, int* n, int* k, double* alpha,
             double* restrict A, int* ldA, double* beta, double* restrict C,
             int* ldC) {
  // ===========================================================================
  // HiGHS BLAS function
  // DSYRK: Double SYmmetric Rank K update
  //
  // Perform one of:
  //  C = beta * C + alpha * A * A^T   [1]
  //  C = beta * C + alpha * A^T * A   [2]
  //
  // Arguments:
  // - uplo  : upper('U','u') or lower('L','l') triangle of C is used.
  // - trans : perform operation [1] if trans = 'N','n';
  //           perform operation [2] if trans = 'T','t','C','c'.
  // - n     : size of square matrix C.
  // - k     : number of columns of A, for operation [1];
  //           number of rows of A, for operation [2].
  // - alpha : scalar.
  // - A     : array of size (lda * k), for operation [1];
  //           array of size (lda * n), for operation [2].
  //           To be accessed by columns.
  // - lda   : leading dimension of A.
  // - beta  : scalar.
  // - C     : array of size (ldc * n), to be accessed by columns.
  // - ldc   : leading dimension of C.
  //
  //
  // Filippo Zanetti, 2024
  // ===========================================================================

  // ===========================================================================
  // Check input
  // ===========================================================================
  if (!uplo || !trans || !n || !k || !alpha || !A || !ldA || !beta || !C ||
      !ldC) {
    printf("Invalid pointer\n");
    return;
  }

  int upper = (*uplo == 'u' || *uplo == 'U');
  int notrans = (*trans == 'n' || *trans == 'N');
  int nrowc = *n;
  int nrowa = notrans ? *n : *k;
  int ldc = *ldC;
  int lda = *ldA;

  if (!upper && !(*uplo == 'l' || *uplo == 'L')) {
    printf("Invalid parameter uplo\n");
    return;
  }
  if (!notrans && !(*trans == 't' || *trans == 'T') &&
      !(*trans == 'c' || *trans == 'C')) {
    printf("Invalid parameter trans\n");
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
  if (ldc < max(1, nrowc)) {
    printf("Invalid parameter ldc\n");
    return;
  }

  // Quick return
  if (*n == 0 || ((*alpha == 0.0 || *k == 0) && *beta == 1.0)) return;

  // ===========================================================================
  // Remark
  // ===========================================================================
  // A and C are stored by columns.
  // To access A(i,j), i.e., entry in row i, column j, the operation is
  //  A[i + lda * j].
  //
  // To perform the operations, there are three loops:
  // - for i, loops over the rows of C
  // - for j, loops over the columns of C
  // - for l, loops over the rows/columns of A (depending on trans)
  //
  // The main operation to compute C = A * A^T or C = A^T * A, in the upper
  // case, is:
  //
  //  for j = 0,...,n-1
  //    for i = 0,...,j
  //      for l = 0,...,k-1
  //        C(i,j) = A(i,l) * A(j,l), case C = A * A^T
  //        C(i,j) = A(l,i) * A(l,j), case C = A^T * A
  //
  // To guarantee cache efficiency, these loops have to be done in the order
  //  j,l,i
  // when computing C = A * A^T, and in the order
  //  j,i,l
  // when computing C = A^T * A.

  // ===========================================================================
  // Case alpha = 0
  // ===========================================================================
  if (*alpha == 0.0) {
    if (upper) {
      if (*beta == 0.0) {
        for (int j = 0; j < nrowc; ++j) {
          for (int i = 0; i < j + 1; ++i) C[ldc * j + i] = 0.0;
        }
      } else {
        for (int j = 0; j < nrowc; ++j) {
          for (int i = 0; i < j + 1; ++i) C[ldc * j + i] *= *beta;
        }
      }
    } else {
      if (*beta == 0.0) {
        for (int j = 0; j < nrowc; ++j) {
          for (int i = j; i < nrowc; ++i) C[ldc * j + i] = 0.0;
        }
      } else {
        for (int j = 0; j < nrowc; ++j) {
          for (int i = j; i < nrowc; ++i) C[ldc * j + i] *= *beta;
        }
      }
    }
    return;
  }

  // ===========================================================================
  // General case
  // ===========================================================================
  double temp;

  if (notrans) {
    // Case C = beta * C + alpha * A * A^T
    // Loops in order jli for locality of reference

    if (upper) {
      for (int j = 0; j < nrowc; ++j) {
        // Prepare C according to beta
        if (*beta == 0.0) {
          for (int i = 0; i < j + 1; ++i) C[ldc * j + i] = 0.0;

        } else if (*beta != 1.0) {
          for (int i = 0; i < j + 1; ++i) C[ldc * j + i] *= *beta;
        }

        // Compute values
        for (int l = 0; l < *k; ++l) {
          if (A[lda * l + j] != 0.0) {
            temp = *alpha * A[lda * l + j];
            for (int i = 0; i < j + 1; ++i)
              C[ldc * j + i] += temp * A[lda * l + i];
          }
        }
      }
    } else {
      for (int j = 0; j < nrowc; ++j) {
        // Prepare C according to beta
        if (*beta == 0.0) {
          for (int i = j; i < nrowc; ++i) C[ldc * j + i] = 0.0;

        } else if (*beta != 1.0) {
          for (int i = j; i < nrowc; ++i) C[ldc * j + i] *= *beta;
        }

        // Compute values
        for (int l = 0; l < *k; ++l) {
          if (A[lda * l + j] != 0.0) {
            temp = *alpha * A[lda * l + j];
            for (int i = j; i < nrowc; ++i)
              C[ldc * j + i] += temp * A[lda * l + i];
          }
        }
      }
    }
  } else {
    // Case C = beta * C + alpha * A^T * A
    // Loops in order jil for locality of reference

    if (upper) {
      for (int j = 0; j < nrowc; ++j) {
        for (int i = 0; i < j + 1; ++i) {
          // Compute value
          temp = 0.0;
          for (int l = 0; l < *k; ++l) temp += A[lda * i + l] * A[lda * j + l];

          // Store in C
          if (*beta == 0.0)
            C[ldc * j + i] = *alpha * temp;
          else
            C[ldc * j + i] = *alpha * temp + *beta * C[ldc * j + i];
        }
      }
    } else {
      for (int j = 0; j < nrowc; ++j) {
        for (int i = j; i < nrowc; ++i) {
          // Compute value
          temp = 0.0;
          for (int l = 0; l < *k; ++l) temp += A[lda * i + l] * A[lda * j + l];

          // Store in C
          if (*beta == 0.0)
            C[ldc * j + i] = *alpha * temp;
          else
            C[ldc * j + i] = *alpha * temp + *beta * C[ldc * j + i];
        }
      }
    }
  }
}