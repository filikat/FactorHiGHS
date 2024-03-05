#include "Hblas.h"

void H_dtrsm(char* side, char* uplo, char* trans, char* diag, int* m, int* n,
             double* alpha, double* A, int* lda, double* B, int* ldb) {
  // ===========================================================================
  // HiGHS BLAS function
  // DTRSM: Double TRiangular Solve with Matrix
  //
  // Perform one of:
  //  B = alpha * A^-1 * B  [1],   for side = 'L' and trans = 'N'
  //  B = alpha * A^-T * B  [2],   for side = 'L' and trans = 'T'
  //  B = alpha * B * A^-1  [3],   for side = 'R' and trans = 'N'
  //  B = alpha * B * A^-T  [4],   for side = 'R' and trans = 'T'
  //
  // Arguments:
  // - side  : matrix A is applied to the left('L','l') or right('R','r') of B.
  // - uplo  : upper('U','u') or lower('L','l') triangle of A is used.
  // - trans : A is transposed if trans = 'N','n';
  //           A is not transposed if trans = 'T','t','C','c'.
  // - diag  : A has unit diagonal ('U','u') or not ('N','n').
  // - m     : number of rows of B.
  // - n     : number of columns of B.
  // - alpha : scalar.
  // - A     : array of size (lda * m), for operation [1] or [2];
  //           array of size (lda * n), for operation [3] or [4].
  //           To be accessed by columns.
  // - lda   : number of rows of A.
  // - B     : array of size (ldb * n), to be accessed by columns.
  // - ldb   : number of rows of B.
  //
  //
  // Filippo Zanetti, 2024
  // ===========================================================================

  // ===========================================================================
  // Check input
  // ===========================================================================
  if (!side || !uplo || !trans || !diag || !m || !n || !alpha || !A || !lda ||
      !B || !ldb) {
    printf("Invalid pointer\n");
    return;
  }

  bool left = (*side == 'l' || *side == 'L');
  bool upper = (*uplo == 'u' || *uplo == 'U');
  bool notrans = (*trans == 'n' || *trans == 'N');
  bool unit = (*diag == 'u' || *diag == 'U');
  int nrowa = left ? *m : *n;
  int nrowb = *m;
  int ncolb = *n;
  if (!left && !(*side == 'r' || *side == 'R')) {
    printf("Invalid parameter side\n");
    return;
  }
  if (!upper && !(*uplo == 'l' || *uplo == 'L')) {
    printf("Invalid parameter uplo\n");
    return;
  }
  if (!notrans && !(*trans == 't' || *trans == 'T') &&
      !(*trans == 'c' || *trans == 'C')) {
    printf("Invalid parameter trans\n");
    return;
  }
  if (!unit && !(*diag == 'n' || *diag == 'N')) {
    printf("Invalid parameter diag\n");
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
  if (*lda < std::max(1, nrowa)) {
    printf("Invalid parameter lda\n");
    return;
  }
  if (*ldb < std::max(1, nrowb)) {
    printf("Invalid parameter ldb\n");
    return;
  }

  // Quick return
  if (*n == 0 || *m == 0) return;

  // ===========================================================================
  // Case alpha = 0
  // ===========================================================================
  if (*alpha == 0.0) {
    for (int j = 0; j < ncolb; ++j) {
      for (int i = 0; i < nrowb; ++i) {
        B[i + nrowb * j] = 0.0;
      }
    }
    return;
  }

  // ===========================================================================
  // General case
  // ===========================================================================
  double temp;

  if (left) {
    if (notrans) {
      // operation [1]
      if (upper) {
        for (int j = 0; j < ncolb; ++j) {
          if (*alpha != 1.0) {
            for (int i = 0; i < nrowb; ++i) B[i + nrowb * j] *= *alpha;
          }
          for (int k = nrowb - 1; k >= 0; --k) {
            if (B[k + nrowb * j] != 0.0) {
              if (!unit) B[k + nrowb * j] /= A[k + nrowa * k];
              for (int i = 0; i < k; ++i)
                B[i + nrowb * j] -= B[k + nrowb * j] * A[i + nrowa * k];
            }
          }
        }
      } else {
        for (int j = 0; j < ncolb; ++j) {
          if (*alpha != 1.0) {
            for (int i = 0; i < nrowb; ++i) B[i + nrowb * j] *= *alpha;
          }
          for (int k = 0; k < nrowb; ++k) {
            if (B[k + nrowb * j] != 0.0) {
              if (!unit) B[k + nrowb * j] /= A[k + nrowa * k];
              for (int i = k + 1; i < nrowb; ++i)
                B[i + nrowb * j] -= B[k + nrowb * j] * A[i + nrowa * k];
            }
          }
        }
      }
    } else {
      // operation [2]
      if (upper) {
        for (int j = 0; j < ncolb; ++j) {
          for (int i = 0; i < nrowb; ++i) {
            temp = *alpha * B[i + nrowb * j];
            for (int k = 0; k < i; ++k)
              temp -= A[k + nrowa * i] * B[k + nrowb * j];
            if (!unit) temp /= A[i + nrowa * i];
            B[i + nrowb * j] = temp;
          }
        }
      } else {
        for (int j = 0; j < ncolb; ++j) {
          for (int i = nrowb - 1; i >= 0; --i) {
            temp = *alpha * B[i + nrowb * j];
            for (int k = i + 1; k < nrowb; ++k)
              temp -= A[k + nrowa * i] * B[k + nrowb * j];
            if (!unit) temp /= A[i + nrowa * i];
            B[i + nrowb * j] = temp;
          }
        }
      }
    }
  } else {
    if (notrans) {
      // operation [3]
      if (upper) {
        for (int j = 0; j < ncolb; ++j) {
          if (*alpha != 1.0) {
            for (int i = 0; i < nrowb; ++i) B[i + nrowb * j] *= *alpha;
          }
          for (int k = 0; k < j; ++k) {
            if (A[k + nrowa * j] != 0.0) {
              for (int i = 0; i < nrowb; ++i)
                B[i + nrowb * j] -= A[k + nrowa * j] * B[i + nrowb * k];
            }
          }
          if (!unit) {
            temp = 1.0 / A[j + nrowa * j];
            for (int i = 0; i < nrowb; ++i) B[i + nrowb * j] *= temp;
          }
        }
      } else {
        for (int j = ncolb - 1; j >= 0; --j) {
          if (*alpha != 1.0) {
            for (int i = 0; i < nrowb; ++i) B[i + nrowb * j] *= *alpha;
          }
          for (int k = j + 1; k < ncolb; ++k) {
            if (A[k + nrowa * j] != 0.0) {
              for (int i = 0; i < nrowb; ++i)
                B[i + nrowb * j] -= A[k + nrowa * j] * B[i + nrowb * k];
            }
          }
          if (!unit) {
            temp = 1.0 / A[j + nrowa * j];
            for (int i = 0; i < nrowb; ++i) B[i + nrowb * j] *= temp;
          }
        }
      }
    } else {
      // operation [4]
      if (upper) {
        for (int k = ncolb - 1; k >= 0; --k) {
          if (!unit) {
            temp = 1.0 / A[k + nrowa * k];
            for (int i = 0; i < nrowb; ++i) B[i + nrowb * k] *= temp;
          }
          for (int j = 0; j < k; ++j) {
            if (A[j + nrowa * k] != 0.0) {
              temp = A[j + nrowa * k];
              for (int i = 0; i < nrowb; ++i)
                B[i + nrowb * j] -= temp * B[i + nrowb * k];
            }
          }
          if (*alpha != 1.0) {
            for (int i = 0; i < nrowb; ++i) B[i + nrowb * k] *= *alpha;
          }
        }
      } else {
        for (int k = 0; k < ncolb; ++k) {
          if (!unit) {
            temp = 1.0 / A[k + nrowa * k];
            for (int i = 0; i < nrowb; ++i) B[i + nrowb * k] *= temp;
          }
          for (int j = k + 1; j < ncolb; ++j) {
            if (A[j + nrowa * k] != 0.0) {
              temp = A[j + nrowa * k];
              for (int i = 0; i < nrowb; ++i)
                B[i + nrowb * j] -= temp * B[i + nrowb * k];
            }
          }
          if (*alpha != 1.0) {
            for (int i = 0; i < nrowb; ++i) B[i + nrowb * k] *= *alpha;
          }
        }
      }
    }
  }
}