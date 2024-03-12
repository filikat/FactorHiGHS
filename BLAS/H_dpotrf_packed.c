#include "Hblas.h"

void dsyrk(char* uplo, char* trans, int* n, int* k, double* alpha, double* a,
           int* lda, double* beta, double* c, int* ldc);

void dtrsm(char* side, char* uplo, char* trans, char* diag, int* m, int* n,
           double* alpha, double* a, int* lda, double* b, int* ldb);

void dgemm(char* transa, char* transb, int* m, int* n, int* k, double* alpha,
           double* A, int* lda, double* B, int* ldb, double* beta, double* C,
           int* ldc);

void dpotrf(char* uplo, int* n, double* A, int* ldA, int* info);

void H_dpotf2(char* uplo, int* n, double* A, int* ldA, int* info);

void printDiag(double* A, int jb) {
  int next = 0;
  for (int i = 0; i < jb; ++i) {
    int j = 0;
    for (; j <= i; ++j) {
      printf("%9.3f ", A[next++]);
    }
    for (; j < jb; ++j) {
      printf("%9.3f ", 0.0);
    }
    printf("\n");
  }
  printf("\n\n");
}

void printFullRow(double* A, int ncol, int nrow) {
  int next = 0;
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      printf("%7.3f ", A[next++]);
    }
    printf("\n");
  }
  printf("\n\n");
}

void H_dpotrf_packed(int* N, double* A, int* NB) {
  // ===========================================================================
  // HiGHS BLAS function
  // DPOTRF: Double POsitive definite TRiangular Factorization
  // This version calls level 3 BLAS functions
  // Packed version
  //
  // Perform one of:
  //  A = L * L^T
  // in packed format.
  //
  // Arguments:
  // - n      : dimension of matrix A.
  // - A      : array of size (lda * n). To be accessed by columns.
  // - nB     : size of the blocks.
  //
  //
  // Filippo Zanetti, 2024
  // ===========================================================================

  // ===========================================================================
  // Check input
  // ===========================================================================
  if (!N || !A || !NB) {
    printf("Invalid pointer\n");
    return;
  }

  int n = *N;
  int nb = *NB;

  if (n < 0) {
    printf("Invalid parameter n\n");
    return;
  }
  if (nb < 0) {
    printf("Invalid parameter nb\n");
    return;
  }

  // Quick return
  if (n == 0) return;

  // number of blocks
  int n_blocks = (n - 1) / nb + 1;

  // start of diagonal blocks
  int* start = malloc(n_blocks * sizeof(double));
  for (int i = 0; i < n_blocks; ++i) {
    start[i] = i * nb * (2 * n + 1 - i * nb) / 2;
  }

  // size of blocks
  int diagSize = nb * (nb + 1) / 2;
  int fullSize = nb * nb;

  char up = 'U';
  char transt = 'T';
  char transn = 'N';
  char side = 'L';
  char unit = 'N';
  double MinusOne = -1.0;
  double One = 1.0;
  int info;

  for (int j = 0; j < n_blocks; ++j) {
    // jb is the size of the block
    int jb = min(nb, n - nb * j);

    // full copy of diagonal block by rows
    double* D = malloc(jb * jb * sizeof(double));
    int offset = 0;
    for (int ii = 0; ii < jb; ++ii) {
      for (int jj = 0; jj <= ii; ++jj) {
        D[ii * jb + jj] = A[start[j] + offset++];
      }
    }

    int M = n - nb * (j + 1);
    int Lij_pos = start[j] + diagSize;

    // update diagonal block and block of columns
    for (int k = 0; k < j; ++k) {
      int Ljk_pos = start[k] + diagSize;
      if (j > k + 1) Ljk_pos += fullSize * (j - k - 1);

      dsyrk(&up, &transt, &jb, &nb, &MinusOne, &A[Ljk_pos], &nb, &One, D, &jb);

      if (j != n_blocks - 1) {
        int Lik_pos = Ljk_pos + fullSize;
        dgemm(&transt, &transn, &jb, &M, &nb, &MinusOne, &A[Ljk_pos], &nb,
              &A[Lik_pos], &nb, &One, &A[Lij_pos], &jb);
      }
    }

    // factorize diagonal block
    dpotrf(&up, &jb, D, &jb, &info);
    if (info != 0) {
      printf("Wrong info\n");
      return;
    }

    if (j != n_blocks - 1) {
      // solve block of columns with diagonal block
      dtrsm(&side, &up, &transt, &unit, &jb, &M, &One, D, &jb, &A[Lij_pos],
            &jb);
    }

    // put D back into packed format
    offset = 0;
    for (int ii = 0; ii < jb; ++ii) {
      for (int jj = 0; jj <= ii; ++jj) {
        A[start[j] + offset++] = D[ii * jb + jj];
      }
    }

    free(D);
  }

  free(start);
}