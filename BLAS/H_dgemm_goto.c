#include "Hblas.h"

void H_dgemm(char* transa, char* transb, int* m, int* n, int* k, double* alpha,
             double* restrict A, int* ldA, double* restrict B, int* ldB,
             double* beta, double* restrict C, int* ldC) {
  // ===========================================================================
  // HiGHS BLAS function
  // DGEMM: Double GEneral Matrix Matrix product
  //
  // Version based on Goto approach
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
      }

      // max index for blocks
      const int K = (*k - 1) / kc;

      for (int j = 0; j <= K; ++j) {
        // size of last block may be smaller than kc
        int k_size = (j < K) ? kc : (*k % kc);
        if (k_size == 0) k_size = kc;

        H_dgepp(transa, transb, m, n, &k_size, alpha, &A[j * kc * lda], ldA,
                &B[j * kc], ldB, C, ldC);
      }

    } else {
      // operation [2]
      printf("Not supported\n");
    }
  } else {
    if (notransa) {
      // operation [3]
      printf("Not supported\n");
    } else {
      // operation [4]
      printf("Not supported\n");
    }
  }

  printf("Time copy A %5.2e\n", time_copy_A);
  printf("Time copy B %5.2e\n", time_copy_B);
  printf("Time kernel %5.2e\n", time_kernel);
}

void H_dgepp(char* transa, char* transb, int* m, int* n, int* k, double* alpha,
             double* restrict A, int* ldA, double* restrict B, int* ldB,
             double* restrict C, int* ldC) {
  // for now works only with transa = 'N', transb = 'N'

  double time = get_time();
  // allocate local copy of B in "special" contiguous format
  double* B_cont = malloc(sizeof(double) * (*n) * (*k));

  const int N = (*n - 1) / nr;

  // go through columns of B
  for (int j = 0; j < *n; ++j) {
    // in which block are we
    const int block = j / nr;

    // size of last block may be smaller than nr
    int n_size = (block < N) ? nr : (*n % nr);
    if (n_size == 0) n_size = nr;

    // go through the rows of B
    for (int i = 0; i < *k; ++i) {
      int pos = block * (*k) * nr + i * n_size + j % nr;
      B_cont[pos] = B[i + j * (*ldB)];
    }
  }
  time_copy_B += get_time() - time;

  // call dgebp
  const int M = (*m - 1) / mc;
  for (int j = 0; j <= M; ++j) {
    int m_size = (j < M) ? mc : (*m % mc);
    if (m_size == 0) m_size = mc;

    H_dgebp(transa, transb, &m_size, n, k, alpha, &A[j * mc], ldA, B_cont,
            &C[j * mc], ldC);
  }

  // free local copy of B
  free(B_cont);
}

void H_dgebp(char* transa, char* transb, int* m, int* n, int* k, double* alpha,
             double* restrict A, int* ldA, double* restrict B,
             double* restrict C, int* ldC) {
  // for now works only with transa = 'N', transb = 'N'

  double time = get_time();
  // allocate local copy of A in "special" contiguous format
  double* A_cont = malloc(sizeof(double) * (*m) * (*k));

  const int M = (*m - 1) / mr;

  // go through column of A
  for (int j = 0; j < *k; ++j) {
    for (int i = 0; i < *m; ++i) {
      // in which block are we
      const int block = i / mr;

      // size of last block may be smaller than mr
      int m_size = (block < M) ? mr : (*m % mr);
      if (m_size == 0) m_size = mr;

      int pos = block * (*k) * mr + j * m_size + i % mr;
      A_cont[pos] = A[i + j * (*ldA)];
    }
  }
  time_copy_A += get_time() - time;

  const int N = (*n - 1) / nr;
  // call kernel
  double* C_aux = calloc(*m * nr, sizeof(double));

  for (int j = 0; j <= N; ++j) {
    int n_size = (j < N) ? nr : (*n % nr);
    if (n_size == 0) n_size = nr;

    time = get_time();
    H_kernel(transa, transb, m, &n_size, k, alpha, A_cont, &B[*k * nr * j],
             C_aux, m);
    time_kernel += get_time() - time;

    // unpack C_aux
    for (int col = 0; col < n_size; ++col) {
      for (int row = 0; row < *m; ++row) {
        C[row + (col + nr * j) * *ldC] += C_aux[row + col * *m];
      }
    }

    memset(C_aux, 0, *m * nr * sizeof(double));
  }

  free(C_aux);
  free(A_cont);
}

void H_kernel(char* transa, char* transb, int* m, int* n, int* k, double* alpha,
              double* restrict A, double* restrict B, double* restrict C,
              int* ldC) {
  const int R = (*m - 1) / mr;

  double* C_aux = calloc(mr * *n, sizeof(double));

  for (int r = 0; r <= R; ++r) {
    int r_size = (r < R) ? mr : (*m % mr);
    if (r_size == 0) r_size = mr;

    for (int col = 0; col < *k; ++col) {
      for (int i = 0; i < r_size; ++i) {
        for (int j = 0; j < *n; ++j) {
          // C[r * mr + i + j * *ldC] +=
          C_aux[i + j * mr] +=
              *alpha * A[r * mr * *k + col * r_size + i] * B[*n * col + j];
        }
      }
    }

    // unpace C_aux into C

    for (int col = 0; col < *n; ++col) {
      for (int row = 0; row < r_size; ++row) {
        C[r * mr + row + col * *ldC] += C_aux[row + col * mr];
      }
    }

    memset(C_aux, 0, mr * *n * sizeof(double));
  }

  free(C_aux);
}