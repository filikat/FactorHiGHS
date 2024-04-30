#include "PartialFact.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "PartialFact_declaration.h"

// ===========================================================================
// Functions to compute dense partial Cholesky or LDL factorizations with
// left-looking approach, with or without blocking.
//
// A is used to access the first k columns, i.e., M(0:n-1,0:k-1).
// B is used to access the remaining lower triangle, i.e., M(k:n-1,k:n-1).
//
// NB: the content of B is discarded.
//
// For indefinite matrices:
// - 2x2 pivoting is not performed. If a zero pivot is found, the code stops.
// - the elements of D are stored as diagonal entries of A; the unit diagonal
//   entries of A are not stored.
//
// These functions are similar to Lapack dpotrf("L", n, A, lda, info)  and
// dpotf2("L", n, A, lda, info), if k >= n.
//
// Arguments:
// - n      : Dimension of matrix M.
// - k      : Number of columns to factorize.
//            If k < n, a partial factorization is computed.
//            If k >= n, a full factorization is computed and B is not used.
// - A      : Array of size (lda * k). To be accessed by columns.
//            On input, it contains the first k columns/rows of M.
//            On output, it contains the trapezoidal factor of the first k
//            columns of M.
// - lda    : Leading dimension of A. It must be at least n. It can be larger
//            if A is stored as a block of a larger matrix.
// - B      : Array of size (ldb * (n-k)). To be accessed by columns.
//            On input, any data is discarded.
//            On output, it contains the Schur complement.
//            It can be null if k >= n.
// - ldb    : Leading dimension of B. It must be at least (n-k), if k < n. It
//            can be larger if B is stored as a block of a larger martix.
//            Not used if k >= n.
//
// Return:
// - 0      : factorization successful
// - i > 0  : i-th pivot illegal
// - -1     : invalid input
//
// ===========================================================================

int PartialFactPosSmall(int n, int k, double* restrict A, int lda,
                        double* restrict B, int ldb) {
  // ===========================================================================
  // Positive definite factorization without blocks.
  // BLAS calls: ddot, dgemv, dscal, dsyrk.
  // ===========================================================================

  // check input
  if (n < 0 || k < 0 || !A || lda < n || (k < n && (!B || ldb < n - k))) {
    printf("Invalid input to PartialFactPosSmall\n");
    return -1;
  }

  // quick return
  if (n == 0) return 0;

  // main operations
  for (int j = 0; j < k; ++j) {
    int N = j;
    int M = n - j - 1;

    // update diagonal element
    double Ajj = A[j + lda * j] - ddot(&N, &A[j], &lda, &A[j], &lda);
    if (Ajj <= 0.0 || isnan(Ajj)) {
      A[j + lda * j] = Ajj;
      return j;
    }

    // compute diagonal element
    Ajj = sqrt(Ajj);
    A[j + lda * j] = Ajj;
    double coeff = 1.0 / Ajj;

    // compute column j
    if (j < n - 1) {
      dgemv(&NN, &M, &N, &d_m_one, &A[j + 1], &lda, &A[j], &lda, &d_one,
            &A[j + 1 + j * lda], &i_one);
      dscal(&M, &coeff, &A[j + 1 + j * lda], &i_one);
    }
  }

  // update Schur complement
  if (k < n) {
    int N = n - k;
    dsyrk(&LL, &NN, &N, &k, &d_m_one, &A[k], &lda, &d_zero, B, &ldb);
  }

  return 0;
}

int PartialFactPosLarge(int n, int k, double* restrict A, int lda,
                        double* restrict B, int ldb, double* times) {
  // ===========================================================================
  // Positive definite factorization with blocks.
  // BLAS calls: dsyrk, dgemm, dtrsm.
  // ===========================================================================

  // check input
  if (n < 0 || k < 0 || !A || lda < n || (k < n && (!B || ldb < n - k))) {
    printf("Invalid input to PartialFactPosLarge\n");
    return -1;
  }

  if (!times) {
    printf("Invalid times\n");
    return -1;
  }

  // quick return
  if (n == 0) return 0;

  // j is the starting col of the block of columns
  for (int j = 0; j < k; j += nb) {
    // jb is the size of the block
    int jb = min(nb, k - j);

    // sizes for blas calls
    int N = jb;
    int K = j;
    int M = n - j - jb;

    // update diagonal block
    double t0 = GetTime();
    dsyrk(&LL, &NN, &N, &K, &d_m_one, &A[j], &lda, &d_one, &A[j + lda * j],
          &lda);
    times[t_dsyrk] += GetTime() - t0;

    // factorize diagonal block
    t0 = GetTime();
    int info = PartialFactPosSmall(N, N, &A[j + lda * j], lda, NULL, 0);
    times[t_fact] += GetTime() - t0;
    if (info != 0) {
      return info + j - 1;
    }
    if (j + jb < n) {
      // update block of columns
      t0 = GetTime();
      dgemm(&NN, &TT, &M, &N, &K, &d_m_one, &A[j + jb], &lda, &A[j], &lda,
            &d_one, &A[j + jb + lda * j], &lda);
      times[t_dgemm] += GetTime() - t0;

      // solve block of columns with diagonal block
      t0 = GetTime();
      dtrsm(&RR, &LL, &TT, &NN, &M, &N, &d_one, &A[j + lda * j], &lda,
            &A[j + jb + lda * j], &lda);
      times[t_dtrsm] += GetTime() - t0;
    }
  }

  // update Schur complement if partial factorization is required
  if (k < n) {
    int N = n - k;
    double t0 = GetTime();
    dsyrk(&LL, &NN, &N, &k, &d_m_one, &A[k], &lda, &d_zero, B, &ldb);
    times[t_dsyrk] += GetTime() - t0;
  }

  return 0;
}

int PartialFactIndSmall(int n, int k, double* restrict A, int lda,
                        double* restrict B, int ldb) {
  // ===========================================================================
  // Infedinite factorization without blocks.
  // BLAS calls: ddot, dgemv, dscal, dcopy, dgemm.
  // ===========================================================================

  // check input
  if (n < 0 || k < 0 || !A || lda < n || (k < n && (!B || ldb < n - k))) {
    printf("Invalid input to PartialFactIndSmall\n");
    return -1;
  }

  // quick return
  if (n == 0) return 0;

  // main operations
  for (int j = 0; j < k; ++j) {
    int N = j;
    int M = n - j - 1;

    // create temporary copy of row j, multiplied by pivots
    double* temp = malloc(j * sizeof(double));
    for (int i = 0; i < j; ++i) {
      temp[i] = A[j + i * lda] * A[i + i * lda];
    }

    // update diagonal element
    double Ajj = A[j + lda * j] - ddot(&N, &A[j], &lda, temp, &i_one);
    if (Ajj == 0.0 || isnan(Ajj)) {
      A[j + lda * j] = Ajj;
      return j;
    }

    // save diagonal element
    A[j + lda * j] = Ajj;
    double coeff = 1.0 / Ajj;

    // compute column j
    if (j < n - 1) {
      dgemv(&NN, &M, &N, &d_m_one, &A[j + 1], &lda, temp, &i_one, &d_one,
            &A[j + 1 + j * lda], &i_one);
      dscal(&M, &coeff, &A[j + 1 + j * lda], &i_one);
    }

    // free temporary copy of row
    free(temp);
  }

  // update Schur complement
  if (k < n) {
    int N = n - k;

    // make temporary copy of M(k:n-1,0:k-1), multiplied by pivots
    double* temp = malloc((n - k) * k * sizeof(double));
    for (int j = 0; j < k; ++j) {
      dcopy(&N, &A[k + j * lda], &i_one, &temp[j * (n - k)], &i_one);
      dscal(&N, &A[j + j * lda], &temp[j * (n - k)], &i_one);
    }

    // update Schur complement using dgemm
    dgemm(&NN, &TT, &N, &N, &k, &d_m_one, &A[k], &lda, temp, &N, &d_zero, B,
          &ldb);

    free(temp);
  }

  return 0;
}

int PartialFactIndLarge(int n, int k, double* restrict A, int lda,
                        double* restrict B, int ldb, double* times) {
  // ===========================================================================
  // Indefinite factorization with blocks.
  // BLAS calls: dcopy, dscal, dgemm, dtrsm, dsyrk
  // ===========================================================================

  // check input
  if (n < 0 || k < 0 || !A || lda < n || (k < n && (!B || ldb < n - k))) {
    printf("Invalid input to PartialFactIndLarge\n");
    return -1;
  }

  // quick return
  if (n == 0) return 0;

  // j is the starting col of the block of columns
  for (int j = 0; j < k; j += nb) {
    // jb is the size of the block
    int jb = min(nb, k - j);

    // sizes for blas calls
    int N = jb;
    int K = j;
    int M = n - j - jb;

    // create temporary copy of block of rows, multiplied by pivots
    double* temp = malloc(j * jb * sizeof(double));
    int ldt = jb;
    for (int i = 0; i < j; ++i) {
      dcopy(&N, &A[j + i * lda], &i_one, &temp[i * ldt], &i_one);
      dscal(&N, &A[i + i * lda], &temp[i * ldt], &i_one);
    }

    // update diagonal block using dgemm
    double t0 = GetTime();
    dgemm(&NN, &TT, &jb, &jb, &j, &d_m_one, &A[j], &lda, temp, &ldt, &d_one,
          &A[j + lda * j], &lda);
    times[t_dgemm] += GetTime() - t0;

    // factorize diagonal block
    t0 = GetTime();
    int info = PartialFactIndSmall(N, N, &A[j + lda * j], lda, NULL, 0);
    times[t_fact] += GetTime() - t0;
    if (info != 0) {
      return info + j - 1;
    }

    if (j + jb < n) {
      // update block of columns
      t0 = GetTime();
      dgemm(&NN, &TT, &M, &N, &K, &d_m_one, &A[j + jb], &lda, temp, &ldt,
            &d_one, &A[j + jb + lda * j], &lda);
      times[t_dgemm] += GetTime() - t0;

      // solve block of columns with L
      t0 = GetTime();
      dtrsm(&RR, &LL, &TT, &UU, &M, &N, &d_one, &A[j + lda * j], &lda,
            &A[j + jb + lda * j], &lda);
      times[t_dtrsm] += GetTime() - t0;

      // solve block of columns with D
      for (int i = 0; i < jb; ++i) {
        double coeff = 1.0 / A[j + i + (j + i) * lda];
        dscal(&M, &coeff, &A[j + jb + lda * (j + i)], &i_one);
      }
    }

    free(temp);
  }

  // update Schur complement
  if (k < n) {
    int N = n - k;

    // count number of positive and negative pivots
    int pos_pivot = 0;
    int neg_pivot = 0;
    for (int i = 0; i < k; ++i) {
      if (A[i + lda * i] >= 0.0) {
        ++pos_pivot;
      } else {
        ++neg_pivot;
      }
    }

    // make temporary copies of positive and negative columns separately
    double* temp_pos = malloc((n - k) * pos_pivot * sizeof(double));
    double* temp_neg = malloc((n - k) * neg_pivot * sizeof(double));
    int ldt = n - k;

    // the copies of the columns are multiplied by sqrt(|Ajj|)
    int start_pos = 0;
    int start_neg = 0;
    for (int j = 0; j < k; ++j) {
      double Ajj = A[j + lda * j];
      if (Ajj >= 0.0) {
        Ajj = sqrt(Ajj);
        dcopy(&N, &A[k + j * lda], &i_one, &temp_pos[start_pos * ldt], &i_one);
        dscal(&N, &Ajj, &temp_pos[start_pos * ldt], &i_one);
        ++start_pos;
      } else {
        Ajj = sqrt(-Ajj);
        dcopy(&N, &A[k + j * lda], &i_one, &temp_neg[start_neg * ldt], &i_one);
        dscal(&N, &Ajj, &temp_neg[start_neg * ldt], &i_one);
        ++start_neg;
      }
    }

    // Update schur complement by subtracting contribution of positive columns
    // and adding contribution of negative columns.
    // In this way, I can use dsyrk instead of dgemm and avoid updating the full
    // square schur complement.
    // First call uses beta = 0.0, to clear content of B.
    // Second call uses beta = 1.0, to not clear the result of the first call.
    double t0 = GetTime();
    dsyrk(&LL, &NN, &N, &pos_pivot, &d_m_one, temp_pos, &ldt, &d_zero, B, &ldb);
    dsyrk(&LL, &NN, &N, &neg_pivot, &d_one, temp_neg, &ldt, &d_one, B, &ldb);
    times[t_dsyrk] += GetTime() - t0;

    free(temp_pos);
    free(temp_neg);
  }

  return 0;
}
