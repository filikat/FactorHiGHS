#include "DenseFact.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "DenseFact_declaration.h"

// #define TIMING

/*
Names:
DenseFact_(pf)(di)(bu)(flh)

pf: Partial or Full factorization
di: (positive) Definite or Indefinite
bu: Blocked or Unblocked
flh: Full format, Lower packed format, or lower-blocked-Hybrid packed format

*/

double t0;

int DenseFact_fduf(char uplo, int n, double* restrict A, int lda) {
  // ===========================================================================
  // Positive definite factorization without blocks.
  // BLAS calls: ddot, dgemv, dscal.
  // ===========================================================================

  // check input
  if (n < 0 || !A || lda < n || (uplo != 'L' && uplo != 'U')) {
    printf("DenseFact_fduf: invalid input\n");
    return ret_invalid_input;
  }

  // quick return
  if (n == 0) return ret_ok;

  // main operations
  if (uplo == 'L') {
    for (int j = 0; j < n; ++j) {
      int N = j;
      int M = n - j - 1;

      // update diagonal element
      double Ajj = A[j + lda * j] - ddot(&N, &A[j], &lda, &A[j], &lda);
      if (Ajj <= 0.0 || isnan(Ajj)) {
        A[j + lda * j] = Ajj;
        printf("DenseFact_fduf: invalid pivot\n");
        return ret_invalid_pivot;
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
  } else {
    for (int j = 0; j < n; ++j) {
      int N = j;
      int M = n - j - 1;

      // update diagonal element
      double Ajj =
          A[j + lda * j] - ddot(&N, &A[lda * j], &i_one, &A[lda * j], &i_one);
      if (Ajj <= 0.0 || isnan(Ajj)) {
        A[j + lda * j] = Ajj;
        printf("DenseFact_fduf: invalid pivot\n");
        return ret_invalid_pivot;
      }

      // compute diagonal element
      Ajj = sqrt(Ajj);
      A[j + lda * j] = Ajj;
      double coeff = 1.0 / Ajj;

      // compute column j
      if (j < n - 1) {
        dgemv(&TT, &N, &M, &d_m_one, &A[lda * (j + 1)], &lda, &A[lda * j],
              &i_one, &d_one, &A[j + (j + 1) * lda], &lda);
        dscal(&M, &coeff, &A[j + (j + 1) * lda], &lda);
      }
    }
  }

  return ret_ok;
}

int DenseFact_fiuf(char uplo, int n, double* restrict A, int lda) {
  // ===========================================================================
  // Infedinite factorization without blocks.
  // BLAS calls: ddot, dgemv, dscal.
  // ===========================================================================

  // check input
  if (n < 0 || !A || lda < n) {
    printf("DenseFact_fiuf: invalid input\n");
    return ret_invalid_input;
  }

  // quick return
  if (n == 0) return ret_ok;

  // main operations
  if (uplo == 'L') {
    // allocate space for copy of col multiplied by pivots
    double* temp = malloc((n - 1) * sizeof(double));
    if (!temp) {
      printf("DenseFact_fiuf: out of memory\n");
      return ret_out_of_memory;
    }

    for (int j = 0; j < n; ++j) {
      int N = j;
      int M = n - j - 1;

      // create temporary copy of row j, multiplied by pivots
      for (int i = 0; i < j; ++i) {
        temp[i] = A[j + i * lda] * A[i + i * lda];
      }

      // update diagonal element
      double Ajj = A[j + lda * j] - ddot(&N, &A[j], &lda, temp, &i_one);
      if (Ajj == 0.0 || isnan(Ajj)) {
        A[j + lda * j] = Ajj;
        printf("DenseFact_fiuf: invalid pivot\n");
        return ret_invalid_pivot;
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
    }
    // free temporary copy of row
    free(temp);
  } else {
    // allocate space for copy of col multiplied by pivots
    double* temp = malloc((n - 1) * sizeof(double));
    if (!temp) {
      printf("DenseFact_fiuf: out of memory\n");
      return ret_out_of_memory;
    }

    for (int j = 0; j < n; ++j) {
      int N = j;
      int M = n - j - 1;

      // create temporary copy of col j, multiplied by pivots
      for (int i = 0; i < j; ++i) {
        temp[i] = A[i + j * lda] * A[i + i * lda];
      }

      // update diagonal element
      double Ajj = A[j + lda * j] - ddot(&N, &A[j * lda], &i_one, temp, &i_one);
      if (Ajj == 0.0 || isnan(Ajj)) {
        A[j + lda * j] = Ajj;
        printf("DenseFact_fiuf: invalid pivot\n");
        return ret_invalid_pivot;
      }

      // save diagonal element
      A[j + lda * j] = Ajj;
      double coeff = 1.0 / Ajj;

      // compute column j
      if (j < n - 1) {
        dgemv(&TT, &N, &M, &d_m_one, &A[(j + 1) * lda], &lda, temp, &i_one,
              &d_one, &A[j + (j + 1) * lda], &lda);
        dscal(&M, &coeff, &A[j + (j + 1) * lda], &lda);
      }
    }

    // free temporary copy of row
    free(temp);
  }

  return ret_ok;
}

// ===========================================================================
// Functions to compute dense partial Cholesky or LDL factorizations with
// left-looking approach, with blocking.
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
// These functions are similar to Lapack dpotrf("L", n, A, lda, info).
//
// Arguments:
// - n      : Dimension of matrix M.
// - k      : Number of columns to factorize.
//            If k < n, a partial factorization is computed.
//            If k >= n, a full factorization is computed and B is not used.
// - nb     : size of the blocks.
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

int DenseFact_pdbf(int n, int k, int nb, double* restrict A, int lda,
                   double* restrict B, int ldb, double* times) {
  // ===========================================================================
  // Positive definite factorization with blocks.
  // BLAS calls: dsyrk, dgemm, dtrsm.
  // ===========================================================================

  // check input
  if (n < 0 || k < 0 || !A || lda < n || (k < n && (!B || ldb < n - k))) {
    printf("DenseFact_pdbf: invalid input\n");
    return ret_invalid_input;
  }

  if (!times) {
    printf("DenseFact_pdbf: invalid input\n");
    return ret_invalid_input;
  }

  // quick return
  if (n == 0) return ret_ok;

  // j is the starting col of the block of columns
  for (int j = 0; j < k; j += nb) {
    // jb is the size of the block
    int jb = min(nb, k - j);

    // sizes for blas calls
    int N = jb;
    int K = j;
    int M = n - j - jb;

    // starting position of matrices for BLAS calls
    double* D = &A[j + lda * j];
    double* P = &A[j];
    double* Q = &A[j + N];
    double* R = &A[j + N + lda * j];

// update diagonal block
#ifdef TIMING
    t0 = GetTime();
#endif

    dsyrk(&LL, &NN, &N, &K, &d_m_one, P, &lda, &d_one, D, &lda);

#ifdef TIMING
    times[t_dsyrk] += GetTime() - t0;
#endif

// factorize diagonal block
#ifdef TIMING
    t0 = GetTime();
#endif

    int info = DenseFact_fduf('L', N, D, lda);

#ifdef TIMING
    times[t_fact] += GetTime() - t0;
#endif
    if (info != 0) return info;

    if (j + jb < n) {
// update block of columns
#ifdef TIMING
      t0 = GetTime();
#endif
      dgemm(&NN, &TT, &M, &N, &K, &d_m_one, Q, &lda, P, &lda, &d_one, R, &lda);
#ifdef TIMING
      times[t_dgemm] += GetTime() - t0;
#endif

// solve block of columns with diagonal block
#ifdef TIMING
      t0 = GetTime();
#endif
      dtrsm(&RR, &LL, &TT, &NN, &M, &N, &d_one, D, &lda, R, &lda);
#ifdef TIMING
      times[t_dtrsm] += GetTime() - t0;
#endif
    }
  }

  // update Schur complement if partial factorization is required
  if (k < n) {
    int N = n - k;
#ifdef TIMING
    t0 = GetTime();
#endif
    dsyrk(&LL, &NN, &N, &k, &d_m_one, &A[k], &lda, &d_zero, B, &ldb);
#ifdef TIMING
    times[t_dsyrk] += GetTime() - t0;
#endif
  }

  return ret_ok;
}

int DenseFact_pibf(int n, int k, int nb, double* restrict A, int lda,
                   double* restrict B, int ldb, double* times) {
  // ===========================================================================
  // Indefinite factorization with blocks.
  // BLAS calls: dcopy, dscal, dgemm, dtrsm, dsyrk
  // ===========================================================================

  // check input
  if (n < 0 || k < 0 || !A || lda < n || (k < n && (!B || ldb < n - k))) {
    printf("DenseFact_pibf: invalid input\n");
    return ret_invalid_input;
  }

  // quick return
  if (n == 0) return ret_ok;

  // j is the starting col of the block of columns
  for (int j = 0; j < k; j += nb) {
    // jb is the size of the block
    int jb = min(nb, k - j);

    // sizes for blas calls
    int N = jb;
    int K = j;
    int M = n - j - jb;

    // starting position of matrices for BLAS calls
    double* D = &A[j + lda * j];
    double* P = &A[j];
    double* Q = &A[j + N];
    double* R = &A[j + N + lda * j];

    // create temporary copy of block of rows, multiplied by pivots
    double* T = malloc(j * jb * sizeof(double));
    if (!T) {
      printf("DenseFact_pibf: out of memory\n");
      return ret_out_of_memory;
    }
    int ldt = jb;
    for (int i = 0; i < j; ++i) {
      dcopy(&N, &A[j + i * lda], &i_one, &T[i * ldt], &i_one);
      dscal(&N, &A[i + i * lda], &T[i * ldt], &i_one);
    }

// update diagonal block using dgemm
#ifdef TIMING
    t0 = GetTime();
#endif
    dgemm(&NN, &TT, &jb, &jb, &j, &d_m_one, P, &lda, T, &ldt, &d_one, D, &lda);
#ifdef TIMING
    times[t_dgemm] += GetTime() - t0;
#endif

// factorize diagonal block
#ifdef TIMING
    t0 = GetTime();
#endif
    int info = DenseFact_fiuf('L', N, D, lda);
#ifdef TIMING
    times[t_fact] += GetTime() - t0;
#endif
    if (info != 0) return info;

    if (j + jb < n) {
// update block of columns
#ifdef TIMING
      t0 = GetTime();
#endif
      dgemm(&NN, &TT, &M, &N, &K, &d_m_one, Q, &lda, T, &ldt, &d_one, R, &lda);
#ifdef TIMING
      times[t_dgemm] += GetTime() - t0;
#endif

// solve block of columns with L
#ifdef TIMING
      t0 = GetTime();
#endif
      dtrsm(&RR, &LL, &TT, &UU, &M, &N, &d_one, D, &lda, R, &lda);
#ifdef TIMING
      times[t_dtrsm] += GetTime() - t0;
#endif

      // solve block of columns with D
      for (int i = 0; i < jb; ++i) {
        double coeff = 1.0 / A[j + i + (j + i) * lda];
        dscal(&M, &coeff, &A[j + jb + lda * (j + i)], &i_one);
      }
    }

    free(T);
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
    if (!temp_pos) {
      printf("DenseFact_pibf: out of memory\n");
      return ret_out_of_memory;
    }
    double* temp_neg = malloc((n - k) * neg_pivot * sizeof(double));
    if (!temp_neg) {
      printf("DenseFact_pibf: out of memory\n");
      return ret_out_of_memory;
    }

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
// In this way, I can use dsyrk instead of dgemm and avoid updating the
// full square schur complement. First call uses beta = 0.0, to clear
// content of B. Second call uses beta = 1.0, to not clear the result of
// the first call.
#ifdef TIMING
    t0 = GetTime();
#endif
    dsyrk(&LL, &NN, &N, &pos_pivot, &d_m_one, temp_pos, &ldt, &d_zero, B, &ldb);
    dsyrk(&LL, &NN, &N, &neg_pivot, &d_one, temp_neg, &ldt, &d_one, B, &ldb);
#ifdef TIMING
    times[t_dsyrk] += GetTime() - t0;
#endif

    free(temp_pos);
    free(temp_neg);
  }

  return ret_ok;
}

int DenseFact_pdbh(int n, int k, int nb, double* restrict A, double* restrict B,
                   double* times) {
  // ===========================================================================
  // Positive definite factorization with blocks in lower-blocked-hybrid
  // format. A should be in lower-blocked-hybrid format. Schur complement is
  // returned in B in lower packed format (not lower-blocked-hybrid).
  // BLAS calls: dsyrk, dgemm, dtrsm, dcopy
  // ===========================================================================

  // check input
  if (n < 0 || k < 0 || !A || (k < n && !B)) {
    printf("DenseFact_pdbh: invalid input\n");
    return ret_invalid_input;
  }

  // quick return
  if (n == 0) return ret_ok;

  // number of blocks of columns
  int n_blocks = (k - 1) / nb + 1;

  // start of diagonal blocks
  int* diag_start = malloc(n_blocks * sizeof(double));
  if (!diag_start) {
    printf("DenseFact_pdbh: out of memory\n");
    return ret_out_of_memory;
  }
  diag_start[0] = 0;
  for (int i = 1; i < n_blocks; ++i) {
    diag_start[i] =
        diag_start[i - 1] + nb * (2 * n - 2 * (i - 1) * nb - nb + 1) / 2;
  }

  // size of blocks
  int diag_size = nb * (nb + 1) / 2;
  int full_size = nb * nb;

  // variables for blas calls
  int info;

  // buffer for full-format diagonal blocks
  double* D = malloc(nb * nb * sizeof(double));
  if (!D) {
    printf("DenseFact_pdbh: out of memory\n");
    return ret_out_of_memory;
  }

  // j is the index of the block column
  for (int j = 0; j < n_blocks; ++j) {
    // jb is the number of columns
    int jb = min(nb, k - nb * j);

    // size of current block could be smaller than diag_size and full_size
    int this_diag_size = jb * (jb + 1) / 2;
    int this_full_size = nb * jb;

// full copy of diagonal block by rows, in D
#ifdef TIMING
    t0 = GetTime();
#endif
    int offset = 0;
    for (int Drow = 0; Drow < jb; ++Drow) {
      int N = Drow + 1;
      dcopy(&N, &A[diag_start[j] + offset], &i_one, &D[Drow * jb], &i_one);
      offset += N;
    }
#ifdef TIMING
    times[t_dcopy] += GetTime() - t0;
#endif

    // number of rows left below block j
    int M = n - nb * j - jb;

    // starting position of block of columns below diagonal block j
    int Lij_pos = diag_start[j] + this_diag_size;

    // update diagonal block and block of columns
    for (int k = 0; k < j; ++k) {
      // starting position of block to input to dsyrk
      int Ljk_pos = diag_start[k] + diag_size;
      if (j > k + 1) Ljk_pos += full_size * (j - k - 1);

#ifdef TIMING
      t0 = GetTime();
#endif
      dsyrk(&UU, &TT, &jb, &nb, &d_m_one, &A[Ljk_pos], &nb, &d_one, D, &jb);
#ifdef TIMING
      times[t_dsyrk] += GetTime() - t0;
#endif

      if (M > 0) {
        int Lik_pos = Ljk_pos + this_full_size;
#ifdef TIMING
        t0 = GetTime();
#endif
        dgemm(&TT, &NN, &jb, &M, &nb, &d_m_one, &A[Ljk_pos], &nb, &A[Lik_pos],
              &nb, &d_one, &A[Lij_pos], &jb);
#ifdef TIMING
        times[t_dgemm] += GetTime() - t0;
#endif
      }
    }

// factorize diagonal block
#ifdef TIMING
    t0 = GetTime();
#endif
    int info = DenseFact_fduf('U', jb, D, jb);
#ifdef TIMING
    times[t_fact] += GetTime() - t0;
#endif
    if (info != 0) return info;

    if (M > 0) {
// solve block of columns with diagonal block
#ifdef TIMING
      t0 = GetTime();
#endif
      dtrsm(&LL, &UU, &TT, &NN, &jb, &M, &d_one, D, &jb, &A[Lij_pos], &jb);
#ifdef TIMING
      times[t_dtrsm] += GetTime() - t0;
#endif
    }

// put D back into packed format
#ifdef TIMING
    t0 = GetTime();
#endif
    offset = 0;
    for (int Drow = 0; Drow < jb; ++Drow) {
      int N = Drow + 1;
      dcopy(&N, &D[Drow * jb], &i_one, &A[diag_start[j] + offset], &i_one);
      offset += N;
    }
#ifdef TIMING
    times[t_dcopy] += GetTime() - t0;
#endif
  }
  free(D);

  // compute Schur complement if partial factorization is required
  if (k < n) {
    // number of rows/columns in the Schur complement
    int ns = n - k;

    // size of last full block (may be smaller than full_size)
    int ncol_last = k % nb;
    int last_full_size = ncol_last == 0 ? full_size : ncol_last * nb;

    double beta = 0.0;

    // buffer for full-format of block of columns of Schur complement
    double* schur_buf = malloc(ns * nb * sizeof(double));
    if (!schur_buf) {
      printf("DenseFact_pdbh: out of memory\n");
      return ret_out_of_memory;
    }

    // number of blocks in Schur complement
    int s_blocks = (ns - 1) / nb + 1;

    // index to write into B
    int B_start = 0;

    // Go through block of columns of Schur complement.
    // Each block will have a full-format copy made in schur_buf
    for (int sb = 0; sb < s_blocks; ++sb) {
      // number of rows of the block
      int nrow = ns - nb * sb;

      // number of columns of the block
      int ncol = min(nb, nrow);

      // number of elements in diagonal block
      int schur_diag_size = ncol * (ncol + 1) / 2;

      beta = 0.0;

      // each block receives contributions from the blocks of the leading part
      // of A
      for (int j = 0; j < n_blocks; ++j) {
        int jb = min(nb, k - nb * j);
        int this_diag_size = jb * (jb + 1) / 2;
        int this_full_size = nb * jb;

        // compute index to access diagonal block in A
        int diag_pos = diag_start[j] + this_diag_size;
        if (j < n_blocks - 1) {
          diag_pos += (n_blocks - j - 2) * full_size + last_full_size;
        }
        diag_pos += sb * this_full_size;

// update diagonal block
#ifdef TIMING
        t0 = GetTime();
#endif
        dsyrk(&UU, &TT, &ncol, &jb, &d_m_one, &A[diag_pos], &jb, &beta,
              schur_buf, &ncol);
#ifdef TIMING
        times[t_dsyrk] += GetTime() - t0;
#endif

        // update subdiagonal part
        int M = nrow - nb;
        if (M > 0) {
#ifdef TIMING
          t0 = GetTime();
#endif
          dgemm(&TT, &NN, &nb, &M, &jb, &d_m_one, &A[diag_pos], &jb,
                &A[diag_pos + this_full_size], &jb, &beta,
                &schur_buf[ncol * ncol], &ncol);
#ifdef TIMING
          times[t_dgemm] += GetTime() - t0;
#endif
        }

        // beta is 0 for the first time (to avoid initializing schur_buf) and
        // 1 for the next calls
        beta = 1.0;
      }

// schur_buf contains Schur complement in hybrid format (with full
// diagonal blocks). Put it in lower-packed format in B.
#ifdef TIMING
      t0 = GetTime();
#endif
      for (int buf_row = 0; buf_row < nrow; ++buf_row) {
        int N = ncol;
        dcopy(&N, &schur_buf[buf_row * ncol], &i_one, &B[B_start + buf_row],
              &nrow);
      }
      B_start += nrow * ncol;

#ifdef TIMING
      times[t_dcopy_schur] += GetTime() - t0;
#endif
    }

    free(schur_buf);
  }

  free(diag_start);

  return ret_ok;
}

int DenseFact_pibh(int n, int k, int nb, double* restrict A, double* restrict B,
                   double* times) {
  // ===========================================================================
  // Indefinite factorization with blocks in lower-blocked-hybrid format.
  // A should be in lower-blocked-hybrid format. Schur complement is returned
  // in B in lower packed format (not lower-blocked-hybrid). BLAS calls:
  // dsyrk, dgemm, dtrsm, dcopy, dscal
  // ===========================================================================

  const int sizeA = n * k - k * (k - 1) / 2;

  // check input
  if (n < 0 || k < 0 || !A || (k < n && !B)) {
    printf("DenseFact_pibh: invalid input\n");
    return ret_invalid_input;
  }

  // quick return
  if (n == 0) return ret_ok;

  // number of blocks of columns
  int n_blocks = (k - 1) / nb + 1;

  // start of diagonal blocks
  int* diag_start = malloc(n_blocks * sizeof(double));
  if (!diag_start) {
    printf("DenseFact_pibh: out of memory\n");
    return ret_out_of_memory;
  }
  diag_start[0] = 0;
  for (int i = 1; i < n_blocks; ++i) {
    diag_start[i] =
        diag_start[i - 1] + nb * (2 * n - 2 * (i - 1) * nb - nb + 1) / 2;
  }

  // size of blocks
  int diag_size = nb * (nb + 1) / 2;
  int full_size = nb * nb;

  // variables for blas calls
  int info;

  // buffer for full-format diagonal blocks
  double* D = malloc(nb * nb * sizeof(double));
  if (!D) {
    printf("DenseFact_pibh: out of memory\n");
    return ret_out_of_memory;
  }

  // buffer for copy of block scaled by pivots
  double* T = malloc(nb * nb * sizeof(double));
  if (!T) {
    printf("DenseFact_pibh: out of memory\n");
    return ret_out_of_memory;
  }

  // j is the index of the block column
  for (int j = 0; j < n_blocks; ++j) {
    // jb is the number of columns
    int jb = min(nb, k - nb * j);

    // size of current block could be smaller than diag_size and full_size
    int this_diag_size = jb * (jb + 1) / 2;
    int this_full_size = nb * jb;

// full copy of diagonal block by rows, in D
#ifdef TIMING
    t0 = GetTime();
#endif
    int offset = 0;
    for (int Drow = 0; Drow < jb; ++Drow) {
      int N = Drow + 1;
      dcopy(&N, &A[diag_start[j] + offset], &i_one, &D[Drow * jb], &i_one);
      offset += N;
    }
#ifdef TIMING
    times[t_dcopy] += GetTime() - t0;
#endif

    // number of rows left below block j
    int M = n - nb * j - jb;

    // starting position of block of columns below diagonal block j
    int Lij_pos = diag_start[j] + this_diag_size;

    // update diagonal block and block of columns
    for (int k = 0; k < j; ++k) {
      // starting position of block to input to dsyrk
      int Ljk_pos = diag_start[k] + diag_size;
      if (j > k + 1) Ljk_pos += full_size * (j - k - 1);

// copy block jk into temp
#ifdef TIMING
      t0 = GetTime();
#endif
      dcopy(&this_full_size, &A[Ljk_pos], &i_one, T, &i_one);
#ifdef TIMING
      times[t_dcopy] += GetTime() - t0;
#endif

// scale temp by pivots
#ifdef TIMING
      t0 = GetTime();
#endif
      int pivot_pos = diag_start[k];
      for (int col = 0; col < nb; ++col) {
        dscal(&jb, &A[pivot_pos], &T[col], &nb);
        pivot_pos += col + 2;
      }
#ifdef TIMING
      times[t_dscal] += GetTime() - t0;
#endif

// update diagonal block with dgemm
#ifdef TIMING
      t0 = GetTime();
#endif
      dgemm(&TT, &NN, &jb, &jb, &nb, &d_m_one, T, &nb, &A[Ljk_pos], &nb, &d_one,
            D, &jb);
#ifdef TIMING
      times[t_dgemm] += GetTime() - t0;
#endif

      // update rectangular block
      if (M > 0) {
        int Lik_pos = Ljk_pos + this_full_size;
#ifdef TIMING
        t0 = GetTime();
#endif
        dgemm(&TT, &NN, &jb, &M, &nb, &d_m_one, T, &nb, &A[Lik_pos], &nb,
              &d_one, &A[Lij_pos], &jb);
#ifdef TIMING
        times[t_dgemm] += GetTime() - t0;
#endif
      }
    }

// factorize diagonal block
#ifdef TIMING
    t0 = GetTime();
#endif
    int info = DenseFact_fiuf('U', jb, D, jb);
#ifdef TIMING
    times[t_fact] += GetTime() - t0;
#endif
    if (info != 0) return info;

    if (M > 0) {
// solve block of columns with diagonal block
#ifdef TIMING
      t0 = GetTime();
#endif
      dtrsm(&LL, &UU, &TT, &UU, &jb, &M, &d_one, D, &jb, &A[Lij_pos], &jb);
#ifdef TIMING
      times[t_dtrsm] += GetTime() - t0;
#endif

// scale columns by pivots
#ifdef TIMING
      t0 = GetTime();
#endif
      int pivot_pos = diag_start[j];
      for (int col = 0; col < jb; ++col) {
        double coeff = 1.0 / A[pivot_pos];
        dscal(&M, &coeff, &A[Lij_pos + col], &jb);
        pivot_pos += col + 2;
      }
#ifdef TIMING
      times[t_dscal] += GetTime() - t0;
#endif
    }

// put D back into packed format
#ifdef TIMING
    t0 = GetTime();
#endif
    offset = 0;
    for (int Drow = 0; Drow < jb; ++Drow) {
      int N = Drow + 1;
      dcopy(&N, &D[Drow * jb], &i_one, &A[diag_start[j] + offset], &i_one);
      offset += N;
    }
#ifdef TIMING
    times[t_dcopy] += GetTime() - t0;
#endif
  }
  free(D);

  // compute Schur complement if partial factorization is required
  if (k < n) {
    // number of rows/columns in the Schur complement
    int ns = n - k;

    // size of last full block (may be smaller than full_size)
    int ncol_last = k % nb;
    int last_full_size = ncol_last == 0 ? full_size : ncol_last * nb;

    double beta = 0.0;

    // buffer for full-format of block of columns of Schur complement
    double* schur_buf = malloc(ns * nb * sizeof(double));
    if (!schur_buf) {
      printf("DenseFact_pibh: out of memory\n");
      return ret_out_of_memory;
    }

    // number of blocks in Schur complement
    int s_blocks = (ns - 1) / nb + 1;

    // index to write into B
    int B_start = 0;

    // Go through block of columns of Schur complement.
    // Each block will have a full-format copy made in schur_buf
    for (int sb = 0; sb < s_blocks; ++sb) {
      // number of rows of the block
      int nrow = ns - nb * sb;

      // number of columns of the block
      int ncol = min(nb, nrow);

      // number of elements in diagonal block
      int schur_diag_size = ncol * (ncol + 1) / 2;

      beta = 0.0;

      // each block receives contributions from the blocks of the leading part
      // of A
      for (int j = 0; j < n_blocks; ++j) {
        int jb = min(nb, k - nb * j);
        int this_diag_size = jb * (jb + 1) / 2;
        int this_full_size = nb * jb;

        // compute index to access diagonal block in A
        int diag_pos = diag_start[j] + this_diag_size;
        if (j < n_blocks - 1) {
          diag_pos += (n_blocks - j - 2) * full_size + last_full_size;
        }
        diag_pos += sb * this_full_size;

        // create copy of block, multiplied by pivots
        int N = ncol * jb;
#ifdef TIMING
        t0 = GetTime();
#endif
        dcopy(&N, &A[diag_pos], &i_one, T, &i_one);
#ifdef TIMING
        times[t_dcopy] += GetTime() - t0;
#endif
        int pivot_pos = diag_start[j];
#ifdef TIMING
        t0 = GetTime();
#endif
        for (int col = 0; col < jb; ++col) {
          dscal(&ncol, &A[pivot_pos], &T[col], &jb);
          pivot_pos += col + 2;
        }
#ifdef TIMING
        times[t_dscal] += GetTime() - t0;
#endif

// update diagonal block using dgemm
#ifdef TIMING
        t0 = GetTime();
#endif

        // printf("%p\n", A);
        assert(diag_pos + ncol * jb - 1 <= sizeA - 1);

        dgemm(&TT, &NN, &ncol, &ncol, &jb, &d_m_one, T, &jb, &A[diag_pos], &jb,
              &beta, schur_buf, &ncol);
#ifdef TIMING
        times[t_dgemm] += GetTime() - t0;
#endif

        // update subdiagonal part
        int M = nrow - nb;
        if (M > 0) {
#ifdef TIMING
          t0 = GetTime();
#endif
          dgemm(&TT, &NN, &ncol, &M, &jb, &d_m_one, T, &jb,
                &A[diag_pos + this_full_size], &jb, &beta,
                &schur_buf[ncol * ncol], &ncol);
#ifdef TIMING
          times[t_dgemm] += GetTime() - t0;
#endif
        }

        // beta is 0 for the first time (to avoid initializing schur_buf) and
        // 1 for the next calls
        beta = 1.0;
      }

// schur_buf contains Schur complement in hybrid format (with full
// diagonal blocks). Put it in lower-packed format in B.
#ifdef TIMING
      t0 = GetTime();
#endif
      for (int buf_row = 0; buf_row < nrow; ++buf_row) {
        int N = ncol;
        dcopy(&N, &schur_buf[buf_row * ncol], &i_one, &B[B_start + buf_row],
              &nrow);
      }
      B_start += nrow * ncol;
#ifdef TIMING
      times[t_dcopy_schur] += GetTime() - t0;
#endif
    }

    free(schur_buf);
  }

  free(T);
  free(diag_start);

  return ret_ok;
}

int DenseFact_pdbh_2(int n, int k, int nb, double* A, double* B,
                     double* times) {
  // ===========================================================================
  // Positive definite factorization with blocks in lower-blocked-hybrid
  // format. A should be in lower-blocked-hybrid format. Schur complement is
  // returned in B in lower-block-hybrid format, with full diagonal blocks.
  // BLAS calls: dsyrk, dgemm, dtrsm, dcopy
  // ===========================================================================

  // check input
  if (n < 0 || k < 0 || !A || (k < n && !B)) {
    printf("DenseFact_pdbh: invalid input\n");
    return ret_invalid_input;
  }

  // quick return
  if (n == 0) return ret_ok;

  // number of blocks of columns
  int n_blocks = (k - 1) / nb + 1;

  // start of diagonal blocks
  int* diag_start = malloc(n_blocks * sizeof(double));
  if (!diag_start) {
    printf("DenseFact_pdbh: out of memory\n");
    return ret_out_of_memory;
  }
  diag_start[0] = 0;
  for (int i = 1; i < n_blocks; ++i) {
    diag_start[i] =
        diag_start[i - 1] + nb * (2 * n - 2 * (i - 1) * nb - nb + 1) / 2;
  }

  // size of blocks
  int diag_size = nb * (nb + 1) / 2;
  int full_size = nb * nb;

  // variables for blas calls
  int info;

  // buffer for full-format diagonal blocks
  double* D = malloc(nb * nb * sizeof(double));
  if (!D) {
    printf("DenseFact_pdbh: out of memory\n");
    return ret_out_of_memory;
  }

  // j is the index of the block column
  for (int j = 0; j < n_blocks; ++j) {
    // jb is the number of columns
    int jb = min(nb, k - nb * j);

    // size of current block could be smaller than diag_size and full_size
    int this_diag_size = jb * (jb + 1) / 2;
    int this_full_size = nb * jb;

// full copy of diagonal block by rows, in D
#ifdef TIMING
    t0 = GetTime();
#endif
    int offset = 0;
    for (int Drow = 0; Drow < jb; ++Drow) {
      int N = Drow + 1;
      dcopy(&N, &A[diag_start[j] + offset], &i_one, &D[Drow * jb], &i_one);
      offset += N;
    }
#ifdef TIMING
    times[t_dcopy] += GetTime() - t0;
#endif

    // number of rows left below block j
    int M = n - nb * j - jb;

    // starting position of block of columns below diagonal block j
    int Lij_pos = diag_start[j] + this_diag_size;

    // update diagonal block and block of columns
    for (int k = 0; k < j; ++k) {
      // starting position of block to input to dsyrk
      int Ljk_pos = diag_start[k] + diag_size;
      if (j > k + 1) Ljk_pos += full_size * (j - k - 1);

#ifdef TIMING
      t0 = GetTime();
#endif
      dsyrk(&UU, &TT, &jb, &nb, &d_m_one, &A[Ljk_pos], &nb, &d_one, D, &jb);
#ifdef TIMING
      times[t_dsyrk] += GetTime() - t0;
#endif

      if (M > 0) {
        int Lik_pos = Ljk_pos + this_full_size;
#ifdef TIMING
        t0 = GetTime();
#endif
        dgemm(&TT, &NN, &jb, &M, &nb, &d_m_one, &A[Ljk_pos], &nb, &A[Lik_pos],
              &nb, &d_one, &A[Lij_pos], &jb);
#ifdef TIMING
        times[t_dgemm] += GetTime() - t0;
#endif
      }
    }

// factorize diagonal block
#ifdef TIMING
    t0 = GetTime();
#endif
    int info = DenseFact_fduf('U', jb, D, jb);
#ifdef TIMING
    times[t_fact] += GetTime() - t0;
#endif
    if (info != 0) return info;

    if (M > 0) {
// solve block of columns with diagonal block
#ifdef TIMING
      t0 = GetTime();
#endif
      dtrsm(&LL, &UU, &TT, &NN, &jb, &M, &d_one, D, &jb, &A[Lij_pos], &jb);
#ifdef TIMING
      times[t_dtrsm] += GetTime() - t0;
#endif
    }

// put D back into packed format
#ifdef TIMING
    t0 = GetTime();
#endif
    offset = 0;
    for (int Drow = 0; Drow < jb; ++Drow) {
      int N = Drow + 1;
      dcopy(&N, &D[Drow * jb], &i_one, &A[diag_start[j] + offset], &i_one);
      offset += N;
    }
#ifdef TIMING
    times[t_dcopy] += GetTime() - t0;
#endif
  }
  free(D);

  // compute Schur complement if partial factorization is required
  if (k < n) {
    // number of rows/columns in the Schur complement
    int ns = n - k;

    // size of last full block (may be smaller than full_size)
    int ncol_last = k % nb;
    int last_full_size = ncol_last == 0 ? full_size : ncol_last * nb;

    double beta = 0.0;

    // number of blocks in Schur complement
    int s_blocks = (ns - 1) / nb + 1;

    // index to write into B
    int B_start = 0;

    // Go through block of columns of Schur complement.
    // Each block will have a full-format copy made in schur_buf
    for (int sb = 0; sb < s_blocks; ++sb) {
      // location in B where to write the next block of columns
      double* schur_buf = &B[B_start];

      // number of rows of the block
      int nrow = ns - nb * sb;

      // number of columns of the block
      int ncol = min(nb, nrow);

      beta = 0.0;

      // each block receives contributions from the blocks of the leading part
      // of A
      for (int j = 0; j < n_blocks; ++j) {
        int jb = min(nb, k - nb * j);
        int this_diag_size = jb * (jb + 1) / 2;
        int this_full_size = nb * jb;

        // compute index to access diagonal block in A
        int diag_pos = diag_start[j] + this_diag_size;
        if (j < n_blocks - 1) {
          diag_pos += (n_blocks - j - 2) * full_size + last_full_size;
        }
        diag_pos += sb * this_full_size;

// update diagonal block
#ifdef TIMING
        t0 = GetTime();
#endif
        dsyrk(&UU, &TT, &ncol, &jb, &d_m_one, &A[diag_pos], &jb, &beta,
              schur_buf, &ncol);
#ifdef TIMING
        times[t_dsyrk] += GetTime() - t0;
#endif

        // update subdiagonal part
        int M = nrow - nb;
        if (M > 0) {
#ifdef TIMING
          t0 = GetTime();
#endif
          dgemm(&TT, &NN, &nb, &M, &jb, &d_m_one, &A[diag_pos], &jb,
                &A[diag_pos + this_full_size], &jb, &beta,
                &schur_buf[ncol * ncol], &ncol);
#ifdef TIMING
          times[t_dgemm] += GetTime() - t0;
#endif
        }

        // beta is 0 for the first time (to avoid initializing schur_buf) and
        // 1 for the next calls
        beta = 1.0;
      }

      B_start += nrow * ncol;
    }
  }

  free(diag_start);

  return ret_ok;
}

int DenseFact_pibh_2(int n, int k, int nb, double* A, double* B,
                     double* times) {
  // ===========================================================================
  // Indefinite factorization with blocks in lower-blocked-hybrid format.
  // A should be in lower-blocked-hybrid format. Schur complement is returned
  // in B in in lower-block-hybrid format, with full diagonal blocks.
  // BLAS calls: dsyrk, dgemm, dtrsm, dcopy, dscal
  // ===========================================================================

  // check input
  if (n < 0 || k < 0 || !A || (k < n && !B)) {
    printf("DenseFact_pibh: invalid input\n");
    return ret_invalid_input;
  }

  // quick return
  if (n == 0) return ret_ok;

  // number of blocks of columns
  int n_blocks = (k - 1) / nb + 1;

  // start of diagonal blocks
  int* diag_start = malloc(n_blocks * sizeof(double));
  if (!diag_start) {
    printf("DenseFact_pibh: out of memory\n");
    return ret_out_of_memory;
  }
  diag_start[0] = 0;
  for (int i = 1; i < n_blocks; ++i) {
    diag_start[i] =
        diag_start[i - 1] + nb * (2 * n - 2 * (i - 1) * nb - nb + 1) / 2;
  }

  // size of blocks
  int diag_size = nb * (nb + 1) / 2;
  int full_size = nb * nb;

  // variables for blas calls
  int info;

  // buffer for full-format diagonal blocks
  double* D = malloc(nb * nb * sizeof(double));
  if (!D) {
    printf("DenseFact_pibh: out of memory\n");
    return ret_out_of_memory;
  }

  // buffer for copy of block scaled by pivots
  double* T = malloc(nb * nb * sizeof(double));
  if (!T) {
    printf("DenseFact_pibh: out of memory\n");
    return ret_out_of_memory;
  }

  // j is the index of the block column
  for (int j = 0; j < n_blocks; ++j) {
    // jb is the number of columns
    int jb = min(nb, k - nb * j);

    // size of current block could be smaller than diag_size and full_size
    int this_diag_size = jb * (jb + 1) / 2;
    int this_full_size = nb * jb;

// full copy of diagonal block by rows, in D
#ifdef TIMING
    t0 = GetTime();
#endif
    int offset = 0;
    for (int Drow = 0; Drow < jb; ++Drow) {
      int N = Drow + 1;
      dcopy(&N, &A[diag_start[j] + offset], &i_one, &D[Drow * jb], &i_one);
      offset += N;
    }
#ifdef TIMING
    times[t_dcopy] += GetTime() - t0;
#endif

    // number of rows left below block j
    int M = n - nb * j - jb;

    // starting position of block of columns below diagonal block j
    int Lij_pos = diag_start[j] + this_diag_size;

    // update diagonal block and block of columns
    for (int k = 0; k < j; ++k) {
      // starting position of block to input to dsyrk
      int Ljk_pos = diag_start[k] + diag_size;
      if (j > k + 1) Ljk_pos += full_size * (j - k - 1);

// copy block jk into temp
#ifdef TIMING
      t0 = GetTime();
#endif
      dcopy(&this_full_size, &A[Ljk_pos], &i_one, T, &i_one);
#ifdef TIMING
      times[t_dcopy] += GetTime() - t0;
#endif

// scale temp by pivots
#ifdef TIMING
      t0 = GetTime();
#endif
      int pivot_pos = diag_start[k];
      for (int col = 0; col < nb; ++col) {
        dscal(&jb, &A[pivot_pos], &T[col], &nb);
        pivot_pos += col + 2;
      }
#ifdef TIMING
      times[t_dscal] += GetTime() - t0;
#endif

// update diagonal block with dgemm
#ifdef TIMING
      t0 = GetTime();
#endif
      dgemm(&TT, &NN, &jb, &jb, &nb, &d_m_one, T, &nb, &A[Ljk_pos], &nb, &d_one,
            D, &jb);
#ifdef TIMING
      times[t_dgemm] += GetTime() - t0;
#endif

      // update rectangular block
      if (M > 0) {
        int Lik_pos = Ljk_pos + this_full_size;
#ifdef TIMING
        t0 = GetTime();
#endif
        dgemm(&TT, &NN, &jb, &M, &nb, &d_m_one, T, &nb, &A[Lik_pos], &nb,
              &d_one, &A[Lij_pos], &jb);
#ifdef TIMING
        times[t_dgemm] += GetTime() - t0;
#endif
      }
    }

// factorize diagonal block
#ifdef TIMING
    t0 = GetTime();
#endif
    int info = DenseFact_fiuf('U', jb, D, jb);
#ifdef TIMING
    times[t_fact] += GetTime() - t0;
#endif
    if (info != 0) return info;

    if (M > 0) {
// solve block of columns with diagonal block
#ifdef TIMING
      t0 = GetTime();
#endif
      dtrsm(&LL, &UU, &TT, &UU, &jb, &M, &d_one, D, &jb, &A[Lij_pos], &jb);
#ifdef TIMING
      times[t_dtrsm] += GetTime() - t0;
#endif

// scale columns by pivots
#ifdef TIMING
      t0 = GetTime();
#endif
      int pivot_pos = diag_start[j];
      for (int col = 0; col < jb; ++col) {
        double coeff = 1.0 / A[pivot_pos];
        dscal(&M, &coeff, &A[Lij_pos + col], &jb);
        pivot_pos += col + 2;
      }
#ifdef TIMING
      times[t_dscal] += GetTime() - t0;
#endif
    }

    // put D back into packed format
#ifdef TIMING
    t0 = GetTime();
#endif
    offset = 0;
    for (int Drow = 0; Drow < jb; ++Drow) {
      int N = Drow + 1;
      dcopy(&N, &D[Drow * jb], &i_one, &A[diag_start[j] + offset], &i_one);
      offset += N;
    }
#ifdef TIMING
    times[t_dcopy] += GetTime() - t0;
#endif
  }
  free(D);

  // compute Schur complement if partial factorization is required
  if (k < n) {
    // number of rows/columns in the Schur complement
    int ns = n - k;

    // size of last full block (may be smaller than full_size)
    int ncol_last = k % nb;
    int last_full_size = ncol_last == 0 ? full_size : ncol_last * nb;

    double beta = 0.0;

    // number of blocks in Schur complement
    int s_blocks = (ns - 1) / nb + 1;

    // index to write into B
    int B_start = 0;

    // Go through block of columns of Schur complement.
    // Each block will have a full-format copy made in schur_buf
    for (int sb = 0; sb < s_blocks; ++sb) {
      // location in B where to write the next block of columns
      double* schur_buf = &B[B_start];

      // number of rows of the block
      int nrow = ns - nb * sb;

      // number of columns of the block
      int ncol = min(nb, nrow);

      beta = 0.0;

      // each block receives contributions from the blocks of the leading part
      // of A
      for (int j = 0; j < n_blocks; ++j) {
        int jb = min(nb, k - nb * j);
        int this_diag_size = jb * (jb + 1) / 2;
        int this_full_size = nb * jb;

        // compute index to access diagonal block in A
        int diag_pos = diag_start[j] + this_diag_size;
        if (j < n_blocks - 1) {
          diag_pos += (n_blocks - j - 2) * full_size + last_full_size;
        }
        diag_pos += sb * this_full_size;

        // create copy of block, multiplied by pivots
        int N = ncol * jb;
#ifdef TIMING
        t0 = GetTime();
#endif
        dcopy(&N, &A[diag_pos], &i_one, T, &i_one);
#ifdef TIMING
        times[t_dcopy] += GetTime() - t0;
#endif
        int pivot_pos = diag_start[j];
#ifdef TIMING
        t0 = GetTime();
#endif
        for (int col = 0; col < jb; ++col) {
          dscal(&ncol, &A[pivot_pos], &T[col], &jb);
          pivot_pos += col + 2;
        }
#ifdef TIMING
        times[t_dscal] += GetTime() - t0;
#endif

// update diagonal block using dgemm
#ifdef TIMING
        t0 = GetTime();
#endif
        dgemm(&TT, &NN, &ncol, &ncol, &jb, &d_m_one, T, &jb, &A[diag_pos], &jb,
              &beta, schur_buf, &ncol);
#ifdef TIMING
        times[t_dgemm] += GetTime() - t0;
#endif

        // update subdiagonal part
        int M = nrow - nb;
        if (M > 0) {
#ifdef TIMING
          t0 = GetTime();
#endif
          dgemm(&TT, &NN, &ncol, &M, &jb, &d_m_one, T, &jb,
                &A[diag_pos + this_full_size], &jb, &beta,
                &schur_buf[ncol * ncol], &ncol);
#ifdef TIMING
          times[t_dgemm] += GetTime() - t0;
#endif
        }

        // beta is 0 for the first time (to avoid initializing schur_buf) and
        // 1 for the next calls
        beta = 1.0;
      }
      B_start += nrow * ncol;
    }
  }

  free(T);
  free(diag_start);

  return ret_ok;
}

int DenseFact_l2h(double* restrict A, int nrow, int ncol, int nb,
                  double* times) {
  // ===========================================================================
  // Takes a matrix in lower-packed format, with nrow rows.
  // Converts the first ncol columns into lower-blocked-hybrid format, with
  // block size nb.
  // BLAS calls: dcopy
  // ===========================================================================

#ifdef TIMING
  t0 = GetTime();
#endif

  double* buf = malloc(nrow * nb * sizeof(double));
  if (!buf) {
    printf("DenseFact_l2h: out of memory\n");
    return ret_out_of_memory;
  }
  int startAtoBuf = 0;
  int startBuftoA = 0;

  for (int k = 0; k <= (ncol - 1) / nb; ++k) {
    // Each block of columns is copied into the buffer, leaving space in
    // between columns, to align rows. This "empty space" contains elements
    // from the previous use of buffer, but they are ignored when copying
    // back.

    // Number of columns in the block. Can be smaller than nb for last block.
    int block_size = min(nb, ncol - k * nb);

    // Number of rows in the block
    int row_size = nrow - k * nb;

    int startBuf = 0;
    for (int j = 0; j < block_size; ++j) {
      // Copy each column in the block
      int N = row_size - j;
      int i_one = 1;
      dcopy(&N, &A[startAtoBuf], &i_one, &buf[startBuf], &i_one);
      startAtoBuf += N;

      // +1 to align rows
      startBuf += row_size + 1;
    }

    // Copy columns back into A, row by row.
    // One call of dcopy for each row of the block of columns.
    for (int i = 0; i < row_size; ++i) {
      int N = min(i + 1, block_size);
      int i_one = 1;
      dcopy(&N, &buf[i], &row_size, &A[startBuftoA], &i_one);
      startBuftoA += N;
    }
  }

  free(buf);

#ifdef TIMING
  times[t_convert] += GetTime() - t0;
#endif

  return ret_ok;
}
