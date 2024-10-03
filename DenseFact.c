#include "DenseFact.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "CallAndTimeBlas.h"
#include "DenseFact_declaration.h"
#include "ReturnValues.h"
#include "timing.h"

/*
Names:
dense_fact_(pf)(di)(bu)(flhs)

pf: Partial or Full factorization
di: (positive) Definite or Indefinite
bu: Blocked or Unblocked
flhs: Full format, Lower packed format, lower-blocked-Hybrid packed format
      with packed schur complement, lower-blocked-hybrid packed format with
      hybrid Schur complement

*/

int dense_fact_fduf(char uplo, int n, double* restrict A, int lda,
                    double thresh, double* regul) {
  // ===========================================================================
  // Positive definite factorization without blocks.
  // BLAS calls: ddot_, dgemv_, dscal_.
  // ===========================================================================

  // check input
  if (n < 0 || !A || lda < n || (uplo != 'L' && uplo != 'U')) {
    printf("\ndense_fact_fduf: invalid input\n");
    return kRetInvalidInput;
  }

  // quick return
  if (n == 0) return kRetOk;

  // main operations
  if (uplo == 'L') {
    for (int j = 0; j < n; ++j) {
      const int N = j;
      const int M = n - j - 1;

      // update diagonal element
      double Ajj = A[j + lda * j] - ddot_(&N, &A[j], &lda, &A[j], &lda);
      if (isnan(Ajj)) {
        A[j + lda * j] = Ajj;
        printf("\ndense_fact_fduf: invalid pivot %e\n", Ajj);
        return kRetInvalidPivot;
      }

      // update column j
      if (j < n - 1) {
        dgemv_(&c_N, &M, &N, &d_m_one, &A[j + 1], &lda, &A[j], &lda, &d_one,
               &A[j + 1 + j * lda], &i_one);
      }

      // compute diagonal element
      if (Ajj <= thresh) {
        // if pivot is not acceptable, push it up to thresh
        // printf("small pivot %e ", Ajj);
        double old_pivot = Ajj;
        Ajj = thresh;

        // compute the minimum pivot required to keep the diagonal of the
        // current block acceptable:
        // b is column below pivot, d is diagonal of block, p is pivot
        // we want d_k - b_k^2 / p \ge thresh
        // i.e.
        // p \ge b_k^2 / (d_k - thresh)
        //
        double required_pivot = 0.0;
        for (int k = j + 1; k < n; ++k) {
          double bk = A[k + j * lda];
          double dk = A[k + lda * k];
          double temp = (dk - thresh);
          temp = (bk * bk) / temp;
          required_pivot = max(required_pivot, temp);
        }

        Ajj = max(Ajj, required_pivot);

        // record regularization used
        regul[j] += Ajj - old_pivot;

        // printf("set to %e\n", Ajj);
      }

      Ajj = sqrt(Ajj);
      A[j + lda * j] = Ajj;
      const double coeff = 1.0 / Ajj;

      // scale column j
      if (j < n - 1) {
        dscal_(&M, &coeff, &A[j + 1 + j * lda], &i_one);
      }
    }
  } else {
    for (int j = 0; j < n; ++j) {
      const int N = j;
      const int M = n - j - 1;

      // update diagonal element
      double Ajj =
          A[j + lda * j] - ddot_(&N, &A[lda * j], &i_one, &A[lda * j], &i_one);
      if (isnan(Ajj)) {
        A[j + lda * j] = Ajj;
        printf("\ndense_fact_fduf: invalid pivot %e\n", Ajj);
        return kRetInvalidPivot;
      }

      // update column j
      if (j < n - 1) {
        dgemv_(&c_T, &N, &M, &d_m_one, &A[lda * (j + 1)], &lda, &A[lda * j],
               &i_one, &d_one, &A[j + (j + 1) * lda], &lda);
      }

      // compute diagonal element
      if (Ajj <= thresh) {
        // if pivot is not acceptable, push it up to thresh
        // printf("small pivot %e ", Ajj);
        double old_pivot = Ajj;
        Ajj = thresh;

        // compute the minimum pivot required to keep the diagonal of the
        // current block acceptable:
        // b is column below pivot, d is diagonal of block, p is pivot
        // we want d_k - b_k^2 / p \ge thresh
        // i.e.
        // p \ge b_k^2 / (d_k - thresh)
        //
        double required_pivot = 0.0;
        for (int k = j + 1; k < n; ++k) {
          double bk = A[j + k * lda];
          double dk = A[k + lda * k];
          double temp = (dk - thresh);
          temp = (bk * bk) / temp;
          required_pivot = max(required_pivot, temp);
        }

        Ajj = max(Ajj, required_pivot);

        // record regularization used
        regul[j] += Ajj - old_pivot;

        // printf("set to %e\n", Ajj);
      }

      Ajj = sqrt(Ajj);
      A[j + lda * j] = Ajj;
      const double coeff = 1.0 / Ajj;

      // scale column j
      if (j < n - 1) {
        dscal_(&M, &coeff, &A[j + (j + 1) * lda], &lda);
      }
    }
  }

  return kRetOk;
}

double pivot_compensated_sum_2(double base, int k, const double* row, int ldr,
                               const double* piv, int ldp) {
  // Use Kahan-Babushka compensated summation, Neumaier version, to compute the
  // pivot. Compute base + \sum_0^k row[i*ldr]^2*piv[i*ldp]

  double sum = base;
  volatile double c = 0.0;

  for (int i = 0; i < k; ++i) {
    double value = -row[ldr * i] * piv[ldp * i] * row[ldr * i];
    double t = sum + value;
    if (fabs(sum) >= fabs(value)) {
      volatile double temp = sum - t;
      c += temp + value;
    } else {
      volatile double temp = value - t;
      c += temp + sum;
    }
    sum = t;
  }

  return sum + c;
}

double pivot_compensated_sum(double base, int k, const double* row, int ldr,
                             const double* piv, int ldp) {
  // Use Kahan-Babushka compensated summation to compute the pivot.
  // Compute base + \sum_0^k row[i*ldr]^2*piv[i*ldp]

  double sum = base;
  volatile double c = 0.0;

  for (int i = 0; i < k; ++i) {
    double value = -row[ldr * i] * piv[ldp * i] * row[ldr * i];
    double y = value - c;
    // use volatile or compiler optimization might use associativity to
    // eliminate operations
    volatile double t = sum + y;
    volatile double z = t - sum;
    c = z - y;
    sum = t;
  }

  return sum;
}

int dense_fact_fiuf(char uplo, int n, double* restrict A, int lda,
                    const int* pivot_sign, double thresh, double* regul) {
  // ===========================================================================
  // Infedinite factorization without blocks.
  // BLAS calls: dgemv_, dscal_.
  // ===========================================================================

  // check input
  if (n < 0 || !A || lda < n) {
    printf("\ndense_fact_fiuf: invalid input\n");
    return kRetInvalidInput;
  }

  // quick return
  if (n == 0) return kRetOk;

  // main operations
  if (uplo == 'L') {
    // allocate space for copy of col multiplied by pivots
    double* temp = malloc((n - 1) * sizeof(double));
    if (!temp) {
      printf("\ndense_fact_fiuf: out of memory\n");
      return kRetOutOfMemory;
    }

    for (int j = 0; j < n; ++j) {
      const int N = j;
      const int M = n - j - 1;

      // create temporary copy of row j, multiplied by pivots
      for (int i = 0; i < j; ++i) {
        temp[i] = A[j + i * lda] * A[i + i * lda];
      }

      // update diagonal element
      double Ajj =
          pivot_compensated_sum(A[j + lda * j], N, &A[j], lda, A, lda + 1);

      if (isnan(Ajj)) {
        A[j + lda * j] = Ajj;
        printf("\ndense_fact_fiuf: invalid pivot %e\n", Ajj);
        return kRetInvalidPivot;
      }

      // update column j
      if (j < n - 1) {
        dgemv_(&c_N, &M, &N, &d_m_one, &A[j + 1], &lda, temp, &i_one, &d_one,
               &A[j + 1 + j * lda], &i_one);
      }

      // sign of pivot
      double sign = (double)pivot_sign[j];
      double old_pivot = Ajj;

      // lift pivots that are too small or have wrong sign
      if (sign * Ajj <= thresh) {
        // if pivot is not acceptable, push it up to thresh
        printf("%d: small pivot %e, with sign %d set to ", j, Ajj,
               pivot_sign[j]);

        if (sign * Ajj <= -thresh * 1e3) {
          Ajj = sign * max(thresh * 1e100, 1e100);
          printf("%e\n", Ajj);
        } else {
          Ajj = thresh * sign;
          printf("%e\n", Ajj);

          // compute the minimum pivot required to keep the diagonal of the
          // current block acceptable:
          // b is column below pivot, d is diagonal of block, p is pivot
          // we want d_k - b_k^2 / p \ge thresh
          // i.e.
          // p \ge b_k^2 / (d_k - thresh)
          //
          double required_pivot = Ajj;
          for (int k = j + 1; k < n; ++k) {
            double bk = A[k + j * lda];
            double dk = A[k + k * lda];

            // if pivot and dk have different sign, skip
            if (sign * (double)pivot_sign[k] < 0) continue;

            double temp = (dk - sign * thresh);
            temp = (bk * bk) / temp;

            if (sign > 0)
              required_pivot = max(required_pivot, temp);
            else
              required_pivot = min(required_pivot, temp);
          }

          if (required_pivot != Ajj) {
            if (sign > 0)
              Ajj = max(Ajj, required_pivot);
            else
              Ajj = min(Ajj, required_pivot);

            printf("%d: pivot %e set to %e\n", j, old_pivot, Ajj);
          }
        }
      }

      // record regularization used
      regul[j] += fabs(Ajj - old_pivot);

      // save diagonal element
      A[j + lda * j] = Ajj;
      const double coeff = 1.0 / Ajj;

      // scale column j
      if (j < n - 1) {
        dscal_(&M, &coeff, &A[j + 1 + j * lda], &i_one);
      }
    }
    // free temporary copy of row
    free(temp);
  } else {
    // allocate space for copy of col multiplied by pivots
    double* temp = malloc((n - 1) * sizeof(double));
    if (!temp) {
      printf("\ndense_fact_fiuf: out of memory\n");
      return kRetOutOfMemory;
    }

    for (int j = 0; j < n; ++j) {
      const int N = j;
      const int M = n - j - 1;

      // create temporary copy of col j, multiplied by pivots
      for (int i = 0; i < j; ++i) {
        temp[i] = A[i + j * lda] * A[i + i * lda];
      }

      // update diagonal element
      double Ajj =
          pivot_compensated_sum(A[j + lda * j], N, &A[lda * j], 1, A, lda + 1);

      if (isnan(Ajj)) {
        A[j + lda * j] = Ajj;
        printf("\ndense_fact_fiuf: invalid pivot %e\n", Ajj);
        return kRetInvalidPivot;
      }

      // update column j
      if (j < n - 1) {
        dgemv_(&c_T, &N, &M, &d_m_one, &A[(j + 1) * lda], &lda, temp, &i_one,
               &d_one, &A[j + (j + 1) * lda], &lda);
      }

      // sign of pivot
      double sign = (double)pivot_sign[j];
      double old_pivot = Ajj;

      // lift pivots that are too small or have wrong sign
      if (sign * Ajj <= thresh) {
        // if pivot is not acceptable, push it up to thresh
        printf("%d: small pivot %e, with sign %d set to ", j, Ajj,
               pivot_sign[j]);

        if (sign * Ajj <= -thresh * 1e3) {
          Ajj = sign * max(thresh * 1e100, 1e100);
          printf("%e\n", Ajj);
        } else {
          Ajj = thresh * sign;
          printf("%e\n", Ajj);

          // compute the minimum pivot required to keep the diagonal of the
          // current block acceptable:
          // b is column below pivot, d is diagonal of block, p is pivot
          // we want d_k - b_k^2 / p \ge thresh
          // i.e.
          // p \ge b_k^2 / (d_k - thresh)
          //
          double required_pivot = Ajj;
          for (int k = j + 1; k < n; ++k) {
            double bk = A[j + k * lda];
            double dk = A[k + lda * k];

            // if pivot and dk have different sign, skip
            if (sign * (double)pivot_sign[k] < 0) continue;

            double temp = (dk - sign * thresh);
            temp = (bk * bk) / temp;

            if (sign > 0)
              required_pivot = max(required_pivot, temp);
            else
              required_pivot = min(required_pivot, temp);
          }

          if (required_pivot != Ajj) {
            if (sign > 0)
              Ajj = max(Ajj, required_pivot);
            else
              Ajj = min(Ajj, required_pivot);

            printf("%d: pivot %e set to %e\n", j, old_pivot, Ajj);
          }
        }
      }

      // record regularization used
      regul[j] += fabs(Ajj - old_pivot);

      // save diagonal element
      A[j + lda * j] = Ajj;
      const double coeff = 1.0 / Ajj;

      // scale column j
      if (j < n - 1) {
        dscal_(&M, &coeff, &A[j + (j + 1) * lda], &lda);
      }
    }

    // free temporary copy of row
    free(temp);
  }

  return kRetOk;
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
// In the following, block D represents the current diagonal block, in the block
// of columns that is being factorized. Block R is the rectangular block below
// D. Block P is the square or rectangular block to the left of D. Block Q is
// the rectangular block below P. See the report for a full explanation.
// ===========================================================================

int dense_fact_pdbf(int n, int k, int nb, double* restrict A, int lda,
                    double* restrict B, int ldb, double thresh, double* regul,
                    double* times) {
  // ===========================================================================
  // Positive definite factorization with blocks.
  // BLAS calls: dsyrk_, dgemm_, dtrsm_.
  // ===========================================================================

  // check input
  if (n < 0 || k < 0 || !A || lda < n || (k < n && (!B || ldb < n - k)) ||
      !times) {
    printf("\ndense_fact_pdbf: invalid input\n");
    return kRetInvalidInput;
  }

  // quick return
  if (n == 0) return kRetOk;

  // j is the starting col of the block of columns
  for (int j = 0; j < k; j += nb) {
    // jb is the size of the block
    const int jb = min(nb, k - j);

    // sizes for blas calls
    const int N = jb;
    const int K = j;
    const int M = n - j - jb;

    // starting position of matrices for BLAS calls
    double* D = &A[j + lda * j];
    const double* P = &A[j];
    const double* Q = &A[j + N];
    double* R = &A[j + N + lda * j];

    // update diagonal block
    callAndTime_dsyrk('L', 'N', N, K, -1.0, P, lda, 1.0, D, lda, times);

    // factorize diagonal block
    double* regul_current = &regul[j];
    int info = callAndTime_fduf('L', N, D, lda, thresh, regul_current, times);
    if (info != 0) return info;

    if (j + jb < n) {
      // update block of columns
      callAndTime_dgemm('N', 'T', M, N, K, -1.0, Q, lda, P, lda, 1.0, R, lda,
                        times);

      // solve block of columns with diagonal block
      callAndTime_dtrsm('R', 'L', 'T', 'N', M, N, 1.0, D, lda, R, lda, times);
    }
  }

  // update Schur complement if partial factorization is required
  if (k < n) {
    const int N = n - k;
    callAndTime_dsyrk('L', 'N', N, k, -1.0, &A[k], lda, 0.0, B, ldb, times);
  }

  return kRetOk;
}

int dense_fact_pibf(int n, int k, int nb, double* restrict A, int lda,
                    double* restrict B, int ldb, const int* pivot_sign,
                    double thresh, double* regul, double* times) {
  // ===========================================================================
  // Indefinite factorization with blocks.
  // BLAS calls: dcopy_, dscal_, dgemm_, dtrsm_, dsyrk_
  // ===========================================================================

  // check input
  if (n < 0 || k < 0 || !A || lda < n || (k < n && (!B || ldb < n - k))) {
    printf("\ndense_fact_pibf: invalid input\n");
    return kRetInvalidInput;
  }

  // quick return
  if (n == 0) return kRetOk;

  // j is the starting col of the block of columns
  for (int j = 0; j < k; j += nb) {
    // jb is the size of the block
    const int jb = min(nb, k - j);

    // sizes for blas calls
    const int N = jb;
    const int K = j;
    const int M = n - j - jb;

    // starting position of matrices for BLAS calls
    double* D = &A[j + lda * j];
    const double* P = &A[j];
    const double* Q = &A[j + N];
    double* R = &A[j + N + lda * j];

    // create temporary copy of block of rows, multiplied by pivots
    double* T = malloc(j * jb * sizeof(double));
    if (!T) {
      printf("\ndense_fact_pibf: out of memory\n");
      return kRetOutOfMemory;
    }
    int ldt = jb;
    for (int i = 0; i < j; ++i) {
      callAndTime_dcopy(N, &A[j + i * lda], 1, &T[i * ldt], 1, times);
      callAndTime_dscal(N, A[i + i * lda], &T[i * ldt], 1, times);
    }

    // update diagonal block using dgemm_
    callAndTime_dgemm('N', 'T', jb, jb, j, -1.0, P, lda, T, ldt, 1.0, D, lda,
                      times);

    // factorize diagonal block
    const int* pivot_sign_current = &pivot_sign[j];
    double* regul_current = &regul[j];
    int info = callAndTime_fiuf('L', N, D, lda, pivot_sign_current, thresh,
                                regul_current, times);
    if (info != 0) return info;

    if (j + jb < n) {
      // update block of columns
      callAndTime_dgemm('N', 'T', M, N, K, -1.0, Q, lda, T, ldt, 1.0, R, lda,
                        times);

      // solve block of columns with L
      callAndTime_dtrsm('R', 'L', 'T', 'U', M, N, 1.0, D, lda, R, lda, times);

      // solve block of columns with D
      for (int i = 0; i < jb; ++i) {
        const double coeff = 1.0 / A[j + i + (j + i) * lda];
        callAndTime_dscal(M, coeff, &A[j + jb + lda * (j + i)], 1, times);
      }
    }

    free(T);
  }

  // update Schur complement
  if (k < n) {
    const int N = n - k;

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
      printf("\ndense_fact_pibf: out of memory\n");
      return kRetOutOfMemory;
    }
    double* temp_neg = malloc((n - k) * neg_pivot * sizeof(double));
    if (!temp_neg) {
      printf("\ndense_fact_pibf: out of memory\n");
      return kRetOutOfMemory;
    }

    const int ldt = n - k;

    // the copies of the columns are multiplied by sqrt(|Ajj|)
    int start_pos = 0;
    int start_neg = 0;
    for (int j = 0; j < k; ++j) {
      double Ajj = A[j + lda * j];
      if (Ajj >= 0.0) {
        Ajj = sqrt(Ajj);
        callAndTime_dcopy(N, &A[k + j * lda], 1, &temp_pos[start_pos * ldt], 1,
                          times);
        callAndTime_dscal(N, Ajj, &temp_pos[start_pos * ldt], 1, times);
        ++start_pos;
      } else {
        Ajj = sqrt(-Ajj);
        callAndTime_dcopy(N, &A[k + j * lda], 1, &temp_neg[start_neg * ldt], 1,
                          times);
        callAndTime_dscal(N, Ajj, &temp_neg[start_neg * ldt], 1, times);
        ++start_neg;
      }
    }

    // Update schur complement by subtracting contribution of positive columns
    // and adding contribution of negative columns.
    // In this way, I can use dsyrk_ instead of dgemm_ and avoid updating the
    // full square schur complement. First call uses beta = 0.0, to clear
    // content of B. Second call uses beta = 1.0, to not clear the result of
    // the first call.

    callAndTime_dsyrk('L', 'N', N, pos_pivot, -1.0, temp_pos, ldt, 0.0, B, ldb,
                      times);
    callAndTime_dsyrk('L', 'N', N, neg_pivot, 1.0, temp_neg, ldt, 1.0, B, ldb,
                      times);

    free(temp_pos);
    free(temp_neg);
  }

  return kRetOk;
}

int dense_fact_pdbh(int n, int k, int nb, double* restrict A,
                    double* restrict B, double thresh, double* regul,
                    double* times) {
  // ===========================================================================
  // Positive definite factorization with blocks in lower-blocked-hybrid
  // format. A should be in lower-blocked-hybrid format. Schur complement is
  // returned in B in lower packed format (not lower-blocked-hybrid).
  // BLAS calls: dsyrk_, dgemm_, dtrsm_, dcopy_
  // ===========================================================================

  // check input
  if (n < 0 || k < 0 || !A || (k < n && !B)) {
    printf("\ndense_fact_pdbh: invalid input\n");
    return kRetInvalidInput;
  }

  // quick return
  if (n == 0) return kRetOk;

  // number of blocks of columns
  const int n_blocks = (k - 1) / nb + 1;

  // start of diagonal blocks
  int* diag_start = malloc(n_blocks * sizeof(double));
  if (!diag_start) {
    printf("\ndense_fact_pdbh: out of memory\n");
    return kRetOutOfMemory;
  }
  diag_start[0] = 0;
  for (int i = 1; i < n_blocks; ++i) {
    diag_start[i] =
        diag_start[i - 1] + nb * (2 * n - 2 * (i - 1) * nb - nb + 1) / 2;
  }

  // size of blocks
  const int diag_size = nb * (nb + 1) / 2;
  const int full_size = nb * nb;

  // variables for blas calls
  int info;

  // buffer for full-format diagonal blocks
  double* D = malloc(nb * nb * sizeof(double));
  if (!D) {
    printf("\ndense_fact_pdbh: out of memory\n");
    return kRetOutOfMemory;
  }

  // j is the index of the block column
  for (int j = 0; j < n_blocks; ++j) {
    // jb is the number of columns
    const int jb = min(nb, k - nb * j);

    // size of current block could be smaller than diag_size and full_size
    const int this_diag_size = jb * (jb + 1) / 2;
    const int this_full_size = nb * jb;

    // full copy of diagonal block by rows, in D
    int offset = 0;
    for (int Drow = 0; Drow < jb; ++Drow) {
      const int N = Drow + 1;
      callAndTime_dcopy(N, &A[diag_start[j] + offset], 1, &D[Drow * jb], 1,
                        times);
      offset += N;
    }

    // number of rows left below block j
    const int M = n - nb * j - jb;

    // block of columns below diagonal block j
    double* R = &A[diag_start[j] + this_diag_size];

    // update diagonal block and block of columns
    for (int k = 0; k < j; ++k) {
      // starting position of block to input to dsyrk_
      int Pk_pos = diag_start[k] + diag_size;
      if (j > k + 1) Pk_pos += full_size * (j - k - 1);
      const double* Pk = &A[Pk_pos];

      callAndTime_dsyrk('U', 'T', jb, nb, -1.0, Pk, nb, 1.0, D, jb, times);

      if (M > 0) {
        const int Qk_pos = Pk_pos + this_full_size;
        const double* Qk = &A[Qk_pos];

        callAndTime_dgemm('T', 'N', jb, M, nb, -1.0, Pk, nb, Qk, nb, 1.0, R, jb,
                          times);
      }
    }

    // factorize diagonal block
    double* regul_current = &regul[j * nb];
    int info = callAndTime_fduf('U', jb, D, jb, thresh, regul_current, times);
    if (info != 0) return info;

    if (M > 0) {
      // solve block of columns with diagonal block
      callAndTime_dtrsm('L', 'U', 'T', 'N', jb, M, 1.0, D, jb, R, jb, times);
    }

    // put D back into packed format
    offset = 0;
    for (int Drow = 0; Drow < jb; ++Drow) {
      const int N = Drow + 1;
      callAndTime_dcopy(N, &D[Drow * jb], 1, &A[diag_start[j] + offset], 1,
                        times);
      offset += N;
    }
  }
  free(D);

  // compute Schur complement if partial factorization is required
  if (k < n) {
    // number of rows/columns in the Schur complement
    const int ns = n - k;

    // size of last full block (may be smaller than full_size)
    const int ncol_last = k % nb;
    const int last_full_size = ncol_last == 0 ? full_size : ncol_last * nb;

    double beta = 0.0;

    // buffer for full-format of block of columns of Schur complement
    double* schur_buf = malloc(ns * nb * sizeof(double));
    if (!schur_buf) {
      printf("\ndense_fact_pdbh: out of memory\n");
      return kRetOutOfMemory;
    }

    // number of blocks in Schur complement
    const int s_blocks = (ns - 1) / nb + 1;

    // index to write into B
    int B_start = 0;

    // Go through block of columns of Schur complement.
    // Each block will have a full-format copy made in schur_buf
    for (int sb = 0; sb < s_blocks; ++sb) {
      // number of rows of the block
      const int nrow = ns - nb * sb;

      // number of columns of the block
      const int ncol = min(nb, nrow);

      // number of elements in diagonal block
      const int schur_diag_size = ncol * (ncol + 1) / 2;

      beta = 0.0;

      // each block receives contributions from the blocks of the leading part
      // of A
      for (int j = 0; j < n_blocks; ++j) {
        const int jb = min(nb, k - nb * j);
        const int this_diag_size = jb * (jb + 1) / 2;
        const int this_full_size = nb * jb;

        // compute index to access diagonal block in A
        int diag_pos = diag_start[j] + this_diag_size;
        if (j < n_blocks - 1) {
          diag_pos += (n_blocks - j - 2) * full_size + last_full_size;
        }
        diag_pos += sb * this_full_size;

        // update diagonal block
        callAndTime_dsyrk('U', 'T', ncol, jb, -1.0, &A[diag_pos], jb, beta,
                          schur_buf, ncol, times);

        // update subdiagonal part
        const int M = nrow - nb;
        if (M > 0) {
          callAndTime_dgemm('T', 'N', nb, M, jb, -1.0, &A[diag_pos], jb,
                            &A[diag_pos + this_full_size], jb, beta,
                            &schur_buf[ncol * ncol], ncol, times);
        }

        // beta is 0 for the first time (to avoid initializing schur_buf) and
        // 1 for the next calls
        beta = 1.0;
      }

      // schur_buf contains Schur complement in hybrid format (with full
      // diagonal blocks). Put it in lower-packed format in B.

      for (int buf_row = 0; buf_row < nrow; ++buf_row) {
        const int N = ncol;
        callAndTime_dcopy_schur(N, &schur_buf[buf_row * ncol], 1,
                                &B[B_start + buf_row], nrow, times);
      }
      B_start += nrow * ncol;
    }

    free(schur_buf);
  }

  free(diag_start);

  return kRetOk;
}

int dense_fact_pibh(int n, int k, int nb, double* restrict A,
                    double* restrict B, const int* pivot_sign, double thresh,
                    double* regul, double* times) {
  // ===========================================================================
  // Indefinite factorization with blocks in lower-blocked-hybrid format.
  // A should be in lower-blocked-hybrid format. Schur complement is returned
  // in B in lower packed format (not lower-blocked-hybrid). BLAS calls:
  // dsyrk_, dgemm_, dtrsm_, dcopy_, dscal_
  // ===========================================================================

  const int sizeA = n * k - k * (k - 1) / 2;

  // check input
  if (n < 0 || k < 0 || !A || (k < n && !B)) {
    printf("\ndense_fact_pibh: invalid input\n");
    return kRetInvalidInput;
  }

  // quick return
  if (n == 0) return kRetOk;

  // number of blocks of columns
  const int n_blocks = (k - 1) / nb + 1;

  // start of diagonal blocks
  int* diag_start = malloc(n_blocks * sizeof(double));
  if (!diag_start) {
    printf("\ndense_fact_pibh: out of memory\n");
    return kRetOutOfMemory;
  }
  diag_start[0] = 0;
  for (int i = 1; i < n_blocks; ++i) {
    diag_start[i] =
        diag_start[i - 1] + nb * (2 * n - 2 * (i - 1) * nb - nb + 1) / 2;
  }

  // size of blocks
  const int diag_size = nb * (nb + 1) / 2;
  const int full_size = nb * nb;

  // variables for blas calls
  int info;

  // buffer for full-format diagonal blocks
  double* D = malloc(nb * nb * sizeof(double));
  if (!D) {
    printf("\ndense_fact_pibh: out of memory\n");
    return kRetOutOfMemory;
  }

  // buffer for copy of block scaled by pivots
  double* T = malloc(nb * nb * sizeof(double));
  if (!T) {
    printf("\ndense_fact_pibh: out of memory\n");
    return kRetOutOfMemory;
  }

  // j is the index of the block column
  for (int j = 0; j < n_blocks; ++j) {
    // jb is the number of columns
    const int jb = min(nb, k - nb * j);

    // size of current block could be smaller than diag_size and full_size
    const int this_diag_size = jb * (jb + 1) / 2;
    const int this_full_size = nb * jb;

    // full copy of diagonal block by rows, in D
    int offset = 0;
    for (int Drow = 0; Drow < jb; ++Drow) {
      const int N = Drow + 1;
      callAndTime_dcopy(N, &A[diag_start[j] + offset], 1, &D[Drow * jb], 1,
                        times);
      offset += N;
    }

    // number of rows left below block j
    const int M = n - nb * j - jb;

    // block of columns below diagonal block j
    const int R_pos = diag_start[j] + this_diag_size;
    double* R = &A[R_pos];

    // update diagonal block and block of columns
    for (int k = 0; k < j; ++k) {
      // starting position of block to input to dsyrk_
      int Pk_pos = diag_start[k] + diag_size;
      if (j > k + 1) Pk_pos += full_size * (j - k - 1);
      const double* Pk = &A[Pk_pos];

      // copy block jk into temp
      callAndTime_dcopy(this_full_size, Pk, 1, T, 1, times);

      // scale temp by pivots
      int pivot_pos = diag_start[k];
      for (int col = 0; col < nb; ++col) {
        callAndTime_dscal(jb, A[pivot_pos], &T[col], nb, times);
        pivot_pos += col + 2;
      }

      // update diagonal block with dgemm_
      callAndTime_dgemm('T', 'N', jb, jb, nb, -1.0, T, nb, Pk, nb, 1.0, D, jb,
                        times);

      // update rectangular block
      if (M > 0) {
        const int Qk_pos = Pk_pos + this_full_size;
        const double* Qk = &A[Qk_pos];
        callAndTime_dgemm('T', 'N', jb, M, nb, -1.0, T, nb, Qk, nb, 1.0, R, jb,
                          times);
      }
    }

    // factorize diagonal block
    double* regul_current = &regul[j * nb];
    const int* pivot_sign_current = &pivot_sign[j * nb];
    int info = callAndTime_fiuf('U', jb, D, jb, pivot_sign_current, thresh,
                                regul_current, times);
    if (info != 0) return info;

    if (M > 0) {
      // solve block of columns with diagonal block
      callAndTime_dtrsm('L', 'U', 'T', 'U', jb, M, 1.0, D, jb, R, jb, times);

      // scale columns by pivots
      int pivot_pos = 0;
      for (int col = 0; col < jb; ++col) {
        const double coeff = 1.0 / D[pivot_pos];
        callAndTime_dscal(M, coeff, &A[R_pos + col], jb, times);
        pivot_pos += jb + 1;
      }
    }

    // put D back into packed format
    offset = 0;
    for (int Drow = 0; Drow < jb; ++Drow) {
      const int N = Drow + 1;
      callAndTime_dcopy(N, &D[Drow * jb], 1, &A[diag_start[j] + offset], 1,
                        times);
      offset += N;
    }
  }
  free(D);

  // compute Schur complement if partial factorization is required
  if (k < n) {
    // number of rows/columns in the Schur complement
    const int ns = n - k;

    // size of last full block (may be smaller than full_size)
    const int ncol_last = k % nb;
    const int last_full_size = ncol_last == 0 ? full_size : ncol_last * nb;

    double beta = 0.0;

    // buffer for full-format of block of columns of Schur complement
    double* schur_buf = malloc(ns * nb * sizeof(double));
    if (!schur_buf) {
      printf("\ndense_fact_pibh: out of memory\n");
      return kRetOutOfMemory;
    }

    // number of blocks in Schur complement
    const int s_blocks = (ns - 1) / nb + 1;

    // index to write into B
    int B_start = 0;

    // Go through block of columns of Schur complement.
    // Each block will have a full-format copy made in schur_buf
    for (int sb = 0; sb < s_blocks; ++sb) {
      // number of rows of the block
      const int nrow = ns - nb * sb;

      // number of columns of the block
      const int ncol = min(nb, nrow);

      // number of elements in diagonal block
      const int schur_diag_size = ncol * (ncol + 1) / 2;

      beta = 0.0;

      // each block receives contributions from the blocks of the leading part
      // of A
      for (int j = 0; j < n_blocks; ++j) {
        const int jb = min(nb, k - nb * j);
        const int this_diag_size = jb * (jb + 1) / 2;
        const int this_full_size = nb * jb;

        // compute index to access diagonal block in A
        int diag_pos = diag_start[j] + this_diag_size;
        if (j < n_blocks - 1) {
          diag_pos += (n_blocks - j - 2) * full_size + last_full_size;
        }
        diag_pos += sb * this_full_size;

        // create copy of block, multiplied by pivots
        const int N = ncol * jb;
        callAndTime_dcopy(N, &A[diag_pos], 1, T, 1, times);

        int pivot_pos = diag_start[j];

        for (int col = 0; col < jb; ++col) {
          callAndTime_dscal(ncol, A[pivot_pos], &T[col], jb, times);
          pivot_pos += col + 2;
        }

        // update diagonal block using dgemm_
        // printf("T: size %d, access %d\n", nb * nb, jb * ncol);
        // printf("A: size %d, access %d\n", sizeA, diag_pos + ncol * jb);
        // printf("S: size %d, acceess %d\n", ns * nb, ncol * ncol);

        callAndTime_dgemm('T', 'N', ncol, ncol, jb, -1.0, &A[diag_pos], jb, T,
                          jb, beta, schur_buf, ncol, times);

        // update subdiagonal part
        const int M = nrow - nb;
        if (M > 0) {
          callAndTime_dgemm('T', 'N', ncol, M, jb, -1.0, T, jb,
                            &A[diag_pos + this_full_size], jb, beta,
                            &schur_buf[ncol * ncol], ncol, times);
        }

        // beta is 0 for the first time (to avoid initializing schur_buf) and
        // 1 for the next calls
        beta = 1.0;
      }

      // schur_buf contains Schur complement in hybrid format (with full
      // diagonal blocks). Put it in lower-packed format in B.
      for (int buf_row = 0; buf_row < nrow; ++buf_row) {
        const int N = ncol;
        callAndTime_dcopy_schur(N, &schur_buf[buf_row * ncol], 1,
                                &B[B_start + buf_row], nrow, times);
      }
      B_start += nrow * ncol;
    }

    free(schur_buf);
  }

  free(T);
  free(diag_start);

  return kRetOk;
}

int dense_fact_pdbs(int n, int k, int nb, double* A, double* B, double thresh,
                    double* regul, double* times) {
  // ===========================================================================
  // Positive definite factorization with blocks in lower-blocked-hybrid
  // format. A should be in lower-blocked-hybrid format. Schur complement is
  // returned in B in lower-block-hybrid format, with full diagonal blocks.
  // BLAS calls: dsyrk_, dgemm_, dtrsm_, dcopy_
  // ===========================================================================

  // check input
  if (n < 0 || k < 0 || !A || (k < n && !B)) {
    printf("\ndense_fact_pdbs: invalid input\n");
    return kRetInvalidInput;
  }

  // quick return
  if (n == 0) return kRetOk;

  // number of blocks of columns
  const int n_blocks = (k - 1) / nb + 1;

  // start of diagonal blocks
  int* diag_start = malloc(n_blocks * sizeof(double));
  if (!diag_start) {
    printf("\ndense_fact_pdbs: out of memory\n");
    return kRetOutOfMemory;
  }
  diag_start[0] = 0;
  for (int i = 1; i < n_blocks; ++i) {
    diag_start[i] =
        diag_start[i - 1] + nb * (2 * n - 2 * (i - 1) * nb - nb + 1) / 2;
  }

  // size of blocks
  const int diag_size = nb * (nb + 1) / 2;
  const int full_size = nb * nb;

  // variables for blas calls
  int info;

  // buffer for full-format diagonal blocks
  double* D = malloc(nb * nb * sizeof(double));
  if (!D) {
    printf("\ndense_fact_pdbs: out of memory\n");
    return kRetOutOfMemory;
  }

  // j is the index of the block column
  for (int j = 0; j < n_blocks; ++j) {
    // jb is the number of columns
    const int jb = min(nb, k - nb * j);

    // size of current block could be smaller than diag_size and full_size
    const int this_diag_size = jb * (jb + 1) / 2;
    const int this_full_size = nb * jb;

    // full copy of diagonal block by rows, in D
    int offset = 0;
    for (int Drow = 0; Drow < jb; ++Drow) {
      const int N = Drow + 1;
      callAndTime_dcopy(N, &A[diag_start[j] + offset], 1, &D[Drow * jb], 1,
                        times);
      offset += N;
    }

    // number of rows left below block j
    const int M = n - nb * j - jb;

    // block of columns below diagonal block j
    double* R = &A[diag_start[j] + this_diag_size];

    // update diagonal block and block of columns
    for (int k = 0; k < j; ++k) {
      // starting position of block to input to dsyrk_
      int Pk_pos = diag_start[k] + diag_size;
      if (j > k + 1) Pk_pos += full_size * (j - k - 1);
      const double* Pk = &A[Pk_pos];

      callAndTime_dsyrk('U', 'T', jb, nb, -1.0, Pk, nb, 1.0, D, jb, times);

      if (M > 0) {
        const int Qk_pos = Pk_pos + this_full_size;
        const double* Qk = &A[Qk_pos];
        callAndTime_dgemm('T', 'N', jb, M, nb, -1.0, Pk, nb, Qk, nb, 1.0, R, jb,
                          times);
      }
    }

    // factorize diagonal block
    double* regul_current = &regul[j * nb];
    int info = callAndTime_fduf('U', jb, D, jb, thresh, regul_current, times);
    if (info != 0) return info;

    if (M > 0) {
      // solve block of columns with diagonal block
      callAndTime_dtrsm('L', 'U', 'T', 'N', jb, M, 1.0, D, jb, R, jb, times);
    }

    // put D back into packed format
    offset = 0;
    for (int Drow = 0; Drow < jb; ++Drow) {
      const int N = Drow + 1;
      callAndTime_dcopy(N, &D[Drow * jb], 1, &A[diag_start[j] + offset], 1,
                        times);
      offset += N;
    }
  }
  free(D);

  // compute Schur complement if partial factorization is required
  if (k < n) {
    // number of rows/columns in the Schur complement
    const int ns = n - k;

    // size of last full block (may be smaller than full_size)
    const int ncol_last = k % nb;
    const int last_full_size = ncol_last == 0 ? full_size : ncol_last * nb;

    double beta = 0.0;

    // number of blocks in Schur complement
    const int s_blocks = (ns - 1) / nb + 1;

    // index to write into B
    int B_start = 0;

    // Go through block of columns of Schur complement.
    // Each block will have a full-format copy made in schur_buf
    for (int sb = 0; sb < s_blocks; ++sb) {
      // location in B where to write the next block of columns
      double* schur_buf = &B[B_start];

      // number of rows of the block
      const int nrow = ns - nb * sb;

      // number of columns of the block
      const int ncol = min(nb, nrow);

      beta = 0.0;

      // each block receives contributions from the blocks of the leading part
      // of A
      for (int j = 0; j < n_blocks; ++j) {
        const int jb = min(nb, k - nb * j);
        const int this_diag_size = jb * (jb + 1) / 2;
        const int this_full_size = nb * jb;

        // compute index to access diagonal block in A
        int diag_pos = diag_start[j] + this_diag_size;
        if (j < n_blocks - 1) {
          diag_pos += (n_blocks - j - 2) * full_size + last_full_size;
        }
        diag_pos += sb * this_full_size;

        // update diagonal block
        callAndTime_dsyrk('U', 'T', ncol, jb, -1.0, &A[diag_pos], jb, beta,
                          schur_buf, ncol, times);

        // update subdiagonal part
        const int M = nrow - nb;
        if (M > 0) {
          callAndTime_dgemm('T', 'N', nb, M, jb, -1.0, &A[diag_pos], jb,
                            &A[diag_pos + this_full_size], jb, beta,
                            &schur_buf[ncol * ncol], ncol, times);
        }

        // beta is 0 for the first time (to avoid initializing schur_buf) and
        // 1 for the next calls
        beta = 1.0;
      }

      B_start += nrow * ncol;
    }
  }

  free(diag_start);

  return kRetOk;
}

int dense_fact_pibs(int n, int k, int nb, double* A, double* B,
                    const int* pivot_sign, double thresh, double* regul,
                    double* times) {
  // ===========================================================================
  // Indefinite factorization with blocks in lower-blocked-hybrid format.
  // A should be in lower-blocked-hybrid format. Schur complement is returned
  // in B in in lower-block-hybrid format, with full diagonal blocks.
  // BLAS calls: dsyrk_, dgemm_, dtrsm_, dcopy_, dscal_
  // ===========================================================================

  // check input
  if (n < 0 || k < 0 || !A || (k < n && !B)) {
    printf("\ndense_fact_pibs: invalid input\n");
    return kRetInvalidInput;
  }

  // quick return
  if (n == 0) return kRetOk;

  // number of blocks of columns
  const int n_blocks = (k - 1) / nb + 1;

  // start of diagonal blocks
  int* diag_start = malloc(n_blocks * sizeof(double));
  if (!diag_start) {
    printf("\ndense_fact_pibs: out of memory\n");
    return kRetOutOfMemory;
  }
  diag_start[0] = 0;
  for (int i = 1; i < n_blocks; ++i) {
    diag_start[i] =
        diag_start[i - 1] + nb * (2 * n - 2 * (i - 1) * nb - nb + 1) / 2;
  }

  // size of blocks
  const int diag_size = nb * (nb + 1) / 2;
  const int full_size = nb * nb;

  // variables for blas calls
  int info;

  // buffer for full-format diagonal blocks
  double* D = malloc(nb * nb * sizeof(double));
  if (!D) {
    printf("\ndense_fact_pibs: out of memory\n");
    return kRetOutOfMemory;
  }

  // buffer for copy of block scaled by pivots
  double* T = malloc(nb * nb * sizeof(double));
  if (!T) {
    printf("\ndense_fact_pibs: out of memory\n");
    return kRetOutOfMemory;
  }

  // j is the index of the block column
  for (int j = 0; j < n_blocks; ++j) {
    // jb is the number of columns
    const int jb = min(nb, k - nb * j);

    // size of current block could be smaller than diag_size and full_size
    const int this_diag_size = jb * (jb + 1) / 2;
    const int this_full_size = nb * jb;

    // full copy of diagonal block by rows, in D
    int offset = 0;
    for (int Drow = 0; Drow < jb; ++Drow) {
      const int N = Drow + 1;
      callAndTime_dcopy(N, &A[diag_start[j] + offset], 1, &D[Drow * jb], 1,
                        times);
      offset += N;
    }

    // number of rows left below block j
    const int M = n - nb * j - jb;

    // starting position of block of columns below diagonal block j
    const int R_pos = diag_start[j] + this_diag_size;
    double* R = &A[R_pos];

    // update diagonal block and block of columns
    for (int k = 0; k < j; ++k) {
      // starting position of block to input to dsyrk_
      int Pk_pos = diag_start[k] + diag_size;
      if (j > k + 1) Pk_pos += full_size * (j - k - 1);
      const double* Pk = &A[Pk_pos];

      // copy block jk into temp
      callAndTime_dcopy(this_full_size, Pk, 1, T, 1, times);

      // scale temp by pivots
      int pivot_pos = diag_start[k];
      for (int col = 0; col < nb; ++col) {
        callAndTime_dscal(jb, A[pivot_pos], &T[col], nb, times);
        pivot_pos += col + 2;
      }

      // update diagonal block with dgemm_
      callAndTime_dgemm('T', 'N', jb, jb, nb, -1.0, T, nb, Pk, nb, 1.0, D, jb,
                        times);

      // update rectangular block
      if (M > 0) {
        const int Qk_pos = Pk_pos + this_full_size;
        double* Qk = &A[Qk_pos];

        callAndTime_dgemm('T', 'N', jb, M, nb, -1.0, T, nb, Qk, nb, 1.0, R, jb,
                          times);
      }
    }

    // factorize diagonal block
    const int* pivot_sign_current = &pivot_sign[j * nb];
    double* regul_current = &regul[j * nb];
    int info = callAndTime_fiuf('U', jb, D, jb, pivot_sign_current, thresh,
                                regul_current, times);
    if (info != 0) return info;

    if (M > 0) {
      // solve block of columns with diagonal block
      callAndTime_dtrsm('L', 'U', 'T', 'U', jb, M, 1.0, D, jb, R, jb, times);

      // scale columns by pivots
      int pivot_pos = 0;
      for (int col = 0; col < jb; ++col) {
        double coeff = 1.0 / D[pivot_pos];
        callAndTime_dscal(M, coeff, &A[R_pos + col], jb, times);
        pivot_pos += jb + 1;
      }
    }

    // put D back into packed format
    offset = 0;
    for (int Drow = 0; Drow < jb; ++Drow) {
      const int N = Drow + 1;
      callAndTime_dcopy(N, &D[Drow * jb], 1, &A[diag_start[j] + offset], 1,
                        times);
      offset += N;
    }
  }
  free(D);

  // compute Schur complement if partial factorization is required
  if (k < n) {
    // number of rows/columns in the Schur complement
    const int ns = n - k;

    // size of last full block (may be smaller than full_size)
    const int ncol_last = k % nb;
    const int last_full_size = ncol_last == 0 ? full_size : ncol_last * nb;

    double beta = 0.0;

    // number of blocks in Schur complement
    const int s_blocks = (ns - 1) / nb + 1;

    // index to write into B
    int B_start = 0;

    // Go through block of columns of Schur complement.
    // Each block will have a full-format copy made in schur_buf
    for (int sb = 0; sb < s_blocks; ++sb) {
      // location in B where to write the next block of columns
      double* schur_buf = &B[B_start];

      // number of rows of the block
      const int nrow = ns - nb * sb;

      // number of columns of the block
      const int ncol = min(nb, nrow);

      beta = 0.0;

      // each block receives contributions from the blocks of the leading part
      // of A
      for (int j = 0; j < n_blocks; ++j) {
        const int jb = min(nb, k - nb * j);
        const int this_diag_size = jb * (jb + 1) / 2;
        const int this_full_size = nb * jb;

        // compute index to access diagonal block in A
        int diag_pos = diag_start[j] + this_diag_size;
        if (j < n_blocks - 1) {
          diag_pos += (n_blocks - j - 2) * full_size + last_full_size;
        }
        diag_pos += sb * this_full_size;

        // create copy of block, multiplied by pivots
        const int N = ncol * jb;
        callAndTime_dcopy(N, &A[diag_pos], 1, T, 1, times);

        int pivot_pos = diag_start[j];
        for (int col = 0; col < jb; ++col) {
          callAndTime_dscal(ncol, A[pivot_pos], &T[col], jb, times);
          pivot_pos += col + 2;
        }

        // update diagonal block using dgemm_
        callAndTime_dgemm('T', 'N', ncol, ncol, jb, -1.0, T, jb, &A[diag_pos],
                          jb, beta, schur_buf, ncol, times);

        // update subdiagonal part
        const int M = nrow - nb;
        if (M > 0) {
          callAndTime_dgemm('T', 'N', ncol, M, jb, -1.0, T, jb,
                            &A[diag_pos + this_full_size], jb, beta,
                            &schur_buf[ncol * ncol], ncol, times);
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

  return kRetOk;
}

int dense_fact_l2h(double* restrict A, int nrow, int ncol, int nb,
                   double* times) {
  // ===========================================================================
  // Takes a matrix in lower-packed format, with nrow rows.
  // Converts the first ncol columns into lower-blocked-hybrid format, with
  // block size nb.
  // BLAS calls: dcopy_
  // ===========================================================================

#ifdef FINEST_TIMING
  double t0 = GetTime();
#endif

  double* buf = malloc(nrow * nb * sizeof(double));
  if (!buf) {
    printf("\ndense_fact_l2h: out of memory\n");
    return kRetOutOfMemory;
  }
  int startAtoBuf = 0;
  int startBuftoA = 0;

  for (int k = 0; k <= (ncol - 1) / nb; ++k) {
    // Each block of columns is copied into the buffer, leaving space in
    // between columns, to align rows. This "empty space" contains elements
    // from the previous use of buffer, but they are ignored when copying
    // back.

    // Number of columns in the block. Can be smaller than nb for last block.
    const int block_size = min(nb, ncol - k * nb);

    // Number of rows in the block
    const int row_size = nrow - k * nb;

    int startBuf = 0;
    for (int j = 0; j < block_size; ++j) {
      // Copy each column in the block
      const int N = row_size - j;
      dcopy_(&N, &A[startAtoBuf], &i_one, &buf[startBuf], &i_one);
      startAtoBuf += N;

      // +1 to align rows
      startBuf += row_size + 1;
    }

    // Copy columns back into A, row by row.
    // One call of dcopy_ for each row of the block of columns.
    for (int i = 0; i < row_size; ++i) {
      const int N = min(i + 1, block_size);
      dcopy_(&N, &buf[i], &row_size, &A[startBuftoA], &i_one);
      startBuftoA += N;
    }
  }

  free(buf);

#ifdef FINEST_TIMING
  times[kTimeDenseFact_convert] += GetTime() - t0;
#endif

  return kRetOk;
}