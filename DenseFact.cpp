#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <vector>

#include "../ProtoIPM/Regularization.h"
#include "Blas_declaration.h"
#include "CallAndTimeBlas.h"
#include "DataCollector.h"
#include "ReturnValues.h"
#include "Timing.h"

/*
Names:
denseFact:
 - K : factorization kernel for diagonal blocks
 - F : blocked factorization in full format
 - HP: blocked factorization in hybrid-packed format
 - HH: blocked factorization in hybrid-hybrid format
*/

double pivotCompensatedSum(double base, int k, const double* row, int ldr,
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

double regularizePivot(double pivot, double thresh, const int* sign,
                       const double* A, int lda, int j, int n, char uplo,
                       DataCollector& DC) {
  // add static regularization
  if (sign[j] == 1)
    pivot += kDualStaticRegularization;
  else
    pivot -= kPrimalStaticRegularization;

  double s = (double)sign[j];
  double old_pivot = pivot;

  double spivot = s * pivot;
  double K = 1e12;

  if (spivot <= thresh && spivot >= -thresh) {
    // small pivot, lift to thresh
    pivot = s * thresh;
    DC.sumRegPiv();
    printf("%3d: small pivot %e, with sign %d, set to %e\n", j, old_pivot,
           sign[j], pivot);
  } else if (spivot < -thresh && spivot >= -thresh * K) {
    // wrong sign, lift more
    pivot = s * thresh * 10;
    DC.sumRegPiv();
    printf("%3d: wrong pivot %e, with sign %d, set to %e\n", j, old_pivot,
           sign[j], pivot);
  } else if (spivot < -thresh * K) {
    // pivot is completely lost
    pivot = s * 1e100;
    DC.sumRegPiv();
    printf("%3d: disaster pivot %e, with sign %d, set to %e\n", j, old_pivot,
           sign[j], pivot);
  }

  // lift pivots that are too small or have wrong sign
  /*if (s * pivot <= thresh) {
    // if pivot is not acceptable, push it up to thresh
    printf("%d: small pivot %e, with sign %d set to ", j, pivot, sign[j]);

    if (s * pivot <= -thresh * 1e12) {
      pivot = s * max(thresh * 1e10, 1e10);
      printf("%e\n", pivot);
    } else {
      pivot = thresh * s;
      printf("%e\n", pivot);

      // compute the minimum pivot required to keep the diagonal of the
      // current block acceptable:
      // b is column below pivot, d is diagonal of block, p is pivot
      // we want d_k - b_k^2 / p \ge thresh
      // i.e.
      // p \ge b_k^2 / (d_k - thresh)
      //
      double required_pivot = pivot;
      for (int k = j + 1; k < n; ++k) {
        double bk = uplo == 'L' ? A[k + j * lda] : A[j + k * lda];
        double dk = A[k + k * lda];

        // if pivot and dk have different sign, skip
        if (s * (double)sign[k] < 0) continue;

        double temp = (dk - s * thresh);
        temp = (bk * bk) / temp;

        if (s > 0)
          required_pivot = max(required_pivot, temp);
        else
          required_pivot = min(required_pivot, temp);
      }

      if (required_pivot != pivot) {
        if (sign > 0)
          pivot = max(pivot, required_pivot);
        else
          pivot = min(pivot, required_pivot);

        printf("%d: pivot %e set to %e\n", j, old_pivot, pivot);
      }
    }
  }*/

  return pivot;
}

int denseFactK(char uplo, int n, double* A, int lda, const int* pivot_sign,
               double thresh, double* regul, DataCollector& DC) {
  // ===========================================================================
  // Factorization kernel.
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
    std::vector<double> temp(n - 1);

    for (int j = 0; j < n; ++j) {
      const int N = j;
      const int M = n - j - 1;

      // create temporary copy of row j, multiplied by pivots
      for (int i = 0; i < j; ++i) {
        temp[i] = A[j + i * lda] * A[i + i * lda];
      }

      // update diagonal element
      double Ajj =
          pivotCompensatedSum(A[j + lda * j], N, &A[j], lda, A, lda + 1);

      if (isnan(Ajj)) {
        A[j + lda * j] = Ajj;
        printf("\ndense_fact_fiuf: invalid pivot %e\n", Ajj);
        return kRetInvalidPivot;
      }

      // update column j
      if (j < n - 1) {
        callAndTime_dgemv('N', M, N, -1.0, &A[j + 1], lda, temp.data(), 1, 1.0,
                          &A[j + 1 + j * lda], 1, DC);
      }

      // add regularization
      double old_pivot = Ajj;
      Ajj = regularizePivot(Ajj, thresh, pivot_sign, A, lda, j, n, uplo, DC);
      regul[j] = fabs(Ajj - old_pivot);
      DC.setMaxReg(regul[j]);

      // save diagonal element
      A[j + lda * j] = Ajj;
      const double coeff = 1.0 / Ajj;

      // scale column j
      if (j < n - 1) {
        callAndTime_dscal(M, coeff, &A[j + 1 + j * lda], 1, DC);
      }
    }
  } else {
    // allocate space for copy of col multiplied by pivots
    std::vector<double> temp(n - 1);

    for (int j = 0; j < n; ++j) {
      const int N = j;
      const int M = n - j - 1;

      // create temporary copy of col j, multiplied by pivots
      for (int i = 0; i < j; ++i) {
        temp[i] = A[i + j * lda] * A[i + i * lda];
      }

      // update diagonal element
      double Ajj =
          pivotCompensatedSum(A[j + lda * j], N, &A[lda * j], 1, A, lda + 1);

      if (isnan(Ajj)) {
        A[j + lda * j] = Ajj;
        printf("\ndense_fact_fiuf: invalid pivot %e\n", Ajj);
        return kRetInvalidPivot;
      }

      // update column j
      if (j < n - 1) {
        callAndTime_dgemv('T', N, M, -1.0, &A[(j + 1) * lda], lda, temp.data(),
                          1, 1.0, &A[j + (j + 1) * lda], lda, DC);
      }

      // add regularization
      double old_pivot = Ajj;
      Ajj = regularizePivot(Ajj, thresh, pivot_sign, A, lda, j, n, uplo, DC);
      regul[j] = fabs(Ajj - old_pivot);
      DC.setMaxReg(regul[j]);

      // save diagonal element
      A[j + lda * j] = Ajj;
      const double coeff = 1.0 / Ajj;

      // scale column j
      if (j < n - 1) {
        callAndTime_dscal(M, coeff, &A[j + (j + 1) * lda], lda, DC);
      }
    }
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
// No pivoting is performed. Regularization is used to deal with small pivots.
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
// ===========================================================================
// In the following, block D represents the current diagonal block, in the block
// of columns that is being factorized. Block R is the rectangular block below
// D. Block P is the square or rectangular block to the left of D. Block Q is
// the rectangular block below P. See the report for a full explanation.
// ===========================================================================

int denseFactF(int n, int k, int nb, double* A, int lda, double* B, int ldb,
               const int* pivot_sign, double thresh, double* regul,
               DataCollector& DC) {
  // ===========================================================================
  // Blocked factorization in full format.
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
    const int jb = std::min(nb, k - j);

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
    std::vector<double> T(j * jb);

    int ldt = jb;
    for (int i = 0; i < j; ++i) {
      callAndTime_dcopy(N, &A[j + i * lda], 1, &T[i * ldt], 1, DC);
      callAndTime_dscal(N, A[i + i * lda], &T[i * ldt], 1, DC);
    }

    // update diagonal block using dgemm_
    callAndTime_dgemm('N', 'T', jb, jb, j, -1.0, P, lda, T.data(), ldt, 1.0, D,
                      lda, DC);

    // factorize diagonal block
    const int* pivot_sign_current = &pivot_sign[j];
    double* regul_current = &regul[j];
    int info = callAndTime_denseFactK('L', N, D, lda, pivot_sign_current,
                                      thresh, regul_current, DC);
    if (info != 0) return info;

    if (j + jb < n) {
      // update block of columns
      callAndTime_dgemm('N', 'T', M, N, K, -1.0, Q, lda, T.data(), ldt, 1.0, R,
                        lda, DC);

      // solve block of columns with L
      callAndTime_dtrsm('R', 'L', 'T', 'U', M, N, 1.0, D, lda, R, lda, DC);

      // solve block of columns with D
      for (int i = 0; i < jb; ++i) {
        const double coeff = 1.0 / A[j + i + (j + i) * lda];
        callAndTime_dscal(M, coeff, &A[j + jb + lda * (j + i)], 1, DC);
      }
    }
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
    std::vector<double> temp_pos((n - k) * pos_pivot);
    std::vector<double> temp_neg((n - k) * neg_pivot);

    const int ldt = n - k;

    // the copies of the columns are multiplied by sqrt(|Ajj|)
    int start_pos = 0;
    int start_neg = 0;
    for (int j = 0; j < k; ++j) {
      double Ajj = A[j + lda * j];
      if (Ajj >= 0.0) {
        Ajj = sqrt(Ajj);
        callAndTime_dcopy(N, &A[k + j * lda], 1, &temp_pos[start_pos * ldt], 1,
                          DC);
        callAndTime_dscal(N, Ajj, &temp_pos[start_pos * ldt], 1, DC);
        ++start_pos;
      } else {
        Ajj = sqrt(-Ajj);
        callAndTime_dcopy(N, &A[k + j * lda], 1, &temp_neg[start_neg * ldt], 1,
                          DC);
        callAndTime_dscal(N, Ajj, &temp_neg[start_neg * ldt], 1, DC);
        ++start_neg;
      }
    }

    // Update schur complement by subtracting contribution of positive columns
    // and adding contribution of negative columns.
    // In this way, I can use dsyrk_ instead of dgemm_ and avoid updating the
    // full square schur complement.

    callAndTime_dsyrk('L', 'N', N, pos_pivot, -1.0, temp_pos.data(), ldt, 1.0,
                      B, ldb, DC);
    callAndTime_dsyrk('L', 'N', N, neg_pivot, 1.0, temp_neg.data(), ldt, 1.0, B,
                      ldb, DC);
  }

  return kRetOk;
}

int denseFactHP(int n, int k, int nb, double* A, double* B,
                const int* pivot_sign, double thresh, double* regul,
                DataCollector& DC) {
  // ===========================================================================
  // Blocked factorization in hybrid-packed format.
  // A should be in lower-blocked-hybrid format. Schur complement is returned
  // in B in lower packed format (not lower-blocked-hybrid).
  // BLAS calls: dgemm_, dtrsm_, dcopy_, dscal_, daxpy_
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
  std::vector<int> diag_start(n_blocks);

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
  std::vector<double> D(nb * nb);

  // buffer for copy of block scaled by pivots
  std::vector<double> T(nb * nb);

  // j is the index of the block column
  for (int j = 0; j < n_blocks; ++j) {
    // jb is the number of columns
    const int jb = std::min(nb, k - nb * j);

    // size of current block could be smaller than diag_size and full_size
    const int this_diag_size = jb * (jb + 1) / 2;
    const int this_full_size = nb * jb;

    // full copy of diagonal block by rows, in D
    int offset = 0;
    for (int Drow = 0; Drow < jb; ++Drow) {
      const int N = Drow + 1;
      callAndTime_dcopy(N, &A[diag_start[j] + offset], 1, &D[Drow * jb], 1, DC);
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
      callAndTime_dcopy(this_full_size, Pk, 1, T.data(), 1, DC);

      // scale temp by pivots
      int pivot_pos = diag_start[k];
      for (int col = 0; col < nb; ++col) {
        callAndTime_dscal(jb, A[pivot_pos], &T[col], nb, DC);
        pivot_pos += col + 2;
      }

      // update diagonal block with dgemm_
      callAndTime_dgemm('T', 'N', jb, jb, nb, -1.0, T.data(), nb, Pk, nb, 1.0,
                        D.data(), jb, DC);

      // update rectangular block
      if (M > 0) {
        const int Qk_pos = Pk_pos + this_full_size;
        const double* Qk = &A[Qk_pos];
        callAndTime_dgemm('T', 'N', jb, M, nb, -1.0, T.data(), nb, Qk, nb, 1.0,
                          R, jb, DC);
      }
    }

    // factorize diagonal block
    double* regul_current = &regul[j * nb];
    const int* pivot_sign_current = &pivot_sign[j * nb];
    int info = callAndTime_denseFactK('U', jb, D.data(), jb, pivot_sign_current,
                                      thresh, regul_current, DC);
    if (info != 0) return info;

    if (M > 0) {
      // solve block of columns with diagonal block
      callAndTime_dtrsm('L', 'U', 'T', 'U', jb, M, 1.0, D.data(), jb, R, jb,
                        DC);

      // scale columns by pivots
      int pivot_pos = 0;
      for (int col = 0; col < jb; ++col) {
        const double coeff = 1.0 / D[pivot_pos];
        callAndTime_dscal(M, coeff, &A[R_pos + col], jb, DC);
        pivot_pos += jb + 1;
      }
    }

    // put D back into packed format
    offset = 0;
    for (int Drow = 0; Drow < jb; ++Drow) {
      const int N = Drow + 1;
      callAndTime_dcopy(N, &D[Drow * jb], 1, &A[diag_start[j] + offset], 1, DC);
      offset += N;
    }
  }

  // compute Schur complement if partial factorization is required
  if (k < n) {
    // number of rows/columns in the Schur complement
    const int ns = n - k;

    // size of last full block (may be smaller than full_size)
    const int ncol_last = k % nb;
    const int last_full_size = ncol_last == 0 ? full_size : ncol_last * nb;

    double beta = 0.0;

    // buffer for full-format of block of columns of Schur complement
    std::vector<double> schur_buf(ns * nb);

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
      const int ncol = std::min(nb, nrow);

      // number of elements in diagonal block
      const int schur_diag_size = ncol * (ncol + 1) / 2;

      beta = 0.0;

      // each block receives contributions from the blocks of the leading part
      // of A
      for (int j = 0; j < n_blocks; ++j) {
        const int jb = std::min(nb, k - nb * j);
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
        callAndTime_dcopy(N, &A[diag_pos], 1, T.data(), 1, DC);

        int pivot_pos = diag_start[j];

        for (int col = 0; col < jb; ++col) {
          callAndTime_dscal(ncol, A[pivot_pos], &T[col], jb, DC);
          pivot_pos += col + 2;
        }

        // update diagonal block using dgemm_
        // printf("T: size %d, access %d\n", nb * nb, jb * ncol);
        // printf("A: size %d, access %d\n", sizeA, diag_pos + ncol * jb);
        // printf("S: size %d, acceess %d\n", ns * nb, ncol * ncol);

        callAndTime_dgemm('T', 'N', ncol, ncol, jb, -1.0, &A[diag_pos], jb,
                          T.data(), jb, beta, schur_buf.data(), ncol, DC);

        // update subdiagonal part
        const int M = nrow - nb;
        if (M > 0) {
          callAndTime_dgemm('T', 'N', ncol, M, jb, -1.0, T.data(), jb,
                            &A[diag_pos + this_full_size], jb, beta,
                            &schur_buf[ncol * ncol], ncol, DC);
        }

        // beta is 0 for the first time (to avoid initializing schur_buf) and
        // 1 for the next calls
        beta = 1.0;
      }

      // schur_buf contains Schur complement in hybrid format (with full
      // diagonal blocks). Sum it in lower-packed format in B.
      for (int buf_row = 0; buf_row < nrow; ++buf_row) {
        const int N = ncol;
        callAndTime_daxpy(N, 1.0, &schur_buf[buf_row * ncol], 1,
                          &B[B_start + buf_row], nrow, DC);
      }
      B_start += nrow * ncol;
    }
  }

  return kRetOk;
}

int denseFactHH(int n, int k, int nb, double* A, double* B,
                const int* pivot_sign, double thresh, double* regul,
                DataCollector& DC) {
  // ===========================================================================
  // Blocked factorization in hybrid-hybrid format.
  // A should be in lower-blocked-hybrid format. Schur complement is returned
  // in B in in lower-block-hybrid format, with full diagonal blocks.
  // BLAS calls: dgemm_, dtrsm_, dcopy_, dscal_
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
  std::vector<int> diag_start(n_blocks);
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
  std::vector<double> D(nb * nb);

  // buffer for copy of block scaled by pivots
  std::vector<double> T(nb * nb);

  // j is the index of the block column
  for (int j = 0; j < n_blocks; ++j) {
    // jb is the number of columns
    const int jb = std::min(nb, k - nb * j);

    // size of current block could be smaller than diag_size and full_size
    const int this_diag_size = jb * (jb + 1) / 2;
    const int this_full_size = nb * jb;

    // full copy of diagonal block by rows, in D
    int offset = 0;
    for (int Drow = 0; Drow < jb; ++Drow) {
      const int N = Drow + 1;
      callAndTime_dcopy(N, &A[diag_start[j] + offset], 1, &D[Drow * jb], 1, DC);
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
      callAndTime_dcopy(this_full_size, Pk, 1, T.data(), 1, DC);

      // scale temp by pivots
      int pivot_pos = diag_start[k];
      for (int col = 0; col < nb; ++col) {
        callAndTime_dscal(jb, A[pivot_pos], &T[col], nb, DC);
        pivot_pos += col + 2;
      }

      // update diagonal block with dgemm_
      callAndTime_dgemm('T', 'N', jb, jb, nb, -1.0, T.data(), nb, Pk, nb, 1.0,
                        D.data(), jb, DC);

      // update rectangular block
      if (M > 0) {
        const int Qk_pos = Pk_pos + this_full_size;
        double* Qk = &A[Qk_pos];

        callAndTime_dgemm('T', 'N', jb, M, nb, -1.0, T.data(), nb, Qk, nb, 1.0,
                          R, jb, DC);
      }
    }

    // factorize diagonal block
    const int* pivot_sign_current = &pivot_sign[j * nb];
    double* regul_current = &regul[j * nb];
    int info = callAndTime_denseFactK('U', jb, D.data(), jb, pivot_sign_current,
                                      thresh, regul_current, DC);
    if (info != 0) return info;

    if (M > 0) {
      // solve block of columns with diagonal block
      callAndTime_dtrsm('L', 'U', 'T', 'U', jb, M, 1.0, D.data(), jb, R, jb,
                        DC);

      // scale columns by pivots
      int pivot_pos = 0;
      for (int col = 0; col < jb; ++col) {
        double coeff = 1.0 / D[pivot_pos];
        callAndTime_dscal(M, coeff, &A[R_pos + col], jb, DC);
        pivot_pos += jb + 1;
      }
    }

    // put D back into packed format
    offset = 0;
    for (int Drow = 0; Drow < jb; ++Drow) {
      const int N = Drow + 1;
      callAndTime_dcopy(N, &D[Drow * jb], 1, &A[diag_start[j] + offset], 1, DC);
      offset += N;
    }
  }

  // compute Schur complement if partial factorization is required
  if (k < n) {
    // number of rows/columns in the Schur complement
    const int ns = n - k;

    // size of last full block (may be smaller than full_size)
    const int ncol_last = k % nb;
    const int last_full_size = ncol_last == 0 ? full_size : ncol_last * nb;

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
      const int ncol = std::min(nb, nrow);

      // each block receives contributions from the blocks of the leading part
      // of A
      for (int j = 0; j < n_blocks; ++j) {
        const int jb = std::min(nb, k - nb * j);
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
        callAndTime_dcopy(N, &A[diag_pos], 1, T.data(), 1, DC);

        int pivot_pos = diag_start[j];
        for (int col = 0; col < jb; ++col) {
          callAndTime_dscal(ncol, A[pivot_pos], &T[col], jb, DC);
          pivot_pos += col + 2;
        }

        // update diagonal block using dgemm_
        callAndTime_dgemm('T', 'N', ncol, ncol, jb, -1.0, T.data(), jb,
                          &A[diag_pos], jb, 1.0, schur_buf, ncol, DC);

        // update subdiagonal part
        const int M = nrow - nb;
        if (M > 0) {
          callAndTime_dgemm('T', 'N', ncol, M, jb, -1.0, T.data(), jb,
                            &A[diag_pos + this_full_size], jb, 1.0,
                            &schur_buf[ncol * ncol], ncol, DC);
        }
      }
      B_start += nrow * ncol;
    }
  }

  return kRetOk;
}

int denseFactL2H(double* A, int nrow, int ncol, int nb, DataCollector& DC) {
  // ===========================================================================
  // Takes a matrix in lower-packed format, with nrow rows.
  // Converts the first ncol columns into lower-blocked-hybrid format, with
  // block size nb.
  // BLAS calls: dcopy_
  // ===========================================================================

#ifdef FINEST_TIMING
  double t0 = GetTime();
#endif

  std::vector<double> buf(nrow * nb);

  int startAtoBuf = 0;
  int startBuftoA = 0;

  for (int k = 0; k <= (ncol - 1) / nb; ++k) {
    // Each block of columns is copied into the buffer, leaving space in
    // between columns, to align rows. This "empty space" contains elements
    // from the previous use of buffer, but they are ignored when copying
    // back.

    // Number of columns in the block. Can be smaller than nb for last block.
    const int block_size = std::min(nb, ncol - k * nb);

    // Number of rows in the block
    const int row_size = nrow - k * nb;

    int startBuf = 0;
    for (int j = 0; j < block_size; ++j) {
      // Copy each column in the block
      const int N = row_size - j;
      callAndTime_dcopy(N, &A[startAtoBuf], 1, &buf[startBuf], 1, DC);
      startAtoBuf += N;

      // +1 to align rows
      startBuf += row_size + 1;
    }

    // Copy columns back into A, row by row.
    // One call of dcopy_ for each row of the block of columns.
    for (int i = 0; i < row_size; ++i) {
      const int N = std::min(i + 1, block_size);
      callAndTime_dcopy(N, &buf[i], row_size, &A[startBuftoA], 1, DC);
      startBuftoA += N;
    }
  }

#ifdef FINEST_TIMING
  DC.sumTime(kTimeDenseFact_convert, GetTime() - t0);
#endif

  return kRetOk;
}