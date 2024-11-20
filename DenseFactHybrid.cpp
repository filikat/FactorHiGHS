#include "Auxiliary.h"
#include "CallAndTimeBlas.h"
#include "DataCollector.h"
#include "ReturnValues.h"

// Factorization with "hybrid formats".

void applySwaps(const int* swaps, int nrow, int ncol, double* R,
                DataCollector& DC) {
  // apply the column swaps to block R

  for (int i = 0; i < ncol; ++i) {
    // current pivot has been accepted, don't swap
    if (swaps[i] == i) continue;

    // swap col i and col swaps[i]
    callAndTime_dswap(nrow, &R[i], ncol, &R[swaps[i]], ncol, DC);
  }
}

int denseFactFH(char format, int n, int k, int nb, double* A, double* B,
                const int* pivot_sign, double thresh, double* regul, int* swaps,
                DataCollector& DC, int sn) {
  // ===========================================================================
  // Partial blocked factorization
  // Matrix A is in format FH
  // Matrix B is in format FP, if format == 'P'
  //                       FH, if format == 'H'
  // BLAS calls: dcopy, dscal, daxpy, dgemm, dtrsm
  // ===========================================================================

#ifdef FINE_TIMING
  Clock clock;
#endif

  // check input
  if (n < 0 || k < 0 || !A || (k < n && !B)) {
    printf("\ndenseFactH: invalid input\n");
    return kRetInvalidInput;
  }

  // quick return
  if (n == 0) return kRetOk;

  // number of blocks of columns
  const int n_blocks = (k - 1) / nb + 1;

  // start of diagonal blocks
  std::vector<int> diag_start(n_blocks);
  getDiagStart(n, k, nb, n_blocks, diag_start);

  // size of blocks
  const int diag_size = nb * nb;
  const int full_size = nb * nb;

  // buffer for copy of block column
  std::vector<double> T(n * nb);

  // ===========================================================================
  // LOOP OVER BLOCKS
  // ===========================================================================
  for (int j = 0; j < n_blocks; ++j) {
    // j is the index of the block column

    // jb is the number of columns
    const int jb = std::min(nb, k - nb * j);

    // size of current block could be smaller than diag_size and full_size
    const int this_diag_size = jb * jb;
    const int this_full_size = nb * jb;

    // diagonal block j
    double* D = &A[diag_start[j]];

    // number of rows left below block j
    const int M = n - nb * j - jb;

    // block of columns below diagonal block j
    const int R_pos = diag_start[j] + this_diag_size;
    double* R = &A[R_pos];

    // ===========================================================================
    // FACTORIZE DIAGONAL BLOCK
    // ===========================================================================
    double* regul_current = &regul[j * nb];
    const int* pivot_sign_current = &pivot_sign[j * nb];
    int* swaps_current = &swaps[j * nb];
    int info =
        callAndTime_denseFactK('U', jb, D, jb, pivot_sign_current, thresh,
                               regul_current, swaps_current, DC, sn, j);
    if (info != 0) return info;

    // swap columns in R
    applySwaps(swaps_current, M, jb, R, DC);

    // ===========================================================================
    // SOLVE COLUMNS
    // ===========================================================================
    if (M > 0) {
      // solve block R with D
      callAndTime_dtrsm('L', 'U', 'T', 'U', jb, M, 1.0, D, jb, R, jb, DC);

      // make copy of partially solved columns
      callAndTime_dcopy(jb * M, R, 1, T.data(), 1, DC);

      // solve block R with pivots
      int pivot_pos = 0;
      for (int col = 0; col < jb; ++col) {
        const double coeff = 1.0 / D[pivot_pos];
        callAndTime_dscal(M, coeff, &R[col], jb, DC);
        pivot_pos += jb + 1;
      }

      // ===========================================================================
      // UPDATE
      // ===========================================================================
      // go through remaining blocks of columns
      for (int jj = j + 1; jj < n_blocks; ++jj) {
        // number of columns in block jj
        const int col_jj = std::min(nb, k - nb * jj);

        // number of rows in block jj
        const int row_jj = n - nb * jj;

        // offset to access T and R
        const int offset = (jj - j - 1) * nb * jb;

        const double* P = &T[offset];
        double* Q = &A[diag_start[jj]];
        const double* Rjj = &R[offset];

        callAndTime_dgemm('T', 'N', col_jj, row_jj, jb, -1.0, P, jb, Rjj, jb,
                          1.0, Q, col_jj, DC);
      }
    }
  }

#ifdef FINE_TIMING
  DC.sumTime(kTimeDenseFact_main, clock.stop());
  clock.start();
#endif

  // ===========================================================================
  // COMPUTE SCHUR COMPLEMENT
  // ===========================================================================
  if (k < n) {
    // number of rows/columns in the Schur complement
    const int ns = n - k;

    // size of last full block (may be smaller than full_size)
    const int ncol_last = k % nb;
    const int last_full_size = ncol_last == 0 ? full_size : ncol_last * nb;

    double beta = 1.0;

    // buffer for full-format of block of columns of Schur complement
    std::vector<double> temp_buf;
    if (format == 'P') temp_buf.resize(ns * nb);

    // number of blocks in Schur complement
    const int s_blocks = (ns - 1) / nb + 1;

    // index to write into B
    int B_start = 0;

    // pointer to access either temp_buf or B
    double* schur_buf;

    // Go through block of columns of Schur complement.
    // Each block will have a full-format copy made in schur_buf
    for (int sb = 0; sb < s_blocks; ++sb) {
      // select where to write the Schur complement based on format
      if (format == 'P')
        schur_buf = temp_buf.data();
      else
        schur_buf = &B[B_start];

      // number of rows of the block
      const int nrow = ns - nb * sb;

      // number of columns of the block
      const int ncol = std::min(nb, nrow);

      if (format == 'P') beta = 0.0;

      // each block receives contributions from the blocks of the leading part
      // of A
      for (int j = 0; j < n_blocks; ++j) {
        const int jb = std::min(nb, k - nb * j);
        const int this_diag_size = jb * jb;
        const int this_full_size = nb * jb;

        // compute index to access diagonal block in A
        int diag_pos = diag_start[j] + this_diag_size;
        if (j < n_blocks - 1) {
          diag_pos += (n_blocks - j - 2) * full_size + last_full_size;
        }
        diag_pos += sb * this_full_size;

        const double* Pj = &A[diag_pos];
        double* D = schur_buf;
        double* R = &schur_buf[ncol * ncol];

        // create copy of block, multiplied by pivots
        const int N = ncol * jb;
        callAndTime_dcopy(N, Pj, 1, T.data(), 1, DC);

        int pivot_pos = diag_start[j];
        for (int col = 0; col < jb; ++col) {
          callAndTime_dscal(ncol, A[pivot_pos], &T[col], jb, DC);
          pivot_pos += jb + 1;
        }

        // update diagonal block
        callAndTime_dgemm('T', 'N', ncol, ncol, jb, -1.0, Pj, jb, T.data(), jb,
                          beta, D, ncol, DC);

        // update subdiagonal part
        const int M = nrow - nb;
        if (M > 0) {
          const double* Qj = &Pj[this_full_size];
          callAndTime_dgemm('T', 'N', ncol, M, jb, -1.0, T.data(), jb, Qj, jb,
                            beta, R, ncol, DC);
        }

        // beta is 0 for the first time (to avoid initializing schur_buf) and
        // 1 for the next calls (if format=='P')
        beta = 1.0;
      }

      if (format == 'P') {
        // schur_buf contains Schur complement in hybrid format (with full
        // diagonal blocks). Store it by columns in B (with full diagonal
        // blocks).
        for (int buf_row = 0; buf_row < nrow; ++buf_row) {
          const int N = ncol;
          callAndTime_daxpy(N, 1.0, &schur_buf[buf_row * ncol], 1,
                            &B[B_start + buf_row], nrow, DC);
        }
      }

      B_start += nrow * ncol;
    }
  }

#ifdef FINE_TIMING
  DC.sumTime(kTimeDenseFact_schur, clock.stop());
#endif

  return kRetOk;
}

int denseFactFP2FH(double* A, int nrow, int ncol, int nb, DataCollector& DC) {
  // ===========================================================================
  // Packed to Hybrid conversion
  // Matrix A on  input is in format FP
  // Matrix A on output is in format FH
  // BLAS calls: dcopy
  // ===========================================================================

#ifdef FINE_TIMING
  Clock clock;
#endif

  std::vector<double> buf(nrow * nb);

  int startAtoBuf = 0;
  int startBuftoA = 0;

  for (int k = 0; k <= (ncol - 1) / nb; ++k) {
    // Number of columns in the block. Can be smaller than nb for last block.
    const int block_size = std::min(nb, ncol - k * nb);

    // Number of rows in the block
    const int row_size = nrow - k * nb;

    // Copy block into buf
    callAndTime_dcopy(row_size * block_size, &A[startAtoBuf], 1, buf.data(), 1,
                      DC);
    startAtoBuf += row_size * block_size;

    // Copy columns back into A, row by row.
    // One call of dcopy_ for each row of the block of columns.
    for (int i = 0; i < row_size; ++i) {
      const int N = block_size;
      callAndTime_dcopy(N, &buf[i], row_size, &A[startBuftoA], 1, DC);
      startBuftoA += N;
    }
  }

#ifdef FINE_TIMING
  DC.sumTime(kTimeDenseFact_convert, clock.stop());
#endif

  return kRetOk;
}
