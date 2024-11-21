#include "Auxiliary.h"
#include "CallAndTimeBlas.h"
#include "DataCollector.h"
#include "FactorHiGHSSettings.h"
#include "ReturnValues.h"

// Factorization with "hybrid formats".

int denseFactFH(char format, int n, int k, int nb, double* A, double* B,
                const int* pivot_sign, double thresh, double* regul, int* swaps,
                double* pivot_2x2, DataCollector& DC, int sn) {
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

  // number of rows/columns in the Schur complement
  const int ns = n - k;

  // number of blocks in Schur complement
  const int s_blocks = (ns - 1) / nb + 1;

  // buffer for full-format of block of columns of Schur complement
  std::vector<double> schur_buf;
  if (format == 'P') schur_buf.resize(ns * nb);

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
    std::vector<int> pivot_sign_current(&pivot_sign[j * nb],
                                        &pivot_sign[j * nb] + jb);
    int* swaps_current = &swaps[j * nb];
    double* pivot_2x2_current = &pivot_2x2[j * nb];
    int info = callAndTime_denseFactK('U', jb, D, jb, pivot_sign_current.data(),
                                      thresh, regul_current, swaps_current,
                                      pivot_2x2_current, DC, sn, j);
    if (info != 0) return info;

#ifdef PIVOTING
    // swap columns in R
    applySwaps(swaps_current, M, jb, R, DC);

    // unswap regularization, to keep it with original ordering
    permuteWithSwaps(regul_current, swaps_current, jb, true);
#endif

    // ===========================================================================
    // SOLVE COLUMNS
    // ===========================================================================
    if (M > 0) {
      // solve block R with D
      callAndTime_dtrsm('L', 'U', 'T', 'U', jb, M, 1.0, D, jb, R, jb, DC);

      // make copy of partially solved columns
      callAndTime_dcopy(jb * M, R, 1, T.data(), 1, DC);

      // solve block R with pivots
      int step = 1;
      for (int col = 0; col < jb; col += step) {
        if (pivot_2x2_current[col] == 0.0) {
          // 1x1 pivots
          step = 1;
          const double coeff = 1.0 / D[col + jb * col];
          callAndTime_dscal(M, coeff, &R[col], jb, DC);
        } else {
          // 2x2 pivots
          step = 2;

          // columns affected
          double* c1 = &R[col];
          double* c2 = &R[col + 1];

          // pivot is [d1 offd; offd d2]
          double d1 = D[col + jb * col];
          double d2 = D[col + 1 + jb * (col + 1)];
          double offd = pivot_2x2_current[col];

          // compute coefficients of 2x2 inverse
          const double denom = d1 * d2 - offd * offd;
          const double i_d1 = d2 / denom;
          const double i_d2 = d1 / denom;
          const double i_off = -offd / denom;

          // copy of original col1
          std::vector<double> c1_temp(M);
          callAndTime_dcopy(M, c1, jb, c1_temp.data(), 1, DC);

          // solve col and col+1
          callAndTime_dscal(M, i_d1, c1, jb, DC);
          callAndTime_daxpy(M, i_off, c2, jb, c1, jb, DC);
          callAndTime_dscal(M, i_d2, c2, jb, DC);
          callAndTime_daxpy(M, i_off, c1_temp.data(), 1, c2, jb, DC);
        }
      }

      // ===========================================================================
      // UPDATE FRONTAL
      // ===========================================================================
      int offset{};

      // go through remaining blocks of columns
      for (int jj = j + 1; jj < n_blocks; ++jj) {
        // number of columns in block jj
        const int col_jj = std::min(nb, k - nb * jj);

        // number of rows in block jj
        const int row_jj = n - nb * jj;

        // offset to access T and R
        // offset = (jj - j - 1) * nb * jb;

        const double* P = &T[offset];
        double* Q = &A[diag_start[jj]];
        const double* Rjj = &R[offset];

        callAndTime_dgemm('T', 'N', col_jj, row_jj, jb, -1.0, P, jb, Rjj, jb,
                          1.0, Q, col_jj, DC);

        offset += jb * col_jj;
      }

#ifdef FINE_TIMING
      DC.sumTime(kTimeDenseFact_main, clock.stop());
      clock.start();
#endif

      // ===========================================================================
      // UPDATE SCHUR COMPLEMENT
      // ===========================================================================
      if (k < n) {
        int B_offset{};

        // go through blocks of columns of the Schur complement
        for (int sb = 0; sb < s_blocks; ++sb) {
          // number of rows of the block
          const int nrow = ns - nb * sb;

          // number of columns of the block
          const int ncol = std::min(nb, nrow);

          const double* P = &T[offset];
          double* Q = format == 'P' ? schur_buf.data() : &B[B_offset];
          const double* Rjj = &R[offset];

          // beta is 0 to avoid initializing schur_buf if format=='P'
          double beta = format == 'P' ? 0.0 : 1.0;

          callAndTime_dgemm('T', 'N', ncol, nrow, jb, -1.0, P, jb, Rjj, jb,
                            beta, Q, ncol, DC);

          if (format == 'P') {
            // schur_buf contains Schur complement in hybrid format (with full
            // diagonal blocks). Store it by columns in B (with full diagonal
            // blocks).
            for (int buf_row = 0; buf_row < nrow; ++buf_row) {
              const int N = ncol;
              callAndTime_daxpy(N, 1.0, &schur_buf[buf_row * ncol], 1,
                                &B[B_offset + buf_row], nrow, DC);
            }
          }

          B_offset += nrow * ncol;
          offset += jb * ncol;
        }
      }
#ifdef FINE_TIMING
      DC.sumTime(kTimeDenseFact_schur, clock.stop());
      clock.start();
#endif
    }
  }

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