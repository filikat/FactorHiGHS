#include "DenseFact.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <vector>

#include "../ProtoIPM/Regularization.h"
#include "Auxiliary.h"
#include "CallAndTimeBlas.h"
#include "DataCollector.h"
#include "ReturnValues.h"
#include "Timing.h"
#include "parallel/HighsParallel.h"

// #define DEBUG

#define PARALLEL_NODE

/*
  Names:
  denseFact:
  - K : factorization kernel for diagonal blocks
  - F : blocked factorization in full format
  - FP: blocked factorization in format FP
  - H : blocked factorization in "hybrid formats"

  Formats used:
  - F : Full format
  - P : lower Packed format
  - FP: lower Packed format with Full diagonal blocks
  - H : lower-blocked-Hybrid format
  - FH: lower-blocked-Hybrid format with Full diagonal blocks

  F, P do not use blocks. FP, H, FH use blocks.
  Blocks are always blocks of columns.
  F, P store by columns.
  FP stores by columns within the blocks. H, FH store by rows within the blocks.
  See report for more details.
*/

double regularizePivot(double pivot, double thresh, const int* sign,
                       const double* A, int lda, int j, int n, char uplo,
                       DataCollector& DC, int sn, int bl) {
  // add static regularization
  if (sign[j] == 1)
    pivot += kDualStaticRegularization;
  else
    pivot -= kPrimalStaticRegularization;

  double s = (double)sign[j];
  double old_pivot = pivot;

  double spivot = s * pivot;
  double K = 1e12;

  bool adjust = false;
  bool modified_pivot = false;

  if (spivot <= thresh && spivot >= -thresh) {
    // small pivot, lift to thresh
    pivot = s * thresh;
    adjust = true;
    modified_pivot = true;
#ifdef DEBUG
    printf("%2d, %2d, %2d: small pivot %e, with sign %d, set to %e\n", sn, bl,
           j, old_pivot, sign[j], pivot);
#endif

  } else if (spivot < -thresh && spivot >= -thresh * K) {
    // wrong sign, lift more
    pivot = s * thresh * 10;
    adjust = true;
    modified_pivot = true;
#ifdef DEBUG
    printf("%2d, %2d, %2d: wrong pivot %e, with sign %d, set to %e\n", sn, bl,
           j, old_pivot, sign[j], pivot);
#endif

  } else if (spivot < -thresh * K) {
    // pivot is completely lost
    pivot = s * 1e100;
    modified_pivot = true;
#ifdef DEBUG
    printf("%2d, %2d, %2d: disaster pivot %e, with sign %d, set to %e\n", sn,
           bl, j, old_pivot, sign[j], pivot);
#endif
  }

  if (adjust) {
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
      double sk = sign[k];
      if (s * sk < 0) continue;

      double temp = (dk - sk * thresh);
      temp = (bk * bk) / temp;

      if (s > 0)
        required_pivot = std::max(required_pivot, temp);
      else
        required_pivot = std::min(required_pivot, temp);
    }

    if (required_pivot != pivot) {
      modified_pivot = true;

      if (s > 0)
        pivot = std::max(pivot, required_pivot);
      else
        pivot = std::min(pivot, required_pivot);

#ifdef DEBUG
      printf("\t%2d, %2d, %2d: adjust %e to %e\n", sn, bl, j, old_pivot, pivot);
#endif
    }
  }

  if (modified_pivot) DC.sumRegPiv();

  return pivot;
}

int denseFactK(char uplo, int n, double* A, int lda, const int* pivot_sign,
               double thresh, double* regul, DataCollector& DC, int sn,
               int bl) {
  // ===========================================================================
  // Factorization kernel
  // Matrix A is in format F
  // BLAS calls: dscal, dcopy, daxpy
  // ===========================================================================

  // check input
  if (n < 0 || !A || lda < n) {
    printf("\ndenseFactK: invalid input\n");
    return kRetInvalidInput;
  }

  // quick return
  if (n == 0) return kRetOk;

  // main operations
  if (uplo == 'L') {
    // allocate space for copy of col
    std::vector<double> temp(n - 1);

    for (int j = 0; j < n; ++j) {
      // diagonal element
      double Ajj = A[j + lda * j];

      if (isnan(Ajj)) {
        printf("\ndenseFactK: invalid pivot %e\n", Ajj);
        return kRetInvalidPivot;
      }

      // add regularization
      double old_pivot = Ajj;
      Ajj = regularizePivot(Ajj, thresh, pivot_sign, A, lda, j, n, uplo, DC, sn,
                            bl);
      regul[j] = std::abs(Ajj - old_pivot);
      DC.setMaxReg(regul[j]);

      // save diagonal element
      A[j + lda * j] = Ajj;

      const int M = n - j - 1;
      if (M > 0) {
        // make copy of column
        callAndTime_dcopy(M, &A[j + 1 + j * lda], 1, temp.data(), 1, DC);

        // scale column j
        callAndTime_dscal(M, 1.0 / Ajj, &A[j + 1 + j * lda], 1, DC);

        // update rest of the matrix
        for (int i = 0; i < M; ++i) {
          callAndTime_daxpy(M - i, -temp[i], &A[j + 1 + i + j * lda], 1,
                            &A[j + 1 + i + (j + 1 + i) * lda], 1, DC);
        }
      }
    }
  } else {
    // allocate space for copy of col
    std::vector<double> temp(n - 1);

    for (int j = 0; j < n; ++j) {
      // diagonal element
      double Ajj = A[j + lda * j];

      if (isnan(Ajj)) {
        A[j + lda * j] = Ajj;
        printf("\ndenseFactK: invalid pivot %e\n", Ajj);
        return kRetInvalidPivot;
      }

      // add regularization
      double old_pivot = Ajj;
      Ajj = regularizePivot(Ajj, thresh, pivot_sign, A, lda, j, n, uplo, DC, sn,
                            bl);
      regul[j] = std::abs(Ajj - old_pivot);
      DC.setMaxReg(regul[j]);

      // save diagonal element
      A[j + lda * j] = Ajj;
      const double coeff = 1.0 / Ajj;

      const int M = n - j - 1;
      if (M > 0) {
        // make copy of row
        callAndTime_dcopy(M, &A[j + (j + 1) * lda], lda, temp.data(), 1, DC);

        // scale row j
        callAndTime_dscal(M, coeff, &A[j + (j + 1) * lda], lda, DC);

        // update rest of the matrix
        for (int i = 0; i < M; ++i) {
          callAndTime_daxpy(M - i, -temp[i], &A[j + (j + 1 + i) * lda], lda,
                            &A[j + 1 + i + (j + 1 + i) * lda], lda, DC);
        }
      }
    }
  }

  return kRetOk;
}

int denseFactF(int n, int k, int nb, double* A, int lda, double* B, int ldb,
               const int* pivot_sign, double thresh, double* regul,
               DataCollector& DC, int sn) {
  // ===========================================================================
  // Partial blocked factorization
  // Matrix A is in format F
  // Matrix B is in format F
  // BLAS calls: dcopy, dscal, dsyrk, dgemm, dtrsm
  // ===========================================================================

#ifdef FINE_TIMING
  Clock clock;
#endif

  // check input
  if (n < 0 || k < 0 || !A || lda < n || (k < n && (!B || ldb < n - k))) {
    printf("\ndenseFactF: invalid input\n");
    return kRetInvalidInput;
  }

  // quick return
  if (n == 0) return kRetOk;

  // create temporary copy of block of rows, multiplied by pivots
  std::vector<double> T(k * nb);

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

    int ldt = jb;
    for (int i = 0; i < j; ++i) {
      callAndTime_dcopy(N, &P[i * lda], 1, &T[i * ldt], 1, DC);
      callAndTime_dscal(N, A[i + i * lda], &T[i * ldt], 1, DC);
    }

    // update diagonal block using dgemm_
    callAndTime_dgemm('N', 'T', jb, jb, j, -1.0, P, lda, T.data(), ldt, 1.0, D,
                      lda, DC);

    // factorize diagonal block
    const int* pivot_sign_current = &pivot_sign[j];
    double* regul_current = &regul[j];
    int bl = j / nb;
    int info = callAndTime_denseFactK('L', N, D, lda, pivot_sign_current,
                                      thresh, regul_current, DC, sn, bl);
    if (info != 0) return info;

    if (j + jb < n) {
      // update block of columns
      callAndTime_dgemm('N', 'T', M, N, K, -1.0, Q, lda, T.data(), ldt, 1.0, R,
                        lda, DC);

      // solve block of columns with L
      callAndTime_dtrsm('R', 'L', 'T', 'U', M, N, 1.0, D, lda, R, lda, DC);

      // solve block of columns with D
      for (int i = 0; i < jb; ++i) {
        const double coeff = 1.0 / D[i + i * lda];
        callAndTime_dscal(M, coeff, &R[lda * i], 1, DC);
      }
    }
  }

#ifdef FINE_TIMING
  DC.sumTime(kTimeDenseFact_main, clock.stop());
  clock.start();
#endif

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

#ifdef FINE_TIMING
  DC.sumTime(kTimeDenseFact_schur, clock.stop());
#endif

  return kRetOk;
}

int denseFactFP(int n, int k, int nb, double* A, double* B,
                const int* pivot_sign, double thresh, double* regul,
                DataCollector& DC, int sn) {
  // ===========================================================================
  // Partial blocked factorization
  // Matrix A is in format FP
  // Matrix B is in format FP
  // BLAS calls: dcopy, dscal, dgemm, dtrsm
  // ===========================================================================

#ifdef FINE_TIMING
  Clock clock;
#endif

  // check input
  if (n < 0 || k < 0 || !A || (k < n && !B)) {
    printf("\ndenseFactFP: invalid input\n");
    return kRetInvalidInput;
  }

  // quick return
  if (n == 0) return kRetOk;

  // number of blocks of columns
  const int n_blocks = (k - 1) / nb + 1;

  // start of diagonal blocks
  std::vector<int> diag_start(n_blocks);
  getDiagStart(n, k, nb, n_blocks, diag_start);

  // buffer for copy of block scaled by pivots
  std::vector<double> T(nb * nb);

  // j is the index of the block column
  for (int j = 0; j < n_blocks; ++j) {
    // jb is the number of columns
    const int jb = std::min(nb, k - nb * j);

    // number of rows left below block j
    const int M = n - nb * j - jb;

    // diagonal block
    double* D = &A[diag_start[j]];

    // block of columns below diagonal block j
    double* R = &A[diag_start[j] + jb];

    // leading dimensions to access arrays
    const int ldD = n - j * nb;
    const int ldR = ldD;

    // update diagonal block and block of columns
    for (int k = 0; k < j; ++k) {
      // starting position of block P
      int Pk_pos = diag_start[k] + nb * (j - k);
      const double* Pk = &A[Pk_pos];

      // leading dimensions
      const int ldP = n - k * nb;
      const int ldQ = ldP;
      const int ldT = jb;

      // copy block Pk into temp and scale by pivots
      const double* Dk = &A[diag_start[k]];
      for (int col = 0; col < nb; ++col) {
        callAndTime_dcopy(jb, &Pk[col * ldP], 1, &T[col * ldT], 1, DC);
        callAndTime_dscal(jb, Dk[col + col * ldP], &T[col * ldT], 1, DC);
      }

// update diagonal block
#ifdef PARALLEL_NODE
      GemmCaller gemm_diag('N', 'T', jb, jb, nb, -1.0, T.data(), ldT, Pk, ldP,
                           1.0, D, ldD, DC);
      highs::parallel::spawn([&]() { gemm_diag.run(); });
#else
      callAndTime_dgemm('N', 'T', jb, jb, nb, -1.0, T.data(), ldT, Pk, ldP, 1.0,
                        D, ldD, DC);
#endif

      // update rectangular block
      if (M > 0) {
        const int Qk_pos = Pk_pos + jb;
        const double* Qk = &A[Qk_pos];

#ifdef PARALLEL_NODE
        // number of vertical blocks in Q and R
        const int qr_blocks = (M - 1) / nb + 1;

        if (qr_blocks > 3) {
          // execute gemm in parallel

          // vector of GemmCallers, to guarantee that the object exists when the
          // scheduler calls it.
          std::vector<GemmCaller> gemm_callers;
          gemm_callers.reserve(qr_blocks);

          for (int jj = 0; jj < qr_blocks; ++jj) {
            const int jjb = std::min(nb, M - jj * nb);
            const double* Qjj = &Qk[nb * jj];
            double* Rjj = &R[nb * jj];

            gemm_callers.push_back(GemmCaller('N', 'T', jjb, jb, nb, -1.0, Qjj,
                                              ldQ, T.data(), ldT, 1.0, Rjj, ldR,
                                              DC));

            highs::parallel::spawn([&, jj]() { gemm_callers[jj].run(); });
          }

          for (int jj = 0; jj < qr_blocks; ++jj) highs::parallel::sync();

        } else {
          // execute gemm in serial
          callAndTime_dgemm('N', 'T', M, jb, nb, -1.0, Qk, ldQ, T.data(), ldT,
                            1.0, R, ldR, DC);
        }

#else
        callAndTime_dgemm('N', 'T', M, jb, nb, -1.0, Qk, ldQ, T.data(), ldT,
                          1.0, R, ldR, DC);
#endif
      }

#ifdef PARALLEL_NODE
      // sync diagonal block
      highs::parallel::sync();
#endif
    }

    // factorize diagonal block
    double* regul_current = &regul[j * nb];
    const int* pivot_sign_current = &pivot_sign[j * nb];
    int info = callAndTime_denseFactK('L', jb, D, ldD, pivot_sign_current,
                                      thresh, regul_current, DC, sn, j);
    if (info != 0) return info;

    // solve block of columns with diagonal block
    if (M > 0) {
#ifdef PARALLEL_NODE
      // number of vertical blocks in Q and R
      const int r_blocks = (M - 1) / nb + 1;

      if (r_blocks > 3) {
        // execute trsm in parallel

        // vector of TrsmCallers, to guarantee that the object exists when the
        // scheduler calls it.
        std::vector<TrsmCaller> trsm_callers;
        trsm_callers.reserve(r_blocks);

        for (int jj = 0; jj < r_blocks; ++jj) {
          const int jjb = std::min(nb, M - jj * nb);
          double* Rjj = &R[nb * jj];

          trsm_callers.push_back(TrsmCaller('R', 'L', 'T', 'U', jjb, jb, 1.0, D,
                                            ldD, Rjj, ldR, DC));

          highs::parallel::spawn([&, jj]() { trsm_callers[jj].run(); });
        }

        for (int jj = 0; jj < r_blocks; ++jj) highs::parallel::sync();

      } else {
        // execute trsm in serial
        callAndTime_dtrsm('R', 'L', 'T', 'U', M, jb, 1.0, D, ldD, R, ldR, DC);
      }
#else
      callAndTime_dtrsm('R', 'L', 'T', 'U', M, jb, 1.0, D, ldD, R, ldR, DC);
#endif

      // scale columns by pivots
      for (int col = 0; col < jb; ++col) {
        const double coeff = 1.0 / D[col + col * ldD];
        callAndTime_dscal(M, coeff, &R[col * ldR], 1, DC);
      }
    }
  }

#ifdef FINE_TIMING
  DC.sumTime(kTimeDenseFact_main, clock.stop());
  clock.start();
#endif

  // compute Schur complement if partial factorization is required
  if (k < n) {
    // number of rows/columns in the Schur complement
    const int ns = n - k;

    // size of last full block
    const int ncol_last = (k % nb == 0 ? nb : k % nb);

    // number of blocks in Schur complement
    const int s_blocks = (ns - 1) / nb + 1;

    int B_start = 0;

    // Go through block of columns of Schur complement
    for (int sb = 0; sb < s_blocks; ++sb) {
      // number of rows of the block
      const int nrow = ns - nb * sb;

      // number of columns of the block
      const int ncol = std::min(nb, nrow);

      double* D = &B[B_start];
      double* R = &B[B_start + ncol];
      const int ldD = nrow;
      const int ldR = ldD;

      // each block receives contributions from the blocks of the leading part
      // of A
      for (int j = 0; j < n_blocks; ++j) {
        const int jb = std::min(nb, k - nb * j);

        // compute index to access block Pj
        const int Pj_pos =
            diag_start[j] + (n_blocks - j - 1) * jb + ncol_last + sb * nb;
        const double* Pj = &A[Pj_pos];
        const int ldP = n - j * nb;
        const int ldT = ncol;

        // copy block Pj into temp and scale by pivots
        const double* Dj = &A[diag_start[j]];
        for (int col = 0; col < jb; ++col) {
          callAndTime_dcopy(ncol, &Pj[col * ldP], 1, &T[col * ldT], 1, DC);
          callAndTime_dscal(ncol, Dj[col + col * ldP], &T[col * ldT], 1, DC);
        }

        const double* Qj = &A[Pj_pos + ncol];
        const int ldQ = ldP;

// update diagonal block
#ifdef PARALLEL_NODE
        GemmCaller gemm_diag('N', 'T', ncol, ncol, jb, -1.0, Pj, ldP, T.data(),
                             ldT, 1.0, D, ldD, DC);
        highs::parallel::spawn([&]() { gemm_diag.run(); });
#else
        callAndTime_dgemm('N', 'T', ncol, ncol, jb, -1.0, Pj, ldP, T.data(),
                          ldT, 1.0, D, ldD, DC);
#endif

        // update subdiagonal part
        const int M = nrow - ncol;
        if (M > 0) {
#ifdef PARALLEL_NODE
          // number of vertical blocks in Q and R
          const int qr_blocks = (M - 1) / nb + 1;

          if (qr_blocks > 3) {
            // execute gemm in parallel

            // vector of GemmCallers, to guarantee that the object exists when
            // the scheduler calls it.
            std::vector<GemmCaller> gemm_callers;
            gemm_callers.reserve(qr_blocks);

            for (int jj = 0; jj < qr_blocks; ++jj) {
              const int jjb = std::min(nb, M - jj * nb);
              const double* Qjj = &Qj[nb * jj];
              double* Rjj = &R[nb * jj];

              gemm_callers.push_back(GemmCaller('N', 'T', jjb, ncol, jb, -1.0,
                                                Qjj, ldQ, T.data(), ldT, 1.0,
                                                Rjj, ldR, DC));

              highs::parallel::spawn([&, jj]() { gemm_callers[jj].run(); });
            }

            for (int jj = 0; jj < qr_blocks; ++jj) highs::parallel::sync();

          } else {
            // execute gemm in serial
            callAndTime_dgemm('N', 'T', M, ncol, jb, -1.0, Qj, ldQ, T.data(),
                              ldT, 1.0, R, ldR, DC);
          }
#else
          callAndTime_dgemm('N', 'T', M, ncol, jb, -1.0, Qj, ldQ, T.data(), ldT,
                            1.0, R, ldR, DC);
#endif
        }

#ifdef PARALLEL_NODE
        // sync diagonal block
        highs::parallel::sync();
#endif
      }

      B_start += nrow * ncol;
    }
  }

#ifdef FINE_TIMING
  DC.sumTime(kTimeDenseFact_schur, clock.stop());
#endif

  return kRetOk;
}

int denseFactFH(char format, int n, int k, int nb, double* A, double* B,
                const int* pivot_sign, double thresh, double* regul,
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

  // buffer for copy of block scaled by pivots
  std::vector<double> T(nb * nb);

  // j is the index of the block column
  for (int j = 0; j < n_blocks; ++j) {
    // jb is the number of columns
    const int jb = std::min(nb, k - nb * j);

    // size of current block could be smaller than diag_size and full_size
    const int this_diag_size = jb * jb;
    const int this_full_size = nb * jb;

    double* D = &A[diag_start[j]];

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
        pivot_pos += nb + 1;
      }

      // update diagonal block
#ifdef PARALLEL_NODE
      GemmCaller gemm_diag('T', 'N', jb, jb, nb, -1.0, T.data(), nb, Pk, nb,
                           1.0, D, jb, DC);
      highs::parallel::spawn([&]() { gemm_diag.run(); });
#else
      callAndTime_dgemm('T', 'N', jb, jb, nb, -1.0, T.data(), nb, Pk, nb, 1.0,
                        D, jb, DC);
#endif

      // update rectangular block
      if (M > 0) {
        const int Qk_pos = Pk_pos + this_full_size;
        const double* Qk = &A[Qk_pos];

#ifdef PARALLEL_NODE
        // number of vertical blocks in Q and R
        const int qr_blocks = (M - 1) / nb + 1;

        if (qr_blocks > 3) {
          // execute gemm in parallel

          // vector of GemmCallers, to guarantee that the object exists when the
          // scheduler calls it.
          std::vector<GemmCaller> gemm_callers;
          gemm_callers.reserve(qr_blocks);

          for (int jj = 0; jj < qr_blocks; ++jj) {
            const int jjb = std::min(nb, M - jj * nb);
            const double* Qjj = &Qk[nb * nb * jj];
            double* Rjj = &R[nb * jb * jj];

            gemm_callers.push_back(GemmCaller('T', 'N', jb, jjb, nb, -1.0,
                                              T.data(), nb, Qjj, nb, 1.0, Rjj,
                                              jb, DC));

            highs::parallel::spawn([&, jj]() { gemm_callers[jj].run(); });
          }

          for (int jj = 0; jj < qr_blocks; ++jj) highs::parallel::sync();

        } else {
          // execute gemm in serial
          callAndTime_dgemm('T', 'N', jb, M, nb, -1.0, T.data(), nb, Qk, nb,
                            1.0, R, jb, DC);
        }
#else
        callAndTime_dgemm('T', 'N', jb, M, nb, -1.0, T.data(), nb, Qk, nb, 1.0,
                          R, jb, DC);
#endif
      }

#ifdef PARALLEL_NODE
      // sync diagonal block
      highs::parallel::sync();
#endif
    }

    // factorize diagonal block
    double* regul_current = &regul[j * nb];
    const int* pivot_sign_current = &pivot_sign[j * nb];
    int info = callAndTime_denseFactK('U', jb, D, jb, pivot_sign_current,
                                      thresh, regul_current, DC, sn, j);
    if (info != 0) return info;

    if (M > 0) {
#ifdef PARALLEL_NODE
      // number of vertical blocks in R
      const int r_blocks = (M - 1) / nb + 1;

      if (r_blocks > 3) {
        // execute trsm in parallel

        // vector of TrsmCallers, to guarantee that the object exists when the
        // scheduler calls it.
        std::vector<TrsmCaller> trsm_callers;
        trsm_callers.reserve(r_blocks);

        for (int jj = 0; jj < r_blocks; ++jj) {
          const int jjb = std::min(nb, M - jj * nb);
          double* Rjj = &R[nb * jb * jj];

          trsm_callers.push_back(
              TrsmCaller('L', 'U', 'T', 'U', jb, jjb, 1.0, D, jb, Rjj, jb, DC));

          highs::parallel::spawn([&, jj]() { trsm_callers[jj].run(); });
        }

        for (int jj = 0; jj < r_blocks; ++jj) highs::parallel::sync();

      } else {
        // execute trsm in serial
        callAndTime_dtrsm('L', 'U', 'T', 'U', jb, M, 1.0, D, jb, R, jb, DC);
      }
#else
      callAndTime_dtrsm('L', 'U', 'T', 'U', jb, M, 1.0, D, jb, R, jb, DC);
#endif

      // scale columns by pivots
      int pivot_pos = 0;
      for (int col = 0; col < jb; ++col) {
        const double coeff = 1.0 / D[pivot_pos];
        callAndTime_dscal(M, coeff, &A[R_pos + col], jb, DC);
        pivot_pos += jb + 1;
      }
    }
  }

#ifdef FINE_TIMING
  DC.sumTime(kTimeDenseFact_main, clock.stop());
  clock.start();
#endif

  // compute Schur complement if partial factorization is required
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
#ifdef PARALLEL_NODE
        GemmCaller gemm_diag('T', 'N', ncol, ncol, jb, -1.0, Pj, jb, T.data(),
                             jb, beta, D, ncol, DC);
        highs::parallel::spawn([&]() { gemm_diag.run(); });
#else
        callAndTime_dgemm('T', 'N', ncol, ncol, jb, -1.0, Pj, jb, T.data(), jb,
                          beta, D, ncol, DC);
#endif

        // update subdiagonal part
        const int M = nrow - nb;
        if (M > 0) {
          const double* Qj = &Pj[this_full_size];

#ifdef PARALLEL_NODE
          // number of vertical blocks in Q and R
          const int qr_blocks = (M - 1) / nb + 1;

          if (qr_blocks > 3) {
            // execute gemm in parallel

            // vector of GemmCallers, to guarantee that the object exists when
            // the scheduler calls it.
            std::vector<GemmCaller> gemm_callers;
            gemm_callers.reserve(qr_blocks);

            for (int jj = 0; jj < qr_blocks; ++jj) {
              const int jjb = std::min(nb, M - jj * nb);
              const double* Qjj = &Qj[nb * jb * jj];
              double* Rjj = &R[nb * ncol * jj];

              gemm_callers.push_back(GemmCaller('T', 'N', ncol, jjb, jb, -1.0,
                                                T.data(), jb, Qjj, jb, beta,
                                                Rjj, ncol, DC));

              highs::parallel::spawn([&, jj]() { gemm_callers[jj].run(); });
            }

            for (int jj = 0; jj < qr_blocks; ++jj) highs::parallel::sync();

          } else {
            // execute gemm in serial
            callAndTime_dgemm('T', 'N', ncol, M, jb, -1.0, T.data(), jb, Qj, jb,
                              beta, R, ncol, DC);
          }
#else
          callAndTime_dgemm('T', 'N', ncol, M, jb, -1.0, T.data(), jb, Qj, jb,
                            beta, R, ncol, DC);
#endif
        }

#ifdef PARALLEL_NODE
        // sync diagonal block
        highs::parallel::sync();
#endif

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

GemmCaller::GemmCaller(char transa, char transb, int m, int n, int k,
                       double alpha, const double* A, int lda, const double* B,
                       int ldb, double beta, double* C, int ldc,
                       DataCollector& DC)
    : transa_{transa},
      transb_{transb},
      m_{m},
      n_{n},
      k_{k},
      alpha_{alpha},
      A_{A},
      lda_{lda},
      B_{B},
      ldb_{ldb},
      beta_{beta},
      C_{C},
      ldc_{ldc},
      DC_{DC} {}

void GemmCaller::run() {
  callAndTime_dgemm(transa_, transb_, m_, n_, k_, alpha_, A_, lda_, B_, ldb_,
                    beta_, C_, ldc_, DC_);
}

TrsmCaller::TrsmCaller(char side, char uplo, char trans, char diag, int m,
                       int n, double alpha, const double* A, int lda, double* B,
                       int ldb, DataCollector& DC)
    : side_{side},
      uplo_{uplo},
      trans_{trans},
      diag_{diag},
      m_{m},
      n_{n},
      alpha_{alpha},
      A_{A},
      lda_{lda},
      B_{B},
      ldb_{ldb},
      DC_{DC} {}

void TrsmCaller::run() {
  callAndTime_dtrsm(side_, uplo_, trans_, diag_, m_, n_, alpha_, A_, lda_, B_,
                    ldb_, DC_);
}