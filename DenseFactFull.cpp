#include "Auxiliary.h"
#include "CallAndTimeBlas.h"
#include "DataCollector.h"
#include "DenseFact.h"
#include "ReturnValues.h"

// Factorization with "full formats".
// Kept only for reference, "hybrid formats" should be preferred.

int denseFactF(int n, int k, int nb, double* A, int lda, double* B, int ldb,
               const int* pivot_sign, double thresh, double* regul, int sn) {
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
      callAndTime_dcopy(N, &P[i * lda], 1, &T[i * ldt], 1);
      callAndTime_dscal(N, A[i + i * lda], &T[i * ldt], 1);
    }

    // update diagonal block using dgemm_
    callAndTime_dgemm('N', 'T', jb, jb, j, -1.0, P, lda, T.data(), ldt, 1.0, D,
                      lda);

    // factorize diagonal block
    std::vector<int> pivot_sign_current(&pivot_sign[j], &pivot_sign[j] + jb);
    double* regul_current = &regul[j];
    int bl = j / nb;
    int info = denseFactK('L', N, D, lda, pivot_sign_current.data(), thresh,
                          regul_current, nullptr, nullptr, sn, bl);
    if (info != 0) return info;

    if (j + jb < n) {
      // update block of columns
      callAndTime_dgemm('N', 'T', M, N, K, -1.0, Q, lda, T.data(), ldt, 1.0, R,
                        lda);

      // solve block of columns with L
      callAndTime_dtrsm('R', 'L', 'T', 'U', M, N, 1.0, D, lda, R, lda);

      // solve block of columns with D
      for (int i = 0; i < jb; ++i) {
        const double coeff = 1.0 / D[i + i * lda];
        callAndTime_dscal(M, coeff, &R[lda * i], 1);
      }
    }
  }

#ifdef FINE_TIMING
  DataCollector::get()->sumTime(kTimeDenseFact_main, clock.stop());
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
        callAndTime_dcopy(N, &A[k + j * lda], 1, &temp_pos[start_pos * ldt], 1);
        callAndTime_dscal(N, Ajj, &temp_pos[start_pos * ldt], 1);
        ++start_pos;
      } else {
        Ajj = sqrt(-Ajj);
        callAndTime_dcopy(N, &A[k + j * lda], 1, &temp_neg[start_neg * ldt], 1);
        callAndTime_dscal(N, Ajj, &temp_neg[start_neg * ldt], 1);
        ++start_neg;
      }
    }

    // Update schur complement by subtracting contribution of positive columns
    // and adding contribution of negative columns.
    // In this way, I can use dsyrk_ instead of dgemm_ and avoid updating the
    // full square schur complement.

    callAndTime_dsyrk('L', 'N', N, pos_pivot, -1.0, temp_pos.data(), ldt, 1.0,
                      B, ldb);
    callAndTime_dsyrk('L', 'N', N, neg_pivot, 1.0, temp_neg.data(), ldt, 1.0, B,
                      ldb);
  }

#ifdef FINE_TIMING
  DataCollector::get()->sumTime(kTimeDenseFact_schur, clock.stop());
#endif

  return kRetOk;
}

int denseFactFP(int n, int k, int nb, double* A, double* B,
                const int* pivot_sign, double thresh, double* regul, int sn) {
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
        callAndTime_dcopy(jb, &Pk[col * ldP], 1, &T[col * ldT], 1);
        callAndTime_dscal(jb, Dk[col + col * ldP], &T[col * ldT], 1);
      }

      // update diagonal block
      callAndTime_dgemm('N', 'T', jb, jb, nb, -1.0, T.data(), ldT, Pk, ldP, 1.0,
                        D, ldD);

      // update rectangular block
      if (M > 0) {
        const int Qk_pos = Pk_pos + jb;
        const double* Qk = &A[Qk_pos];
        callAndTime_dgemm('N', 'T', M, jb, nb, -1.0, Qk, ldQ, T.data(), ldT,
                          1.0, R, ldR);
      }
    }

    // factorize diagonal block
    double* regul_current = &regul[j * nb];
    std::vector<int> pivot_sign_current(&pivot_sign[j * nb],
                                        &pivot_sign[j * nb] + jb);
    int info = denseFactK('L', jb, D, ldD, pivot_sign_current.data(), thresh,
                          regul_current, nullptr, nullptr, sn, j);
    if (info != 0) return info;

    // solve block of columns with diagonal block
    if (M > 0) {
      callAndTime_dtrsm('R', 'L', 'T', 'U', M, jb, 1.0, D, ldD, R, ldR);

      // scale columns by pivots
      for (int col = 0; col < jb; ++col) {
        const double coeff = 1.0 / D[col + col * ldD];
        callAndTime_dscal(M, coeff, &R[col * ldR], 1);
      }
    }
  }

#ifdef FINE_TIMING
  DataCollector::get()->sumTime(kTimeDenseFact_main, clock.stop());
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
          callAndTime_dcopy(ncol, &Pj[col * ldP], 1, &T[col * ldT], 1);
          callAndTime_dscal(ncol, Dj[col + col * ldP], &T[col * ldT], 1);
        }

        const double* Qj = &A[Pj_pos + ncol];
        const int ldQ = ldP;

        // update diagonal block
        callAndTime_dgemm('N', 'T', ncol, ncol, jb, -1.0, Pj, ldP, T.data(),
                          ldT, 1.0, D, ldD);

        // update subdiagonal part
        const int M = nrow - ncol;
        if (M > 0) {
          callAndTime_dgemm('N', 'T', M, ncol, jb, -1.0, Qj, ldQ, T.data(), ldT,
                            1.0, R, ldR);
        }
      }

      B_start += nrow * ncol;
    }
  }

#ifdef FINE_TIMING
  DataCollector::get()->sumTime(kTimeDenseFact_schur, clock.stop());
#endif

  return kRetOk;
}
