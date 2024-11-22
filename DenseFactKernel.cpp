#include "../ProtoIPM/Regularization.h"
#include "Auxiliary.h"
#include "CallAndTimeBlas.h"
#include "DataCollector.h"
#include "FactorHiGHSSettings.h"
#include "ReturnValues.h"
#include "util/HighsRandom.h"

// Dense Factorization kernel

// #define DEBUG

double regularizePivot(double pivot, double thresh, const int* sign,
                       const double* A, int lda, int j, int n, char uplo,
                       int sn, int bl) {
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

      if (s > 0) {
        required_pivot = std::max(required_pivot, temp);
        required_pivot = std::min(required_pivot, 1e100);
      } else {
        required_pivot = std::min(required_pivot, temp);
        required_pivot = std::max(required_pivot, -1e100);
      }
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

  if (modified_pivot) DataCollector::get()->sumRegPiv();

  return pivot;
}

int denseFactK(char uplo, int n, double* A, int lda, int* pivot_sign,
               double thresh, double* regul, int* swaps, double* pivot_2x2,
               int sn, int bl) {
  // ===========================================================================
  // Factorization kernel
  // Matrix A is in format F
  // BLAS calls: dscal, dcopy, dger
  // ===========================================================================

  // check input
  if (n < 0 || !A || lda < n) {
    printf("\ndenseFactK: invalid input\n");
    return kRetInvalidInput;
  }

  // quick return
  if (n == 0) return kRetOk;

  // ===========================================================================
  // LOWER TRIANGULAR
  // ===========================================================================
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
      Ajj =
          regularizePivot(Ajj, thresh, pivot_sign, A, lda, j, n, uplo, sn, bl);
      regul[j] = std::abs(Ajj - old_pivot);
      DataCollector::get()->setMaxReg(regul[j]);

      // save diagonal element
      A[j + lda * j] = Ajj;

      const int M = n - j - 1;
      if (M > 0) {
        // make copy of column
        callAndTime_dcopy(M, &A[j + 1 + j * lda], 1, temp.data(), 1);

        // scale column j
        callAndTime_dscal(M, 1.0 / Ajj, &A[j + 1 + j * lda], 1);

        // update rest of the matrix
        callAndTime_dger(M, M, -1.0, temp.data(), 1, &A[j + 1 + j * lda], 1,
                         &A[j + 1 + (j + 1) * lda], lda);
      }
    }
  }

  // ===========================================================================
  // UPPER TRIANGULAR
  // ===========================================================================
  else {
    if (!swaps || !pivot_2x2) {
      printf("\ndenseFactK: invalid input\n");
      return kRetInvalidInput;
    }

    // initialize order of pivots
    for (int i = 0; i < n; ++i) swaps[i] = i;

    // allocate space for copy of col(s)
    std::vector<double> temp(n - 1);
    std::vector<double> temp2(n - 1);

    int step = 1;

    for (int j = 0; j < n; j += step) {
      bool flag_2x2 = false;

#ifdef PIVOTING
      // some random pivoting to test
      if (n > 20) {
        HighsRandom hr;
        hr.initialise(j);
        int col = hr.integer(n);
        if (col > j) {
          swapCols('U', n, A, lda, j, col, swaps, pivot_sign);
        }

        if (j < n - 1 && col > n / 2) flag_2x2 = true;
      }
#endif

      // cannot do 2x2 pivoting on last column
      assert(j < n - 1 || flag_2x2 == false);

      if (!flag_2x2) {
        // 1x1 pivots
        step = 1;

        // diagonal element
        double Ajj = A[j + lda * j];

        if (isnan(Ajj)) {
          printf("\ndenseFactK: invalid pivot %e\n", Ajj);
          return kRetInvalidPivot;
        }

        // add regularization
        double old_pivot = Ajj;
        Ajj = regularizePivot(Ajj, thresh, pivot_sign, A, lda, j, n, uplo, sn,
                              bl);
        regul[j] = std::abs(Ajj - old_pivot);
        DataCollector::get()->setMaxReg(regul[j]);

        // save reciprocal of pivot
        A[j + lda * j] = 1.0 / Ajj;

        const int M = n - j - 1;
        if (M > 0) {
          // make copy of row
          callAndTime_dcopy(M, &A[j + (j + 1) * lda], lda, temp.data(), 1);

          // scale row j
          callAndTime_dscal(M, 1.0 / Ajj, &A[j + (j + 1) * lda], lda);

          // update rest of the matrix
          callAndTime_dger(M, M, -1.0, temp.data(), 1, &A[j + (j + 1) * lda],
                           lda, &A[j + 1 + (j + 1) * lda], lda);
        }
      } else {
        // 2x2 pivots
        step = 2;

        // diagonal pivot elements
        const double d1 = A[j + lda * j];
        const double d2 = A[j + 1 + lda * (j + 1)];

        if (isnan(d1) || isnan(d2)) {
          printf("\ndenseFactK: invalid pivot %e %e\n", d1, d2);
          return kRetInvalidPivot;
        }

        // off-diagonal pivot element
        const double offd = A[j + lda * (j + 1)];
        A[j + lda * (j + 1)] = 0.0;

        // compute coefficients of 2x2 inverse
        const double denom = d1 * d2 - offd * offd;
        const double i_d1 = d2 / denom;
        const double i_d2 = d1 / denom;
        const double i_off = -offd / denom;

        // save them in place of pivots
        A[j + lda * j] = i_d1;
        A[j + 1 + lda * (j + 1)] = i_d2;
        pivot_2x2[j] = i_off;

        const int M = n - j - 2;
        if (M > 0) {
          double* r1 = &A[j + (j + 2) * lda];
          double* r2 = &A[j + 1 + (j + 2) * lda];

          // make a copy of first row
          callAndTime_dcopy(M, r1, lda, temp.data(), 1);

          // make a copy of second row
          callAndTime_dcopy(M, r2, lda, temp2.data(), 1);

          // solve rows j,j+1
          callAndTime_dscal(M, i_d1, r1, lda);
          callAndTime_daxpy(M, i_off, temp2.data(), 1, r1, lda);
          callAndTime_dscal(M, i_d2, r2, lda);
          callAndTime_daxpy(M, i_off, temp.data(), 1, r2, lda);

          // update rest of the matrix
          callAndTime_dger(M, M, -1.0, temp.data(), 1, r1, lda,
                           &A[j + 2 + (j + 2) * lda], lda);
          callAndTime_dger(M, M, -1.0, temp2.data(), 1, r2, lda,
                           &A[j + 2 + (j + 2) * lda], lda);
        }
      }
    }
  }

  return kRetOk;
}