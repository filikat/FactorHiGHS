#include "../ProtoIPM/Regularization.h"
#include "CallAndTimeBlas.h"
#include "DataCollector.h"
#include "ReturnValues.h"

// Dense Factorization kernel

// #define DEBUG

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

  if (modified_pivot) DC.sumRegPiv();

  return pivot;
}

int denseFactK(char uplo, int n, double* A, int lda, const int* pivot_sign,
               double thresh, double* regul, DataCollector& DC, int sn,
               int bl) {
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
        callAndTime_dger(M, M, -1.0, temp.data(), 1, &A[j + 1 + j * lda], 1,
                         &A[j + 1 + (j + 1) * lda], lda, DC);
      }
    }
  }

  // ===========================================================================
  // UPPER TRIANGULAR
  // ===========================================================================
  else {
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

      const int M = n - j - 1;
      if (M > 0) {
        // make copy of row
        callAndTime_dcopy(M, &A[j + (j + 1) * lda], lda, temp.data(), 1, DC);

        // scale row j
        callAndTime_dscal(M, 1.0 / Ajj, &A[j + (j + 1) * lda], lda, DC);

        // update rest of the matrix
        callAndTime_dger(M, M, -1.0, temp.data(), 1, &A[j + (j + 1) * lda], lda,
                         &A[j + 1 + (j + 1) * lda], lda, DC);
      }
    }
  }

  return kRetOk;
}

void swapCols(char uplo, int n, double* A, int lda, int i, int j,
              std::vector<int>& swaps, DataCollector& DC) {
  // Exchange rows/cols i and j of symmetric matrix A

  // make sure that i < j
  if (i == j) return;
  if (i > j) std::swap(i, j);

  // swap diagonal elements
  std::swap(A[i + i * lda], A[j + j * lda]);

  // swap rest of rows/cols
  if (uplo == 'L') {
    callAndTime_dswap(i, &A[i], lda, &A[j], lda, DC);
    callAndTime_dswap(n - j - 1, &A[j + 1 + i * lda], 1, &A[j + 1 + j * lda], 1,
                      DC);
    callAndTime_dswap(j - i - 1, &A[i + 1 + i * lda], 1, &A[j + (i + 1) * lda],
                      lda, DC);
  } else {
    callAndTime_dswap(i, &A[i * lda], 1, &A[j * lda], 1, DC);
    callAndTime_dswap(n - j - 1, &A[i + (j + 1) * lda], lda,
                      &A[j + (j + 1) * lda], lda, DC);
    callAndTime_dswap(j - i - 1, &A[i + (i + 1) * lda], lda,
                      &A[i + 1 + j * lda], 1, DC);
  }

  // keep track of order of swaps
  swaps[i] = j;
}