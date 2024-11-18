#include "Numeric.h"

#include "Auxiliary.h"
#include "CallAndTimeBlas.h"
#include "Timing.h"
#include "FormatHandler.h"

Numeric::Numeric(const Symbolic& S, DataCollector& DC) : S_{S}, DC_{DC} {}

void Numeric::forwardSolve(std::vector<double>& x) const {
  // Forward solve.
  // Blas calls: dtrsv, dtpsv, dgemv

  // variables for BLAS calls
  const char c_L = 'L';
  const char c_N = 'N';
  const char c_T = 'T';
  const char c_U = 'U';
  const int i_one = 1;
  const double d_one = 1.0;
  const double d_zero = 0.0;

  if (S_.formatType() == FormatType::HybridPacked ||
      S_.formatType() == FormatType::HybridHybrid) {
    // supernode columns in format FH

    const int nb = S_.blockSize();

    for (int sn = 0; sn < S_.sn(); ++sn) {
      // leading size of supernode
      const int ldSn = S_.ptr(sn + 1) - S_.ptr(sn);

      // number of columns in the supernode
      const int sn_size = S_.snStart(sn + 1) - S_.snStart(sn);

      // first colums of the supernode
      const int sn_start = S_.snStart(sn);

      // index to access S->rows for this supernode
      const int start_row = S_.ptr(sn);

      // number of blocks of columns
      const int n_blocks = (sn_size - 1) / nb + 1;

      // index to access snColumns[sn]
      int SnCol_ind{};

      // go through blocks of columns for this supernode
      for (int j = 0; j < n_blocks; ++j) {
        // number of columns in the block
        const int jb = std::min(nb, sn_size - nb * j);

        // number of entries in diagonal part
        const int diag_entries = jb * jb;

        // index to access vector x
        const int x_start = sn_start + nb * j;

        callAndTime_dtrsv('U', 'T', 'U', jb, &sn_columns_[sn][SnCol_ind], jb,
                          &x[x_start], 1, DC_);

        SnCol_ind += diag_entries;

        // temporary space for gemv
        const int gemv_space = ldSn - nb * j - jb;
        std::vector<double> y(gemv_space);

        callAndTime_dgemv('T', jb, gemv_space, 1.0, &sn_columns_[sn][SnCol_ind],
                          jb, &x[x_start], 1, 0.0, y.data(), 1, DC_);
        SnCol_ind += jb * gemv_space;

        // scatter solution of gemv
        for (int i = 0; i < gemv_space; ++i) {
          const int row = S_.rows(start_row + nb * j + jb + i);
          x[row] -= y[i];
        }
      }
    }

  } else if (S_.formatType() == FormatType::PackedPacked) {
    // supernode columns in format FP

    const int nb = S_.blockSize();

    for (int sn = 0; sn < S_.sn(); ++sn) {
      // leading size of supernode
      const int ldSn = S_.ptr(sn + 1) - S_.ptr(sn);

      // number of columns in the supernode
      const int sn_size = S_.snStart(sn + 1) - S_.snStart(sn);

      // first colums of the supernode
      const int sn_start = S_.snStart(sn);

      // index to access S->rows for this supernode
      const int start_row = S_.ptr(sn);

      // number of blocks of columns
      const int n_blocks = (sn_size - 1) / nb + 1;

      // index to access snColumns[sn]
      int SnCol_ind{};

      // go through blocks of columns for this supernode
      for (int j = 0; j < n_blocks; ++j) {
        // number of columns in the block
        const int jb = std::min(nb, sn_size - nb * j);

        // index to access vector x
        const int x_start = sn_start + nb * j;

        const int lda = ldSn - j * nb;
        callAndTime_dtrsv('L', 'N', 'U', jb, &sn_columns_[sn][SnCol_ind], lda,
                          &x[x_start], 1, DC_);

        // temporary space for gemv
        const int gemv_space = ldSn - nb * j - jb;
        std::vector<double> y(gemv_space);

        callAndTime_dgemv('N', gemv_space, jb, 1.0,
                          &sn_columns_[sn][SnCol_ind + jb], lda, &x[x_start], 1,
                          0.0, y.data(), 1, DC_);

        SnCol_ind += jb * jb + jb * gemv_space;

        // scatter solution of gemv
        for (int i = 0; i < gemv_space; ++i) {
          const int row = S_.rows(start_row + nb * j + jb + i);
          x[row] -= y[i];
        }
      }
    }
  } else {
    // supernode columns in format F

    for (int sn = 0; sn < S_.sn(); ++sn) {
      // leading size of supernode
      const int ldSn = S_.ptr(sn + 1) - S_.ptr(sn);

      // number of columns in the supernode
      const int sn_size = S_.snStart(sn + 1) - S_.snStart(sn);

      // first colums of the supernode
      const int sn_start = S_.snStart(sn);

      // size of clique of supernode
      const int clique_size = ldSn - sn_size;

      // index to access S->rows for this supernode
      const int start_row = S_.ptr(sn);

      callAndTime_dtrsv('L', 'N', 'U', sn_size, sn_columns_[sn].data(), ldSn,
                        &x[sn_start], 1, DC_);

      // temporary space for gemv
      std::vector<double> y(clique_size);

      callAndTime_dgemv('N', clique_size, sn_size, 1.0,
                        &sn_columns_[sn][sn_size], ldSn, &x[sn_start], 1, 0.0,
                        y.data(), 1, DC_);

      // scatter solution of gemv
      for (int i = 0; i < clique_size; ++i) {
        const int row = S_.rows(start_row + sn_size + i);
        x[row] -= y[i];
      }
    }
  }
}

void Numeric::backwardSolve(std::vector<double>& x) const {
  // Backward solve.
  // Blas calls: dtrsv, dtpsv, dgemv

  // variables for BLAS calls
  const char c_L = 'L';
  const char c_N = 'N';
  const char c_T = 'T';
  const char c_U = 'U';
  const int i_one = 1;
  const double d_m_one = -1.0;
  const double d_one = 1.0;
  const double d_zero = 0.0;

  if (S_.formatType() == FormatType::HybridPacked ||
      S_.formatType() == FormatType::HybridHybrid) {
    // supernode columns in format FH

    const int nb = S_.blockSize();

    // go through the sn in reverse order
    for (int sn = S_.sn() - 1; sn >= 0; --sn) {
      // leading size of supernode
      const int ldSn = S_.ptr(sn + 1) - S_.ptr(sn);

      // number of columns in the supernode
      const int sn_size = S_.snStart(sn + 1) - S_.snStart(sn);

      // first colums of the supernode
      const int sn_start = S_.snStart(sn);

      // index to access S->rows for this supernode
      const int start_row = S_.ptr(sn);

      // number of blocks of columns
      const int n_blocks = (sn_size - 1) / nb + 1;

      // index to access snColumns[sn]
      // initialized with the total number of entries of snColumns[sn]
      int SnCol_ind = sn_columns_[sn].size() - extra_space;

      // go through blocks of columns for this supernode in reverse order
      for (int j = n_blocks - 1; j >= 0; --j) {
        // number of columns in the block
        const int jb = std::min(nb, sn_size - nb * j);

        // number of entries in diagonal part
        const int diag_entries = jb * jb;

        // index to access vector x
        const int x_start = sn_start + nb * j;

        // temporary space for gemv
        const int gemv_space = ldSn - nb * j - jb;
        std::vector<double> y(gemv_space);

        // scatter entries into y
        for (int i = 0; i < gemv_space; ++i) {
          const int row = S_.rows(start_row + nb * j + jb + i);
          y[i] = x[row];
        }

        SnCol_ind -= jb * gemv_space;
        callAndTime_dgemv('N', jb, gemv_space, -1.0,
                          &sn_columns_[sn][SnCol_ind], jb, y.data(), 1, 1.0,
                          &x[x_start], 1, DC_);

        SnCol_ind -= diag_entries;
        callAndTime_dtrsv('U', 'N', 'U', jb, &sn_columns_[sn][SnCol_ind], jb,
                          &x[x_start], 1, DC_);
      }
    }
  } else if (S_.formatType() == FormatType::PackedPacked) {
    // supernode columns in format FP

    const int nb = S_.blockSize();

    // go through the sn in reverse order
    for (int sn = S_.sn() - 1; sn >= 0; --sn) {
      // leading size of supernode
      const int ldSn = S_.ptr(sn + 1) - S_.ptr(sn);

      // number of columns in the supernode
      const int sn_size = S_.snStart(sn + 1) - S_.snStart(sn);

      // first colums of the supernode
      const int sn_start = S_.snStart(sn);

      // index to access S->rows for this supernode
      const int start_row = S_.ptr(sn);

      // number of blocks of columns
      const int n_blocks = (sn_size - 1) / nb + 1;

      // index to access snColumns[sn]
      // initialized with the total number of entries of snColumns[sn]
      std::vector<int> diag_start;
      getDiagStart(ldSn, sn_size, nb, n_blocks, diag_start);

      // go through blocks of columns for this supernode in reverse order
      for (int j = n_blocks - 1; j >= 0; --j) {
        // number of columns in the block
        const int jb = std::min(nb, sn_size - nb * j);

        // number of entries in diagonal part
        const int diag_entries = jb * (jb + 1) / 2;

        // index to access vector x
        const int x_start = sn_start + nb * j;

        // temporary space for gemv
        const int gemv_space = ldSn - nb * j - jb;
        std::vector<double> y(gemv_space);

        // scatter entries into y
        for (int i = 0; i < gemv_space; ++i) {
          const int row = S_.rows(start_row + nb * j + jb + i);
          y[i] = x[row];
        }

        const int lda = ldSn - nb * j;
        callAndTime_dgemv('T', gemv_space, jb, -1.0,
                          &sn_columns_[sn][diag_start[j] + jb], lda, y.data(),
                          1, 1.0, &x[x_start], 1, DC_);

        callAndTime_dtrsv('L', 'T', 'U', jb, &sn_columns_[sn][diag_start[j]],
                          lda, &x[x_start], 1, DC_);
      }
    }
  } else {
    // supernode columns in format F

    // go through the sn in reverse order
    for (int sn = S_.sn() - 1; sn >= 0; --sn) {
      // leading size of supernode
      const int ldSn = S_.ptr(sn + 1) - S_.ptr(sn);

      // number of columns in the supernode
      const int sn_size = S_.snStart(sn + 1) - S_.snStart(sn);

      // first colums of the supernode
      const int sn_start = S_.snStart(sn);

      // size of clique of supernode
      const int clique_size = ldSn - sn_size;

      // index to access S->rows for this supernode
      const int start_row = S_.ptr(sn);

      // temporary space for gemv
      std::vector<double> y(clique_size);

      // scatter entries into y
      for (int i = 0; i < clique_size; ++i) {
        const int row = S_.rows(start_row + sn_size + i);
        y[i] = x[row];
      }

      callAndTime_dgemv('T', clique_size, sn_size, -1.0,
                        &sn_columns_[sn][sn_size], ldSn, y.data(), 1, 1.0,
                        &x[sn_start], 1, DC_);

      callAndTime_dtrsv('L', 'T', 'U', sn_size, sn_columns_[sn].data(), ldSn,
                        &x[sn_start], 1, DC_);
    }
  }
}

void Numeric::diagSolve(std::vector<double>& x) const {
  // Diagonal solve

  if (S_.formatType() == FormatType::HybridPacked ||
      S_.formatType() == FormatType::HybridHybrid) {
    // supernode columns in format FH

    const int nb = S_.blockSize();

    for (int sn = 0; sn < S_.sn(); ++sn) {
      // leading size of supernode
      const int ldSn = S_.ptr(sn + 1) - S_.ptr(sn);

      // number of columns in the supernode
      const int sn_size = S_.snStart(sn + 1) - S_.snStart(sn);

      // first colums of the supernode
      const int sn_start = S_.snStart(sn);

      // number of blocks of columns
      const int n_blocks = (sn_size - 1) / nb + 1;

      // index to access diagonal part of block
      int diag_start{};

      // go through blocks of columns for this supernode
      for (int j = 0; j < n_blocks; ++j) {
        // number of columns in the block
        const int jb = std::min(nb, sn_size - nb * j);

        // go through columns of block
        for (int col = 0; col < jb; ++col) {
          const double d = sn_columns_[sn][diag_start + (jb + 1) * col];
          x[sn_start + nb * j + col] /= d;
        }

        // move diag_start forward by number of diagonal entries in block
        diag_start += jb * jb;

        // move diag_start forward by number of sub-diagonal entries in block
        diag_start += (ldSn - nb * j - jb) * jb;
      }
    }
  } else if (S_.formatType() == FormatType::PackedPacked) {
    // supernode columns in format FP

    const int nb = S_.blockSize();

    for (int sn = 0; sn < S_.sn(); ++sn) {
      // leading size of supernode
      const int ldSn = S_.ptr(sn + 1) - S_.ptr(sn);

      // number of columns in the supernode
      const int sn_size = S_.snStart(sn + 1) - S_.snStart(sn);

      // first colums of the supernode
      const int sn_start = S_.snStart(sn);

      // number of blocks of columns
      const int n_blocks = (sn_size - 1) / nb + 1;

      // index to access diagonal part of block
      int diag_start{};

      // go through blocks of columns for this supernode
      for (int j = 0; j < n_blocks; ++j) {
        // number of columns in the block
        const int jb = std::min(nb, sn_size - nb * j);

        const int lda = ldSn - nb * j;
        // go through columns of block
        for (int col = 0; col < jb; ++col) {
          const double d = sn_columns_[sn][diag_start + col + lda * col];
          x[sn_start + nb * j + col] /= d;
        }

        // move diag_start forward by number of diagonal entries in block
        diag_start += jb * jb;

        // move diag_start forward by number of sub-diagonal entries in block
        diag_start += (ldSn - nb * j - jb) * jb;
      }
    }
  } else {
    // supernode columns in format F

    for (int sn = 0; sn < S_.sn(); ++sn) {
      // leading size of supernode
      const int ldSn = S_.ptr(sn + 1) - S_.ptr(sn);

      for (int col = S_.snStart(sn); col < S_.snStart(sn + 1); ++col) {
        // relative index of column within supernode
        const int j = col - S_.snStart(sn);

        // diagonal entry of column j
        const double d = sn_columns_[sn][j + j * ldSn];

        x[col] /= d;
      }
    }
  }
}

void Numeric::solve(std::vector<double>& x) const {
#ifdef COARSE_TIMING
  Clock clock{};
  clock.start();
#endif

  // permute rhs
  permuteVectorInverse(x, S_.iperm());

  // solve
  forwardSolve(x);
  diagSolve(x);
  backwardSolve(x);

  // unpermute solution
  permuteVector(x, S_.iperm());

#ifdef COARSE_TIMING
  DC_.sumTime(kTimeSolve, clock.stop());
#endif
}
