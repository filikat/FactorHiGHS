#include "Numeric.h"

void Numeric::Lsolve(std::vector<double>& x) const {
  // Forward solve.
  // Uses Blas 2

  // variables for BLAS calls
  char LL = 'L';
  char NN = 'N';
  int i_one = 1;
  double d_one = 1.0;
  double d_zero = 0.0;

  for (int sn = 0; sn < S->Sn(); ++sn) {
    // leading size of supernode
    int ldSn = S->Ptr(sn + 1) - S->Ptr(sn);

    // number of columns in the supernode
    int sn_size = S->SnStart(sn + 1) - S->SnStart(sn);

    // first colums of the supernode
    int sn_start = S->SnStart(sn);

    // size of clique of supernode
    int clique_size = ldSn - sn_size;

    // index to access S->rows for this supernode
    int start_row = S->Ptr(sn);

    dtrsv(&LL, &NN, &NN, &sn_size, SnColumns[sn].data(), &ldSn, &x[sn_start],
          &i_one);

    // temporary space for gemv
    std::vector<double> y(clique_size);

    dgemv(&NN, &clique_size, &sn_size, &d_one, &SnColumns[sn][sn_size], &ldSn,
          &x[sn_start], &i_one, &d_zero, y.data(), &i_one);

    // scatter solution of gemv
    for (int i = 0; i < clique_size; ++i) {
      int row = S->Rows(start_row + sn_size + i);
      x[row] -= y[i];
    }
  }
}

void Numeric::Ltsolve(std::vector<double>& x) const {
  // Backward solve.
  // Uses Blas 2

  // variables for BLAS calls
  char LL = 'L';
  char NN = 'N';
  char TT = 'T';
  int i_one = 1;
  double d_m_one = -1.0;
  double d_one = 1.0;

  // go through the sn in reverse order
  for (int sn = S->Sn() - 1; sn >= 0; --sn) {
    // leading size of supernode
    int ldSn = S->Ptr(sn + 1) - S->Ptr(sn);

    // number of columns in the supernode
    int sn_size = S->SnStart(sn + 1) - S->SnStart(sn);

    // first colums of the supernode
    int sn_start = S->SnStart(sn);

    // size of clique of supernode
    int clique_size = ldSn - sn_size;

    // index to access S->rows for this supernode
    int start_row = S->Ptr(sn);

    // temporary space for gemv
    std::vector<double> y(clique_size);

    // scatter entries into y
    for (int i = 0; i < clique_size; ++i) {
      int row = S->Rows(start_row + sn_size + i);
      y[i] = x[row];
    }

    dgemv(&TT, &clique_size, &sn_size, &d_m_one, &SnColumns[sn][sn_size], &ldSn,
          y.data(), &i_one, &d_one, &x[sn_start], &i_one);

    dtrsv(&LL, &TT, &NN, &sn_size, SnColumns[sn].data(), &ldSn, &x[sn_start],
          &i_one);
  }
}

void Numeric::Solve(std::vector<double>& x) const {
  PermuteVector(x, S->Perm());
  Lsolve(x);
  Ltsolve(x);
  PermuteVector(x, S->Iperm());
}
