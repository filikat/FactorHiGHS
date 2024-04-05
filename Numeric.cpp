#include "Numeric.h"

void Numeric::Lsolve(std::vector<double>& x) const {
  // Forward solve
  // (row,col) indicate absolute indices within L.
  // (i,j) indicate relative indices within supernode.

  for (int sn = 0; sn < S->Sn(); ++sn) {
    // leading size of supernode
    int ldSn = S->Ptr(sn + 1) - S->Ptr(sn);

    for (int col = S->Sn_start(sn); col < S->Sn_start(sn + 1); ++col) {
      // relative index of column within supernode
      int j = col - S->Sn_start(sn);

      // scale entry x_col
      x[col] /= SnColumns[sn][j + ldSn * j];

      // index to access S->rows for this supernode
      int startRow = S->Ptr(sn);

      // update x with entries in column j
      for (int i = j + 1; i < ldSn; ++i) {
        // absolute index of row to update
        int row = S->Rows(startRow + i);

        x[row] -= x[col] * SnColumns[sn][i + ldSn * j];
      }
    }
  }
}

void Numeric::Ltsolve(std::vector<double>& x) const {
  // Backward solve
  // (row,col) indicate absolute indices within L.
  // (i,j) indicate relative indices within supernode.

  // go thourh the sn in reverse order
  for (int sn = S->Sn() - 1; sn >= 0; --sn) {
    // leading size of supernode
    int ldSn = S->Ptr(sn + 1) - S->Ptr(sn);

    // go through the columns of the sn in reverse order
    for (int col = S->Sn_start(sn + 1) - 1; col >= S->Sn_start(sn); --col) {
      // relative index of column within supernode
      int j = col - S->Sn_start(sn);

      // index to access S->rows for this supernode
      int startRow = S->Ptr(sn);

      // update x with entries in column j
      for (int i = j + 1; i < ldSn; ++i) {
        // absolute index of row to update
        int row = S->Rows(startRow + i);

        x[col] -= x[row] * SnColumns[sn][i + ldSn * j];
      }

      // scale entry x_col
      x[col] /= SnColumns[sn][j + ldSn * j];
    }
  }
}

void Numeric::Solve(std::vector<double>& x) const {
  PermuteVector(x, S->Perm());
  Lsolve(x);
  Ltsolve(x);
  PermuteVector(x, S->Iperm());
}
