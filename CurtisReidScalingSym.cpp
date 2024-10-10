#include "CurtisReidScalingSym.h"

#include "../ProtoIPM/VectorOperations.h"

void product(const std::vector<double>& x, std::vector<double>& y,
             const std::vector<int>& ptr, const std::vector<int>& rows,
             const std::vector<double>& N) {
  // Multiply by matrix E, i.e. matrix A with all entries equal to one, in lower
  // triangular form, and sum component-wise product of N with x.
  // E * x + N .* x= y

  int n = x.size();

  y.assign(n, 0.0);

  // multiply by E
  for (int col = 0; col < n; ++col) {
    for (int el = ptr[col]; el < ptr[col + 1]; ++el) {
      int row = rows[el];
      y[row] += x[col];
      if (row != col) y[col] += x[row];
    }
  }

  // multiply by N
  for (int i = 0; i < n; ++i) y[i] += N[i] * x[i];
}

void CG_for_CR_scaling(const std::vector<double>& b, std::vector<double>& x,
                       const std::vector<double>& N,
                       const std::vector<int>& ptr,
                       const std::vector<int>& rows) {
  int n = N.size();

  // initial residual
  std::vector<double> r = b;

  // initial approximation
  x.assign(n, 0.0);

  // direction
  std::vector<double> p = r;

  int iter{};
  std::vector<double> Ap(n);

  while (iter < 50) {
    product(p, Ap, ptr, rows, N);

    double norm_r = dotProd(r, r);
    double alpha = norm_r / dotProd(p, Ap);

    // x = x + alpha * p
    vectorAdd(x, p, alpha);

    // r = r - alpha * Ap;
    vectorAdd(r, Ap, -alpha);

    // exit test
    if (norm2(r) / norm2(b) < 1e-6) break;

    double beta = dotProd(r, r) / norm_r;

    // p = r + beta * p
    vectorAdd(p, r, 1.0, beta);

    ++iter;
  }
}

void CurtisReidScalingSym(const std::vector<int>& ptr,
                          const std::vector<int>& rows,
                          const std::vector<double>& val,
                          std::vector<int>& colexp) {
  // Takes as input the CSC matrix A.
  // Computes Curtis-Reid scaling exponents of the matrix, using powers of 2.

  int n = colexp.size();

  // rhs for CG
  std::vector<double> logsumcol(n, 0.0);

  // number of entries in each column
  std::vector<double> col_entries(n, 0.0);

  // log A_ij
  for (int col = 0; col < n; ++col) {
    for (int el = ptr[col]; el < ptr[col + 1]; ++el) {
      int row = rows[el];
      if (val[el] != 0.0) {
        double temp = log2(std::abs(val[el]));
        logsumcol[col] += temp;
        col_entries[col] += 1.0;

        // only lower triangle is used, so add components corresponding to the
        // upper triangle
        if (col != row) {
          logsumcol[row] += temp;
          col_entries[row] += 1.0;
        }
      }
    }
  }

  // solve linear system with CG
  std::vector<double> exponents(n);
  CG_for_CR_scaling(logsumcol, exponents, col_entries, ptr, rows);

  // round exponents
  for (int j = 0; j < n; ++j) colexp[j] = -std::round(exponents[j]);
}
