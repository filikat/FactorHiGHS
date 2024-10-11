#ifndef NUMERIC_H
#define NUMERIC_H

#include <vector>

#include "Auxiliary.h"
#include "Blas_declaration.h"
#include "DenseFact_declaration.h"
#include "Symbolic.h"
#include "timing.h"

class Numeric {
  // columns of factorization, stored by supernode
  std::vector<std::vector<double>> sn_columns_{};

  // symbolic object
  const Symbolic* S_;

  // scaling applied to the matrix
  std::vector<double> colscale_{};
  std::vector<int> colexp_{};

  friend class Factorise;

 public:
  // dynamic regularization applied to the matrix
  std::vector<double> total_reg_{};

  // number of pivots that received dynamic regularization
  int n_reg_piv_{};

  // Forward solve with single right hand side
  void forwardSolve(std::vector<double>& x) const;

  // Backward solve with single right hand side
  void backwardSolve(std::vector<double>& x) const;

  // Diagonal solve for LDL
  void diagSolve(std::vector<double>& x) const;

  // Full solve
  void solve(std::vector<double>& x) const;
};

#endif