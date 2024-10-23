#ifndef NUMERIC_H
#define NUMERIC_H

#include <vector>

#include "Auxiliary.h"
#include "Blas_declaration.h"
#include "Symbolic.h"
#include "Timing.h"
#include "DataCollector.h"

class Numeric {
  // columns of factorization, stored by supernode
  std::vector<std::vector<double>> sn_columns_{};

  // symbolic object
  const Symbolic& S_;

  // object to handle times and statistics
  DataCollector& DC_;

  // scaling applied to the matrix
  std::vector<double> colscale_{};

  friend class Factorise;

 public:
  // dynamic regularization applied to the matrix
  std::vector<double> total_reg_{};

  Numeric(const Symbolic& S, DataCollector& DC);

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