#ifndef NUMERIC_H
#define NUMERIC_H

#include <vector>

#include "Auxiliary.h"
#include "Blas_declaration.h"
#include "DenseFact_declaration.h"
#include "Symbolic.h"
#include "timing.h"

class Numeric {
  std::vector<std::vector<double>> sn_columns_{};
  const Symbolic* S_;

  friend class Factorise;

 public:
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