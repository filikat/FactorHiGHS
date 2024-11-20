#ifndef NUMERIC_H
#define NUMERIC_H

#include <vector>

#include "DataCollector.h"
#include "SolveHandler.h"
#include "Symbolic.h"

class Numeric {
  // columns of factorization, stored by supernode
  std::vector<std::vector<double>> sn_columns_{};

  // swaps of columns for each supernode, ordered locally within a block
  std::vector<std::vector<int>> swaps_{};

  // symbolic object
  const Symbolic& S_;

  // object to handle times and statistics
  DataCollector& DC_;

  // object to handle solve phase in different formats
  std::unique_ptr<SolveHandler> SH_;

  friend class Factorise;

 public:
  // dynamic regularization applied to the matrix
  std::vector<double> total_reg_{};

  Numeric(const Symbolic& S, DataCollector& DC);

  // Full solve
  void solve(std::vector<double>& x) const;
};

#endif