#ifndef NUMERIC_H
#define NUMERIC_H

#include <vector>

#include "Auxiliary.h"
#include "Symbolic.h"

class Numeric {
  std::vector<std::vector<double>> SnColumns{};
  const Symbolic* S;

  friend class Factorize;

 public:
  // Forward solve with single right hand side
  void Lsolve(std::vector<double>& x) const;

  // Backward solve with single right hand side
  void Ltsolve(std::vector<double>& x) const;

  // Full solve
  void Solve(std::vector<double>& x) const;
};

#endif