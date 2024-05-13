#ifndef NUMERIC_H
#define NUMERIC_H

#include <vector>

#include "Auxiliary.h"
#include "Blas_declaration.h"
#include "DenseFact_declaration.h"
#include "Symbolic.h"

class Numeric {
  std::vector<std::vector<double>> SnColumns{};
  const Symbolic* S;

  int nb = hybridBlockSize;

  friend class Factorise;

 public:
  // Forward solve with single right hand side
  void Lsolve(std::vector<double>& x) const;

  // Backward solve with single right hand side
  void Ltsolve(std::vector<double>& x) const;

  // Diagonal solve for LDL
  void Dsolve(std::vector<double>& x) const;

  // Full solve
  void Solve(std::vector<double>& x) const;
};

#endif