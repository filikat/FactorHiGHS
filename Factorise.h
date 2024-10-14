#ifndef FACTORIZE_H
#define FACTORIZE_H

#include <cmath>

#include "Auxiliary.h"
#include "Blas_declaration.h"
#include "CliqueStack.h"
#include "CurtisReidScalingSym.h"
#include "DenseFact_declaration.h"
#include "FormatHandler.h"
#include "Numeric.h"
#include "ReturnValues.h"
#include "Symbolic.h"

class Factorise {
 public:
  // matrix to factorise
  std::vector<int> rowsA_{};
  std::vector<int> ptrA_{};
  std::vector<double> valA_{};
  int n_{};
  int nzA_{};

  // symbolic factorisation
  const Symbolic& S_;

  // interface to specific format used for dense matrices
  FormatHandler* FH_;

  // children in supernodal elimination tree
  std::vector<int> first_children_{};
  std::vector<int> next_children_{};

  // stack of Schur complements
  CliqueStack& clique_stack_;

  // columns of L, stored as dense supernodes
  std::vector<std::vector<double>> sn_columns_{};

  // largest diagonal element in the original matrix
  double max_diag_{};
  double min_diag_{};

  // symmetric scaling to apply to the original matrix
  std::vector<double> colscale_{};
  std::vector<int> colexp_{};

  // regularization
  std::vector<double> total_reg_{};

  // number of pivots that received dynamic regularization
  int n_reg_piv_{};

 public:
  void permute(const std::vector<int>& iperm);
  int processSupernode(int sn);
  bool check() const;
  void equilibrate();
  void scale();

  // extreme values of the factorisation
  double maxD_ = 0.0;
  double minD_ = std::numeric_limits<double>::max();
  double maxoffD_ = 0.0;
  double minoffD_ = std::numeric_limits<double>::max();

 public:
  Factorise(const Symbolic& S, CliqueStack& stack,
            const std::vector<int>& rowsA, const std::vector<int>& ptrA,
            const std::vector<double>& valA);

  int run(Numeric& num);
};

#endif