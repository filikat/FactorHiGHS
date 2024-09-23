#ifndef FACTORIZE_H
#define FACTORIZE_H

#include <cmath>

#include "Auxiliary.h"
#include "Blas_declaration.h"
#include "DenseFact_declaration.h"
#include "FormatHandler.h"
#include "Numeric.h"
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

  // generated elements, aka Schur complements.
  // NB: they cannot be stored as std::vector, since they need to be
  // instantiated without being initialized to zero.
  std::vector<double*> schur_contribution_{};

  // columns of L, stored as dense supernodes
  std::vector<std::vector<double>> sn_columns_{};

  // starting point of the diagonal blocks for hybrid format
  std::vector<std::vector<int>> clique_block_start_{};

  // largest diagonal element in the original matrix
  double max_diag{};

 public:
  void permute(const std::vector<int>& iperm);
  int processSupernode(int sn);
  bool check() const;

 public:
  Factorise(const Symbolic& S, const std::vector<int>& rowsA,
            const std::vector<int>& ptrA, const std::vector<double>& valA);

  int run(Numeric& num);
};

#endif