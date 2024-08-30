#ifndef FACTORIZE_H
#define FACTORIZE_H

#include <cmath>

#include "Auxiliary.h"
#include "Blas_declaration.h"
#include "DenseFact_declaration.h"
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

  // children in supernodal elimination tree
  std::vector<int> first_children_{};
  std::vector<int> next_children_{};

  // generated elements
  std::vector<double*> schur_contribution_{};

  // columns of L, stored as dense supernodes
  std::vector<std::vector<double>> sn_columns_{};

  std::vector<std::vector<int>> clique_block_start_{};

 public:
  void permute(const std::vector<int>& iperm);
  int processSupernode(int sn);
  bool check() const;
  void printTimes() const;

 public:
  Factorise(const Symbolic& S, const std::vector<int>& rowsA,
            const std::vector<int>& ptrA, const std::vector<double>& valA);

  int run(Numeric& num);

  std::vector<double> time_per_sn_{};

  double time_prepare_{};
  double time_assemble_original_{};
  double time_assemble_children_F_{};
  double time_assemble_children_C_{};
  double time_factorise_{};
  double time_total_{};
  std::vector<double> times_dense_fact_;
};

#endif