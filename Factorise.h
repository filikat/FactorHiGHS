#ifndef FACTORIZE_H
#define FACTORIZE_H

#include <cmath>

#include "DataCollector.h"
#include "FormatHandler.h"
#include "Numeric.h"
#include "Symbolic.h"
#include "FactorHiGHSSettings.h"

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

  // object to handle times and statistics
  DataCollector& DC_;

  // children in supernodal elimination tree
  std::vector<int> first_child_{};
  std::vector<int> next_child_{};

#ifdef PARALLEL_TREE
  // reverse linked lists of chidlren
  std::vector<int> first_child_reverse_{};
  std::vector<int> next_child_reverse_{};

  // std::vector<int> thr_per_sn{};
#endif

  // generated elements, aka Schur complements.
  std::vector<std::vector<double>> schur_contribution_{};

  // columns of L, stored as dense supernodes
  std::vector<std::vector<double>> sn_columns_{};

  // swaps of columns for each supernode, ordered locally within a block
  std::vector<std::vector<int>> swaps_{};

  // largest diagonal element in the original matrix
  double max_diag_{};
  double min_diag_{};

  // regularization
  std::vector<double> total_reg_{};

  // flag to stop computation
  bool flag_stop_ = false;

 public:
  void permute(const std::vector<int>& iperm);
  void processSupernode(int sn);

 public:
  Factorise(const Symbolic& S, DataCollector& DC, const std::vector<int>& rowsA,
            const std::vector<int>& ptrA, const std::vector<double>& valA);

  bool run(Numeric& num);
};

#endif