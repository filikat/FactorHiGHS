#ifndef ANALYSE_H
#define ANALYSE_H

#include <algorithm>
#include <vector>

#include "Auxiliary.h"
#include "GKlib.h"
#include "ReturnValues.h"
#include "Symbolic.h"
#include "metis.h"

// parameters for supernode amalgamation
const int kStartThreshRelax = 256;
const double kUpperRatioRelax = 0.02;
const double kLowerRatioRelax = 0.01;
const int kMaxIterRelax = 10;

// Class to perform the analyse phase of the factorization.
// The final symbolic factorization is stored in an object of type Symbolic.
class Analyse {
  bool ready_ = false;

  // Matrix to be factorized, stored in upper and lower triangular format
  std::vector<int> rows_upper_{};
  std::vector<int> ptr_upper_{};
  std::vector<int> rows_lower_{};
  std::vector<int> ptr_lower_{};
  int n_{};
  int nz_{};
  double nz_factor_{};
  double operations_{};
  double operations_no_relax_{};
  double operations_assembly_{};
  int negative_pivots_{};

  // Permutation and inverse permutation from Metis
  std::vector<int> perm_{};
  std::vector<int> iperm_{};

  // Elimination tree
  std::vector<int> parent_{};

  // postorder of the elimination tree
  std::vector<int> postorder_{};

  // number of entries in each column of L
  std::vector<int> col_count_{};

  // sparsity pattern of supernodes of L
  std::vector<int> rows_sn_{};
  std::vector<int> ptr_sn_{};

  std::vector<int> sn_indices_{};

  // fundamental supernodes information
  int sn_count_{};
  int artificial_nz_{};
  std::vector<int> sn_belong_{};
  std::vector<int> sn_start_{};
  std::vector<int> sn_parent_{};

  // temporary storage for relaxing supernodes
  std::vector<int> fake_nz_{};
  std::vector<int> merged_into_{};
  int merged_sn_{};

  // relative indices of original columns wrt L columns
  std::vector<int> relind_cols_{};

  // relative indices of clique wrt parent
  std::vector<std::vector<int>> relind_clique_{};

  // information about consecutive indices in relindClique
  std::vector<std::vector<int>> consecutive_sums_{};

  // estimate of maximum storage
  double max_storage_{};
  double max_stack_entries_{};

  // data for parallelization
  std::vector<std::vector<int>> subtrees_per_thread_{};
  std::vector<double> ops_per_thread_{};
  int sn_above_layer0_{};

  // symbolic object where to store result
  Symbolic& S_;

  // Functions to perform analyse phase
  int getPermutation();
  void permute(const std::vector<int>& iperm);
  void eTree();
  void postorder();
  void colCount();
  void fundamentalSupernodes();
  void relaxSupernodes();
  void relaxSupernodes2();
  void afterRelaxSn();
  void snPattern();
  void relativeIndCols();
  void relativeIndClique();
  bool check() const;
  void generateLayer0(double imbalance_ratio);
  void reorderChildren();
  void computeStorage(int fr, int sz, double& fr_entries,
                      double& cl_entries) const;

 public:
  // Constructor: matrix must be in lower triangular format
  Analyse(const std::vector<int>& rows, const std::vector<int>& ptr,
          Symbolic& S, const std::vector<int>& order = {},
          int negative_pivots = 0);

  // Run analyse phase and save the result in Symbolic object S_
  int run();

  // save metis iperm to be used by hsl codes for comparison
  std::vector<int> metis_order_{};
};

#endif