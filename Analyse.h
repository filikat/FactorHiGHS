#ifndef ANALYSE_H
#define ANALYSE_H

#include <vector>

#include "Auxiliary.h"
#include "GKlib.h"
#include "Symbolic.h"
#include "metis.h"

// parameters for supernode amalgamation
const int max_artificial_nz = 1024;
const int small_sn_thresh = 16;

// Class to perform the analyse phase of the factorization.
// The final symbolic factorization is stored in an object of type Symbolic.
class Analyse {
  bool ready = false;

 public:
  // Matrix to be factorized, stored in upper and lower triangular format
  std::vector<int> rowsUpper{};
  std::vector<int> ptrUpper{};
  std::vector<int> rowsLower{};
  std::vector<int> ptrLower{};
  int n{};
  int nz{};
  int nzL{};
  double operations{};
  double operations_norelax{};
  double operations_assembly{};

  // Permutation and inverse permutation from Metis
  std::vector<int> perm{};
  std::vector<int> iperm{};

  // Elimination tree
  std::vector<int> parent{};

  // postorder of the elimination tree
  std::vector<int> postorder{};

  // number of entries in each column of L
  std::vector<int> colcount{};
  std::vector<int> rowcount{};

  // sparsity pattern of supernodes of L
  std::vector<int> rowsLsn{};
  std::vector<int> ptrLsn{};

  std::vector<int> sn_indices{};

  // fundamental supernodes information
  int sn_count{};
  int artificialNz{};
  std::vector<int> sn_belong{};
  std::vector<int> sn_start{};
  std::vector<int> sn_parent{};

  // temporary storage for relaxing supernodes
  std::vector<int> fake_nonzeros{};
  std::vector<int> mergedInto{};
  int merged_sn{};

  // relative indices of original columns wrt L columns
  std::vector<int> relind_cols{};

  // relative indices of clique wrt parent
  std::vector<std::vector<int>> relind_clique{};

  // information about consecutive indices in relind_clique
  std::vector<std::vector<int>> consecutiveSums{};

  void GetPermutation();
  void Permute(const std::vector<int>& iperm);
  void ETree();
  void Postorder();
  void RowColCount();
  void ColCount();
  void FundamentalSupernodes();
  void RelaxSupernodes();
  void AfterRelaxSn();
  void RelaxSupernodes_2();
  void RelaxSupernodes_3();
  void SnPattern();
  void RelativeInd_cols();
  void RelativeInd_clique();
  void Clear();
  bool Check() const;

  void PrintTimes() const;

 public:
  // Constructor: matrix must be in upper triangular format
  Analyse(const std::vector<int>& rows_input, const std::vector<int>& ptr_input,
          const std::vector<int>& order = {});

  // Run analyse phase and save the result in Symbolic object S
  void Run(Symbolic& S);

  // times
  double time_metis{};
  double time_tree{};
  double time_count{};
  double time_pattern{};
  double time_sn{};
  double time_relind{};
  double time_total{};

  std::vector<int> metis_order{};
};

#endif