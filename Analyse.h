#ifndef ANALYSE_H
#define ANALYSE_H

#include <vector>

#include "Auxiliary.h"
#include "GKlib.h"
#include "Symbolic.h"
#include "metis.h"

// parameters for supernode amalgamation
const int k_start_thresh_relax = 256;
const double k_upper_ratio_relax = 0.02;
const double k_lower_ratio_relax = 0.01;
const int k_max_iter_relax = 10;

// Class to perform the analyse phase of the factorization.
// The final symbolic factorization is stored in an object of type Symbolic.
class Analyse {
  bool ready = false;

  // Matrix to be factorized, stored in upper and lower triangular format
  std::vector<int> rowsUpper{};
  std::vector<int> ptrUpper{};
  std::vector<int> rowsLower{};
  std::vector<int> ptrLower{};
  int n{};
  int nz{};
  double nzL{};
  double operations{};
  double operationsNorelax{};
  double operationsAssembly{};
  FactType type{};

  // Permutation and inverse permutation from Metis
  std::vector<int> perm{};
  std::vector<int> iperm{};

  // Elimination tree
  std::vector<int> parent{};

  // postorder of the elimination tree
  std::vector<int> postorder{};

  // number of entries in each column of L
  std::vector<int> colCount{};

  // sparsity pattern of supernodes of L
  std::vector<int> rowsLsn{};
  std::vector<int> ptrLsn{};

  std::vector<int> snIndices{};

  // fundamental supernodes information
  int snCount{};
  int artificialNz{};
  std::vector<int> snBelong{};
  std::vector<int> snStart{};
  std::vector<int> snParent{};

  // temporary storage for relaxing supernodes
  std::vector<int> fakeNonzeros{};
  std::vector<int> mergedInto{};
  int mergedSn{};

  // relative indices of original columns wrt L columns
  std::vector<int> relindCols{};

  // relative indices of clique wrt parent
  std::vector<std::vector<int>> relindClique{};

  // information about consecutive indices in relindClique
  std::vector<std::vector<int>> consecutiveSums{};

  void GetPermutation();
  void Permute(const std::vector<int>& iperm);
  void ETree();
  void Postorder();
  void ColCount();
  void FundamentalSupernodes();
  void RelaxSupernodes();
  void AfterRelaxSn();
  void SnPattern();
  void RelativeIndCols();
  void RelativeIndClique();
  bool Check() const;

  void PrintTimes() const;

 public:
  // Constructor: matrix must be in upper triangular format
  Analyse(const std::vector<int>& rows_input, const std::vector<int>& ptr_input,
          FactType type_input, const std::vector<int>& order = {});

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

  // save metis iperm to be used by hsl codes for comparison
  std::vector<int> metis_order{};
};

#endif