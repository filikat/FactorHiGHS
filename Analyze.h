#ifndef ANALYZE_H
#define ANALYZE_H

#include <vector>

#include "Auxiliary.h"
#include "GKlib.h"
#include "Symbolic.h"
#include "metis.h"

// Class to perform the analyze phase of the factorization.
// The final symbolic factorization is stored in an object of type Symbolic.

class Analyze {
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

  // sparsity pattern of L
  std::vector<int> rowsL{};
  std::vector<int> ptrL{};

  // fundamental supernodes information
  int fsn{};
  std::vector<int> fsn_belong{};
  std::vector<int> fsn_start{};
  std::vector<int> fsn_parent{};

  // relative indices of original columns wrt L columns
  std::vector<int> relind_cols{};

  // relative indices of clique wrt parent
  std::vector<std::vector<int>> relind_clique{};

  // information about consecutive indices in relind_clique
  std::vector<std::vector<int>> consecutiveSums{};

  double time_metis{};
  double time_tree{};
  double time_count{};
  double time_pattern{};
  double time_fsn{};
  double time_relind{};

  void GetPermutation();
  void Permute(const std::vector<int>& iperm);
  void ETree();
  void Postorder();
  void DFS_post(int node, int& start, std::vector<int>& head,
                const std::vector<int>& next);
  void RowColCount();
  void ColPattern();
  void FundamentalSupernodes();
  void RelativeInd_cols();
  void RelativeInd_clique();
  void Clear();
  bool Check() const;

 public:
  // Constructor: matrix must be in upper triangular format
  Analyze(const std::vector<int>& rows_input,
          const std::vector<int>& ptr_input);

  // Run analyze phase and save the result in Symbolic object S
  void Run(Symbolic& S);
};

#endif