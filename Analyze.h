#ifndef ANALYZE_H
#define ANALYZE_H

#include <vector>

#include "Aux_analyze.h"
#include "GKlib.h"
#include "Symbolic.h"
#include "metis.h"

// Class to perform the analyze phase of the factorization.
// The final symbolic factorization is stored in an object of type Symbolic.

class Analyze {
  // Information about the original matrix.
  const int* original_rows{};
  const int* original_ptr{};
  int original_nz{};
  bool original_upper = false;
  bool ready = false;

  // Matrix to be factorized, stored in upper and lower triangular format
  std::vector<int> rows{};
  std::vector<int> ptr{};
  std::vector<int> rowsLower{};
  std::vector<int> ptrLower{};
  int n{};
  int nz{};
  int nzL{};

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
  std::vector<int> fsn_ptr{};
  std::vector<int> fsn_parent{};

  // relative indices for frontal matrices
  std::vector<std::vector<int>> relind{};

  void GetPermutation();
  void Permute(const std::vector<int>& iperm);
  void ETree();
  void Postorder();
  void DFS_post(int node, int& start, std::vector<int>& head,
                const std::vector<int>& next);
  void Transpose(std::vector<int>& rowsT, std::vector<int>& ptrT) const;
  void RowColCount();
  void ColPattern();
  void FundamentalSupernodes();
  void RelativeInd();
  void clear();

 public:
  Analyze(const int* row_ind, const int* col_ptr, int size, int nonzeros,
          bool is_upper);
  void Run(Symbolic& S, bool runSymbolic, bool runSupernodes);
};

#endif