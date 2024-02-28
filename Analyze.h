#ifndef ANALYZE_H
#define ANALYZE_H

#include <vector>

#include "Aux_analyze.h"
#include "GKlib.h"
#include "metis.h"

class Analyze {
 public:
  // Information about the original matrix.
  const int* original_rows{};
  const int* original_ptr{};
  int original_nz{};
  bool original_upper = false;

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

 public:
  Analyze(const int* row_ind, const int* col_ptr, int size, int nonzeros,
          bool is_upper);
  void GetPermutation();
  void Permute(const std::vector<int>& iperm);
  void ETree();
  void Postorder();
  void DFS_post(int node, int& start, std::vector<int>& head,
                const std::vector<int>& next);
  void Transpose(std::vector<int>& rowsT, std::vector<int>& ptrT) const;
  void RowColCount();
  void ColPattern(std::vector<int>& rowsL, std::vector<int>& ptrL) const;
};

#endif