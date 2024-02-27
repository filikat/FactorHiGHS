#ifndef ANALYZE_H
#define ANALYZE_H

#include <vector>

#include "GKlib.h"
#include "metis.h"
#include "Aux_analyze.h"

class Analyze {
 public:
  // Information about the original matrix.
  const int* original_rows{};
  const int* original_ptr{};
  int original_nz{};
  bool original_upper = false;

  // Matrix to be factorized, stored in upper triangular format
  std::vector<int> rows{};
  std::vector<int> ptr{};
  int n{};
  int nz{};

  // Permutation and inverse permutation from Metis
  std::vector<int> perm{};
  std::vector<int> iperm{};

  // Elimination tree
  std::vector<int> parent{};

  std::vector<int> postorder{};

 public:
  Analyze(const int* row_ind, const int* col_ptr, int size, int nonzeros,
          bool is_upper);
  void GetPermutation();
  void Permute(const std::vector<int>& iperm);
  void ETree();
  void Postorder();
  void DFS_post(int node, int& start, std::vector<int>& head,
                const std::vector<int>& next);
};

#endif