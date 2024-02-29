#ifndef SYMBOLIC_H
#define SYMBOLIC_H

#include <vector>

class Symbolic {
 public:
  // general information
  int n{};
  int nz{};
  int fsn{};

  // permutation
  std::vector<int> perm{};
  std::vector<int> iperm{};

  // elimination tree (postordered)
  std::vector<int> parent{};

  // row and column counts
  std::vector<int> rowcount{};
  std::vector<int> colcount{};

  // sparsity pattern of L
  std::vector<int> rows{};
  std::vector<int> ptr{};

  // supernodal elimination tree
  std::vector<int> fsn_parent{};

  // supernodal pointers
  std::vector<int> fsn_ptr{};

  // relative indices for frontal matrices
  std::vector<std::vector<int>> relind{};

  friend class Analyze;

 public:
  int sn_begin(int sn) const;
  int sn_end(int sn) const;
  int clique_begin(int sn) const;
  int clique_end(int sn) const;

  void clique_info(int sn, int& position, int& snsize, int& cliquesize) const;
};

#endif