#ifndef SYMBOLIC_H
#define SYMBOLIC_H

#include <vector>

class Symbolic {
 public:
  // Size of the matrix L
  int n{};

  // Number of nonzeros in L
  int nz{};

  // Number of fundamental supernodes
  int fsn{};

  // Permutation and inverse permutation
  std::vector<int> perm{};
  std::vector<int> iperm{};

  // Elimination tree (postordered):
  // - parent[i] gives the parent of node i in the elimination tree
  std::vector<int> parent{};

  // Row and column counts:
  // - rowcount[i] gives the number of nonzero entries in row i of L
  // - colcount[i] gives the number of nonzero entries in column i of L
  std::vector<int> rowcount{};
  std::vector<int> colcount{};

  // Sparsity pattern of L in CSC format
  std::vector<int> rows{};
  std::vector<int> ptr{};

  // Supernodal elimination tree:
  // - fsn_parent[i] gives the parent of supernode i in the supernodal
  //   elimination tree
  std::vector<int> fsn_parent{};

  // Supernode initial node:
  // - fsn_start[i] gives the first node in supernode i.
  //   Supernode i is made of nodes from fsn_start[i] to fsn_start[i+1]-1
  std::vector<int> fsn_start{};

  // Relative indices for frontal matrices:
  // - relind[i] contains the local indices of the nonzero rows of the clique of
  //   the current supernode with respect to the numbering of the parent
  //   supernode.
  // - relind[i][j] = k implies that the row in position j in the clique of
  //   supernode i corresponds to the row in position k in the frontal matrix
  //   of supernode fsn_parent[i].
  //   This is useful when summing the generated elements from supernode i into
  //   supernode fsn_parent[i].
  std::vector<std::vector<int>> relind{};

  friend class Analyze;

 public:
  int sn_begin(int sn) const;
  int sn_end(int sn) const;
  int clique_begin(int sn) const;
  int clique_end(int sn) const;

  void clique_info(int sn, int& position, int& snsize, int& cliquesize) const;
};

// Each supernode i corresponds to a frontal matrix Fi.
// The indices of the rows of Fi are called Ri.
// Ri contains the indices of the supernode
//  {fsn_start[i],...,fsn_start[i+1]-1}
// and then the indices of the clique, or generated element
// (i.e., the entries of the Schur complement that are modified).
//
// E.g., supernode i has the following structure:
//
//        2 3 4
//
//  2     x
//  3     x x
//  4     x x x
// ...
//  7     x x x
// ...
// 15     x x x
//
// The supernode is made of nodes {2,3,4}.
// The clique is made of indices {7,15}.
// The frontal matrix Fi has 5 rows which correspond to the rows of L given by
//  the indices Ri = {2,3,4,7,15}.
//
// The parent of supernode i is fsn_parent[i] = p.
// Supernode p has the following structure:
//
//        7 8 9
//
// 7      x
// 8      x x
// 9      x x x
// ...
// 14     x x x
// 15     x x x
// 16
// 17     x x x
// 18
// 19     x x x
//
// The supernode is made of nodes {7,8,9}.
// The clique is made of indices {14,15,17,19}.
// The frontal matrix Fp has 7 rows which correspond to the rows of L given by
//  the indices Rp = {7,8,9,14,15,17,19}.
//
// relind[i] contains the relative position of the indices of the clique of
//  supernode i {7,15} with respect to Rp {7,8,9,14,15,17,19},
//  i.e., relind[i] = {0,4}.

#endif