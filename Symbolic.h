#ifndef SYMBOLIC_H
#define SYMBOLIC_H

#include <vector>

// Type of factorization:
// normal equations or augmented system
enum FactType { NormEq, AugSys };
enum PackType { Full, HybridPacked, HybridHybrid };

class Symbolic {
  // Type of factorization
  FactType type_{};

  // Packed or full format
  PackType pack_format_ = PackType::HybridPacked;

  // Size of blocks for dense factorization
  const int block_size_ = 128;

  // Size of the matrix L
  int n_{};

  // Number of nonzeros in L
  double nz_{};
  double fillin_{};
  double max_storage_{};

  // Number of dense operations and assembly operations
  double dense_ops_{};
  double assembly_ops_{};

  // Number of supernodes
  int sn_{};

  // Number of artificial nonzero entries introduced to merge supernodes
  int artificial_nz_{};
  double artificial_ops_{};

  // size of the largest frontal matrix and largest sn
  int largest_front_{};
  int largest_sn_{};

  // Permutation and inverse permutation
  std::vector<int> perm_{};
  std::vector<int> iperm_{};

  // Sparsity pattern of each supernode of L
  std::vector<int> rows_{};
  std::vector<int> ptr_{};

  // Supernodal elimination tree:
  // - sn_parent_[i] gives the parent of supernode i in the supernodal
  //   elimination tree
  std::vector<int> sn_parent_{};

  // Supernode initial node:
  // - sn_start_[i] gives the first node in supernode i.
  //   Supernode i is made of nodes from sn_start_[i] to sn_start_[i+1]-1
  std::vector<int> sn_start_{};

  // Relative indices of original columns wrt columns of L.
  // - relind_cols_[i] contains the relative indices of entry i, with respect to
  //   the numbering of the frontal matrix of the corresponding supernode.
  // - Given the row indices of the original matrix, rowsA:
  //   relind_cols_[i] = k implies that the i-th entry of the original matrix
  //   (which has original row index given by rowsA[i]) corresponds to the row
  //   in position k in the frontal matrix of the supernode corresponding to the
  //   column to which the i-th entry belongs.
  //   This is useful when assemblying the entries of the original matrix into
  //   the frontal matrix.
  std::vector<int> relind_cols_{};

  // Relative indices of clique wrt parent supernode.
  // - relind_clique_[i] contains the local indices of the nonzero rows of the
  //   clique of the current supernode with respect to the numbering of the
  //   parent supernode.
  // - relind_clique_[i][j] = k implies that the row in position j in the clique
  //   of supernode i corresponds to the row in position k in the frontal matrix
  //   of supernode sn_parent_[i].
  //   This is useful when summing the generated elements from supernode i into
  //   supernode sn_parent_[i].
  std::vector<std::vector<int>> relind_clique_{};

  // Number of consecutive sums that can be done with one BLAS call.
  // - consecutive_sums_[i] contains information about the assembly of supernode
  //   i into the frontal matrix of its parent.
  // - consecutive_sums_[i][j] = k implies that, when summing contributions from
  //   row j of the clique of supernodes i into the frontal matrix of its
  //   parent, k consecutive indices are found. This means that instead of doing
  //   k individual sums, we can use one single call to daxpy, with k entries
  //   and increment equal to one.
  std::vector<std::vector<int>> consecutive_sums_{};

  friend class Analyse;

 public:
  // print information to screen
  void print() const;

  // provide const access to symbolic factorization
  FactType type() const;
  PackType packFormat() const;
  int blockSize() const;
  int size() const;
  int nz() const;
  double ops() const;
  double assemblyOps() const;
  int sn() const;
  int rows(int i) const;
  int ptr(int i) const;
  int snStart(int i) const;
  int relindCols(int i) const;
  int relindClique(int i, int j) const;
  int consecutiveSums(int i, int j) const;
  const std::vector<int>& ptr() const;
  const std::vector<int>& perm() const;
  const std::vector<int>& iperm() const;
  const std::vector<int>& snParent() const;
  const std::vector<int>& snStart() const;
};

// Explanation of relative indices:
// Each supernode i corresponds to a frontal matrix Fi.
// The indices of the rows of Fi are called Ri.
// Ri contains the indices of the supernode
//  {sn_start_[i],...,sn_start_[i+1]-1}
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
// The original matrix has the following structure instead
//
//        2 3 4
//
//  2     x
//  3     x x
//  4     x 0 x
// ...
//  7     x 0 0
// ...
// 15     x x 0
//
// The parent of supernode i is snParent[i] = p.
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
// The original matrix, for columns 2,3,4, has indices
//  {2,3,4,7,15,3,15,4}.
// relind_cols_, for entries corresponding to columns 2,3,4, has the relative
// position of these indices wrt the indices in Ri {2,3,4,7,15}, i.e.,
// {0,1,2,3,4,1,4,2}.
//
// relind_clique_[i] contains the relative position of the indices of the clique
// of supernode i {7,15} with respect to Rp {7,8,9,14,15,17,19}, i.e.,
// relind_clique_[i] = {0,4}.

// Explanation of consecutive sums:
// if relind_clique_[i] = {2,5,8,9,10,11,12,14}, there are (up to) 8 entries
// that need to be summed for each column of the clique. However, 5 of these
// indices are consecutive {8,9,10,11,12}. Summing these consecutive entries can
// be done using daxpy with increment equal to one, which is more efficient that
// summing one by one. consecutive_sums_[i] would contain {1,1,5,4,3,2,1,1},
// which means that, if we start from a given row, we can find out how many
// consecutive copies can be done. E.g., starting from row 4,
// consecutive_sums_[i][4] = 3, which means that the next 3 indices need not be
// summed by hand, but they can be done using daxpy.

#endif