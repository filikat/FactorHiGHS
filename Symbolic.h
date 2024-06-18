#ifndef SYMBOLIC_H
#define SYMBOLIC_H

#include <vector>

// Type of factorization:
// normal equations or augmented system
enum class FactType { NormEq, AugSys };
enum class PackType { Full, Hybrid, Hybrid2 };

class Symbolic {
  // Type of factorization
  FactType type{};

  // Packed or full format
  PackType packed = PackType::Hybrid;

  // Size of blocks for dense factorization
  const int blockSize = 128;

  // Size of the matrix L
  int n{};

  // Number of nonzeros in L
  double nz{};
  double fillin{};
  double maxStorage{};

  // Number of dense operations and assembly operations
  double operations{};
  double assemblyOp{};

  // Number of supernodes
  int sn{};

  // Number of artificial nonzero entries introduced to merge supernodes
  int artificialNz{};
  double artificialOp{};

  // size of the largest frontal matrix and largest sn
  int largestFront{};
  int largestSn{};

  // Permutation and inverse permutation
  std::vector<int> perm{};
  std::vector<int> iperm{};

  // Sparsity pattern of each supernode of L
  std::vector<int> rows{};
  std::vector<int> ptr{};

  // Supernodal elimination tree:
  // - snParent[i] gives the parent of supernode i in the supernodal
  //   elimination tree
  std::vector<int> snParent{};

  // Supernode initial node:
  // - snStart[i] gives the first node in supernode i.
  //   Supernode i is made of nodes from snStart[i] to snStart[i+1]-1
  std::vector<int> snStart{};

  // Relative indices of original columns wrt columns of L.
  // - relindCols[i] contains the relative indices of entry i, with respect to
  //   the numbering of the frontal matrix of the corresponding supernode.
  // - Given the row indices of the original matrix, rowsA:
  //   relindCols[i] = k implies that the i-th entry of the original matrix
  //   (which has original row index given by rowsA[i]) corresponds to the row
  //   in position k in the frontal matrix of the supernode corresponding to the
  //   column to which the i-th entry belongs.
  //   This is useful when assemblying the entries of the original matrix into
  //   the frontal matrix.
  std::vector<int> relindCols{};

  // Relative indices of clique wrt parent supernode.
  // - relindClique[i] contains the local indices of the nonzero rows of the
  //   clique of the current supernode with respect to the numbering of the
  //   parent supernode.
  // - relindClique[i][j] = k implies that the row in position j in the clique
  //   of supernode i corresponds to the row in position k in the frontal matrix
  //   of supernode snParent[i].
  //   This is useful when summing the generated elements from supernode i into
  //   supernode snParent[i].
  std::vector<std::vector<int>> relindClique{};

  // Number of consecutive sums that can be done with one BLAS call.
  // - consecutiveSums[i] contains information about the assembly of supernode i
  //   into the frontal matrix of its parent.
  // - consecutiveSums[i][j] = k implies that, when summing contributions from
  //   row j of the clique of supernodes i into the frontal matrix of its
  //   parent, k consecutive indices are found. This means that instead of doing
  //   k individual sums, we can use one single call to daxpy, with k entries
  //   and increment equal to one.
  std::vector<std::vector<int>> consecutiveSums{};

  friend class Analyse;

 public:
  // print information to screen
  void Print() const;

  // provide const access to symbolic factorization
  FactType Type() const;
  PackType Packed() const;
  int BlockSize() const;
  int Size() const;
  int Nz() const;
  double Ops() const;
  double AssemblyOps() const;
  int Sn() const;
  int Rows(int i) const;
  int Ptr(int i) const;
  int SnStart(int i) const;
  int RelindCols(int i) const;
  int RelindClique(int i, int j) const;
  int ConsecutiveSums(int i, int j) const;
  const std::vector<int>& Ptr() const;
  const std::vector<int>& Perm() const;
  const std::vector<int>& Iperm() const;
  const std::vector<int>& SnParent() const;
  const std::vector<int>& SnStart() const;
};

// Explanation of relative indices:
// Each supernode i corresponds to a frontal matrix Fi.
// The indices of the rows of Fi are called Ri.
// Ri contains the indices of the supernode
//  {snStart[i],...,snStart[i+1]-1}
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
// relindCols, for entries corresponding to columns 2,3,4, has the relative
// position of these indices wrt the indices in Ri {2,3,4,7,15}, i.e.,
// {0,1,2,3,4,1,4,2}.
//
// relindClique[i] contains the relative position of the indices of the clique
// of supernode i {7,15} with respect to Rp {7,8,9,14,15,17,19}, i.e.,
// relindClique[i] = {0,4}.

// Explanation of consecutive sums:
// if relindClique[i] = {2,5,8,9,10,11,12,14}, there are (up to) 8 entries that
// need to be summed for each column of the clique.
// However, 5 of these indices are consecutive {8,9,10,11,12}. Summing these
// consecutive entries can be done using daxpy with increment equal to one,
// which is more efficient that summing one by one.
// consecutiveSums[i] would contain {1,1,5,4,3,2,1,1}, which means that, if we
// start from a given row, we can find out how many consecutive copies can be
// done.
// E.g., starting from row 4, consecutiveSums[i][4] = 3, which means that the
// next 3 indices need not be summed by hand, but they can be done using daxpy.

#endif