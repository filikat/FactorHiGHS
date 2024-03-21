#include "Analyze.h"

#include <fstream>
#include <iostream>
#include <random>
#include <stack>

Analyze::Analyze(const int* row_ind, const int* col_ptr, int size,
                 int nonzeros) {
  // Input the symmetric matrix to be factorized in CSC format.
  // row_ind contains the row indices.
  // col_ptr contains the starting points of each column.
  // size is the number of rows/columns.
  // nonzeros is the number of nonzero entries.
  // Only the upper triangular part is used.

  n = size;
  rowsUpper = std::vector<int>(row_ind, row_ind + nonzeros);
  ptrUpper = std::vector<int>(col_ptr, col_ptr + n + 1);

  // Permute the matrix with identical permutation, to extract upper triangular
  // part, if the input is not upper triangular.
  std::vector<int> iperm(n);
  for (int i = 0; i < n; ++i) iperm[i] = i;
  Permute(iperm);

  nz = ptrUpper.back();
  ready = true;
}

void Analyze::GetPermutation() {
  // Use Metis to compute a nested dissection permutation of the original matrix

  perm.resize(n);
  iperm.resize(n);

  // Build temporary full copy of the matrix, to be used for Metis.
  // NB: Metis adjacency list should not contain the vertex itself, so diagonal
  // element is skipped.

  std::vector<int> work(n, 0);

  // go through the columns to count nonzeros
  for (int j = 0; j < n; ++j) {
    for (int el = ptrUpper[j]; el < ptrUpper[j + 1]; ++el) {
      int i = rowsUpper[el];

      // skip diagonal entries
      if (i == j) continue;

      // nonzero in column j
      ++work[j];

      // duplicated on the lower part of column i
      ++work[i];
    }
  }

  // compute column pointers from column counts
  std::vector<int> temp_ptr(n + 1, 0);
  Counts2Ptr(temp_ptr, work);

  std::vector<int> temp_rows(temp_ptr.back(), 0);

  for (int j = 0; j < n; ++j) {
    for (int el = ptrUpper[j]; el < ptrUpper[j + 1]; ++el) {
      int i = rowsUpper[el];

      if (i == j) continue;

      // insert row i in column j
      temp_rows[work[j]++] = i;

      // insert row j in column i
      temp_rows[work[i]++] = j;
    }
  }

  int options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);
  int status = METIS_NodeND(&n, temp_ptr.data(), temp_rows.data(), NULL,
                            options, perm.data(), iperm.data());
  assert(status == METIS_OK);
}

void Analyze::Permute(const std::vector<int>& iperm) {
  // Symmetric permutation of the upper triangular matrix based on inverse
  // permutation iperm.
  // The resulting matrix is upper triangular, regardless of the inpur matrix.

  std::vector<int> work(n, 0);

  // go through the columns to count the nonzeros
  for (int j = 0; j < n; ++j) {
    // get new index of column
    int col = iperm[j];

    // go through elements of column
    for (int el = ptrUpper[j]; el < ptrUpper[j + 1]; ++el) {
      int i = rowsUpper[el];

      // ignore potential entries in lower triangular part
      if (i > j) continue;

      // get new index of row
      int row = iperm[i];

      // since only upper triangular part is used, column is larger than row
      int actual_col = std::max(row, col);
      ++work[actual_col];
    }
  }

  std::vector<int> new_ptr(n + 1);

  // get column pointers by summing the count of nonzeros in each column.
  // copy column pointers into work
  Counts2Ptr(new_ptr, work);

  std::vector<int> new_rows(new_ptr.back());

  // go through the columns to assign row indices
  for (int j = 0; j < n; ++j) {
    // get new index of column
    int col = iperm[j];

    // go through elements of column
    for (int el = ptrUpper[j]; el < ptrUpper[j + 1]; ++el) {
      int i = rowsUpper[el];

      // ignore potential entries in lower triangular part
      if (i > j) continue;

      // get new index of row
      int row = iperm[i];

      // since only upper triangular part is used, column is larger than row
      int actual_col = std::max(row, col);
      int actual_row = std::min(row, col);

      int pos = work[actual_col]++;
      new_rows[pos] = actual_row;
    }
  }

  ptrUpper = std::move(new_ptr);
  rowsUpper = std::move(new_rows);
}

void Analyze::ETree() {
  // Find elimination tree.
  // It works only for upper triangular matrices.
  // The tree is stored in the vector parent:
  //  parent[i] = j
  // means that j is the parent of i in the tree.
  // For the root(s) of the tree, parent[root] = -1.

  parent.resize(n);
  std::vector<int> ancestor(n);
  int next{};

  for (int j = 0; j < n; ++j) {
    // initialize parent and ancestor, which are still unknown
    parent[j] = -1;
    ancestor[j] = -1;

    for (int el = ptrUpper[j]; el < ptrUpper[j + 1]; ++el) {
      for (int i = rowsUpper[el]; i != -1 && i < j; i = next) {
        // next is used to move up the tree
        next = ancestor[i];

        // ancestor keeps track of the known part of the tree, to avoid
        // repeating (aka path compression): from j there is a known path to i
        ancestor[i] = j;

        if (next == -1) parent[i] = j;
      }
    }
  }
}

void Analyze::Postorder() {
  // Find a postordering of the elimination tree using depth first search

  postorder.resize(n);

  // Create linked lists of children:
  // head[node] is the first child of node,
  // next[head[node]] is the second child,
  // next[next[head[node]]] is the third child...
  // until -1 is reached.
  std::vector<int> head(n, -1);
  std::vector<int> next(n);
  for (int node = n - 1; node >= 0; --node) {
    if (parent[node] == -1) continue;
    next[node] = head[parent[node]];
    head[parent[node]] = node;
  }

  // Execute depth first search only for root node(s)
  int start{};
  for (int node = 0; node < n; ++node) {
    if (parent[node] == -1) {
      DFS_post(node, start, head, next);
    }
  }

  // Permute elimination tree based on postorder
  std::vector<int> ipost(n);
  InversePerm(postorder, ipost);
  std::vector<int> new_parent(n);
  for (int i = 0; i < n; ++i) {
    if (parent[i] != -1) {
      new_parent[ipost[i]] = ipost[parent[i]];
    } else {
      new_parent[ipost[i]] = -1;
    }
  }
  parent = std::move(new_parent);

  // Permute matrix based on postorder
  Permute(ipost);

  // create the lower triangular part of the matrix
  rowsLower.resize(nz);
  ptrLower.resize(n + 1);
  Transpose(rowsLower, ptrLower);

  // Update perm and iperm
  PermuteVector(perm, postorder);
  InversePerm(perm, iperm);
}

void Analyze::DFS_post(int node, int& start, std::vector<int>& head,
                       const std::vector<int>& next) {
  // Perform depth first search starting from root node and order the nodes
  // starting from the value start. head and next contain the linked list of
  // children.

  std::stack<int> stack;
  stack.push(node);

  while (!stack.empty()) {
    int current = stack.top();
    int child = head[current];

    if (child == -1) {
      // no children left to order,
      // remove from the stack and order
      stack.pop();
      postorder[start++] = current;
    } else {
      // at least one child left to order,
      // add it to the stack and remove it from the list of children
      stack.push(child);
      head[current] = next[child];
    }
  }
}

void Analyze::Transpose(std::vector<int>& rowsT, std::vector<int>& ptrT) const {
  // Compute the transpose of the matrix and return it in rowsT and ptrT

  std::vector<int> work(n);

  // count the entries in each row into work
  for (int i = 0; i < ptrUpper.back(); ++i) {
    ++work[rowsUpper[i]];
  }

  // sum row sums to obtain pointers
  Counts2Ptr(ptrT, work);

  for (int j = 0; j < n; ++j) {
    for (int el = ptrUpper[j]; el < ptrUpper[j + 1]; ++el) {
      int i = rowsUpper[el];

      // entry (i,j) becomes entry (j,i)
      int pos = work[i]++;
      rowsT[pos] = j;
    }
  }
}

void Analyze::RowColCount() {
  // Compute the number of nonzero entries for each row and column of the
  // Cholesky factor.

  rowcount.resize(n);
  colcount.resize(n);

  // keep track of visited columns
  std::vector<int> mark(n, -1);

  // consider each row
  for (int i = 0; i < n; ++i) {
    // mark diagonal entry
    mark[i] = i;
    int current_row_count = 1;
    ++colcount[i];

    // for all entries in the row of lower triangle
    for (int el = ptrUpper[i]; el < ptrUpper[i + 1]; ++el) {
      int j = rowsUpper[el];
      if (j == i) continue;

      // while columns are not yet considered
      while (mark[j] != i) {
        mark[j] = i;
        ++current_row_count;
        ++colcount[j];

        // go up the elimination tree
        j = parent[j];
      }
    }

    rowcount[i] = current_row_count;
  }

  // compute nonzeros of L and number of operations
  for (int j = 0; j < n; ++j) {
    nzL += colcount[j];
    double ccj = (double)colcount[j];
    operations += ccj * ccj;
  }
}

void Analyze::ColPattern() {
  // Compute sparsity pattern of Cholesky factor L.
  // This is very similar to RowColCount(), but it needs the information about
  // the column counts.

  rowsL.resize(nzL);
  ptrL.resize(n + 1);

  // keep track of visited columns
  std::vector<int> mark(n, -1);

  // compute column pointers of L
  std::vector<int> work(colcount);
  Counts2Ptr(ptrL, work);

  // consider each row
  for (int i = 0; i < n; ++i) {
    // mark diagonal entry
    mark[i] = i;

    // there is a nonzero entry in column i at row i
    rowsL[work[i]++] = i;

    // for all entries in the row of lower triangle
    for (int el = ptrUpper[i]; el < ptrUpper[i + 1]; ++el) {
      int j = rowsUpper[el];
      if (j == i) continue;

      // while columns are not yet considered
      while (mark[j] != i) {
        mark[j] = i;

        // there is a nonzero entry in column j at row i
        rowsL[work[j]++] = i;

        // go up the elimination tree
        j = parent[j];
      }
    }
  }
}

void Analyze::FundamentalSupernodes() {
  // Find fundamental supernodes.
  // The end result is a vector fsn, such that
  //  fsn[i] = j
  // means that node i belongs to supernode j.

  // isSN[i] is true if node i is the start of a fundamental supernode
  std::vector<bool> isSN(n, false);

  std::vector<int> prev_nonz(n, -1);

  // compute sizes of subtrees
  std::vector<int> subtreeSizes(n);
  SubtreeSize(parent, subtreeSizes);

  for (int j = 0; j < n; ++j) {
    for (int el = ptrLower[j]; el < ptrLower[j + 1]; ++el) {
      int i = rowsLower[el];
      int k = prev_nonz[i];

      // mark as fundamental sn, nodes which are leaf of subtrees
      if (k < j - subtreeSizes[j] + 1) {
        isSN[j] = true;
      }

      // mark as fundamental sn, nodes which have more than one child
      if (parent[i] != -1 && subtreeSizes[i] + 1 != subtreeSizes[parent[i]]) {
        isSN[parent[i]] = true;
      }

      prev_nonz[i] = j;
    }
  }

  // create information about fundamental supernodes
  fsn_belong.resize(n);
  int sn_number = -1;
  for (int i = 0; i < n; ++i) {
    // if isSN[i] is true, then node i is the start of a new supernode
    if (isSN[i]) ++sn_number;

    // mark node i as belonging to the current supernode
    fsn_belong[i] = sn_number;
  }

  // number of supernodes found
  fsn = fsn_belong.back() + 1;

  // fsn_ptr contains pointers to the starting node of each supernode
  fsn_start.resize(fsn + 1);
  int next = 0;
  for (int i = 0; i < n; ++i) {
    if (isSN[i]) {
      fsn_start[next] = i;
      ++next;
    }
  }
  fsn_start[next] = n;

  // build supernodal elimination tree
  fsn_parent.resize(fsn);
  for (int i = 0; i < fsn - 1; ++i) {
    int j = parent[fsn_start[i + 1] - 1];
    if (j != -1) {
      fsn_parent[i] = fsn_belong[j];
    } else {
      fsn_parent[i] = -1;
    }
  }
  fsn_parent.back() = -1;
}

void Analyze::RelativeInd_cols() {
  // Find the relative indices of the original column wrt the frontal matrix of
  // the corresponding supernode

  relind_cols.resize(nz);

  // go through the supernodes
  for (int sn = 0; sn < fsn; ++sn) {
    const int ptrL_start = ptrL[fsn_start[sn]];
    const int ptrL_end = ptrL[fsn_start[sn] + 1];

    // go through the columns of the supernode
    for (int col = fsn_start[sn]; col < fsn_start[sn + 1]; ++col) {
      // go through original column and supernodal column
      int ptrA = ptrLower[col];
      int ptrL = ptrL_start;

      // offset wrt ptrLower[col]
      int index{};

      // size of the column of the original matrix
      int col_size = ptrLower[col + 1] - ptrLower[col];

      while (ptrL < ptrL_end) {
        // if found all the relative indices that are needed, stop
        if (index == col_size) {
          break;
        }

        // check if indices coincide
        if (rowsL[ptrL] == rowsLower[ptrA]) {
          // yes: save relative index and move pointers forward
          relind_cols[ptrLower[col] + index] = ptrL - ptrL_start;
          ++index;
          ++ptrL;
          ++ptrA;
        } else {
          // no: move pointer of L forward
          ++ptrL;
        }
      }
    }
  }
}

void Analyze::RelativeInd_clique() {
  // Find the relative indices of the child clique wrt the frontal matrix of the
  // parent supernode

  relind_clique.resize(fsn);

  for (int sn = 0; sn < fsn; ++sn) {
    // number of nodes in the supernode
    int sn_size = fsn_start[sn + 1] - fsn_start[sn];

    // column of the first node in the supernode
    int j = fsn_start[sn];

    // size of the first column of the supernode
    int sn_column_size = ptrL[j + 1] - ptrL[j];

    // size of the clique of the supernode
    int sn_clique_size = sn_column_size - sn_size;

    relind_clique[sn].resize(sn_clique_size);

    // iterate through the clique of sn
    int ptr_current = ptrL[fsn_start[sn]] + sn_size;

    // iterate through the full column of parent sn
    int ptr_parent = ptrL[fsn_start[fsn_parent[sn]]];

    // keep track of start and end of parent sn column
    const int ptr_parent_start = ptr_parent;
    const int ptr_parent_end = ptrL[fsn_start[fsn_parent[sn]] + 1];

    // where to write into relind
    int index{};

    // iterate though the column of the parent sn
    while (ptr_parent < ptr_parent_end) {
      // if found all the relative indices that are needed, stop
      if (index == sn_clique_size) {
        break;
      }

      // check if indices coincide
      if (rowsL[ptr_current] == rowsL[ptr_parent]) {
        // yes: save relative index and move pointers forward
        relind_clique[sn][index] = ptr_parent - ptr_parent_start;
        ++index;
        ++ptr_parent;
        ++ptr_current;
      } else {
        // no: move pointer of parent forward
        ++ptr_parent;
      }
    }
  }
}

void Analyze::Clear() {
  ready = false;
  rowsUpper.clear();
  ptrUpper.clear();
  rowsLower.clear();
  ptrLower.clear();
  n = 0;
  nz = 0;
  nzL = 0;
  operations = 0.0;
  perm.clear();
  iperm.clear();
  parent.clear();
  postorder.clear();
  colcount.clear();
  rowcount.clear();
  rowsL.clear();
  ptrL.clear();
  fsn = 0;
  fsn_belong.clear();
  fsn_start.clear();
  fsn_parent.clear();
  relind_cols.clear();
  relind_clique.clear();
}

bool Analyze::Check() const {
  // Check that the symbolic factorization is correct, by using dense linear
  // algebra operations.
  // Return true if check is successful, or if matrix is too large.
  // To be used for debug.

  int wrong_entries{};

  // Check relative indices cols
  for (int sn = 0; sn < fsn; ++sn) {
    const int fsn_col_start = ptrL[fsn_start[sn]];

    for (int col = fsn_start[sn]; col < fsn_start[sn + 1]; ++col) {
      for (int el = ptrLower[col]; el < ptrLower[col + 1]; ++el) {
        int row = rowsLower[el];

        if (row != rowsL[fsn_col_start + relind_cols[el]]) {
          printf("==> Found wrong relind_cols, col %d, row %d, rowsL %d\n", col,
                 row, rowsL[fsn_col_start + relind_cols[el]]);
          ++wrong_entries;
        }
      }
    }
  }

  // Check symbolic factorization
  if (n > 1000) {
    printf("==> Matrix is too large for dense checking\n");
    return true;
  }

  // initialize random number generator (to avoid numerical cancellation)
  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_real_distribution<double> distr(0.1, 10.0);

  // assemble sparse matrix into dense matrix
  std::vector<double> M(n * n);
  for (int col = 0; col < n; ++col) {
    for (int el = ptrUpper[col]; el < ptrUpper[col + 1]; ++el) {
      int row = rowsUpper[el];

      // insert random element in position (row,col)
      M[row + col * n] = distr(rng);

      // guarantee matrix is diagonally dominant (thus positive definite)
      if (row == col) {
        M[row + col * n] += n * 10;
      }
    }
  }

  // use Lapack to factorize the dense matrix
  char uplo = 'U';
  int N = n;
  int info;
  dpotrf(&uplo, &N, M.data(), &N, &info);
  if (info != 0) {
    printf("==> dpotrf failed\n");
    return false;
  }

  // Check that expected nonzeros are nonzeros
  for (int col = 0; col < n; ++col) {
    for (int el = ptrL[col]; el < ptrL[col + 1]; ++el) {
      int row = rowsL[el];

      // Check that nonzeros are nonzeros.
      // Numerical cancellation is highly unlikely with random data in M.
      if (M[col + row * n] == 0.0) {
        printf("==> (%d,%d) Found zero, expected nonzero\n", row, col);
        ++wrong_entries;
      }

      // set nonzeros to zero
      M[col + row * n] = 0.0;
    }
  }

  // Check that there are no nonzeros outside of expected sparsity pattern
  for (int col = 0; col < n; ++col) {
    for (int row = 0; row <= col; ++row) {
      if (M[row + col * n] != 0.0) {
        printf("==> (%d,%d) Found nonzero, expected zero\n", row, col);
        ++wrong_entries;
      }
    }
  }

  return wrong_entries == 0;
}

void Analyze::Run(Symbolic& S) {
  // Perform analyze phase and store the result into the symbolic object S.
  // After Run returns, the current Analyze object is cleared.

  if (!ready) return;

  GetPermutation();
  Permute(iperm);
  ETree();
  Postorder();
  RowColCount();
  ColPattern();
  FundamentalSupernodes();
  RelativeInd_cols();
  RelativeInd_clique();

  if (!Check()) {
    printf("\n==> Analyze check failed\n\n");
  } else {
    printf("\n==> Analyze check successful\n\n");
  }

  // move relevant stuff into S
  S.n = n;
  S.nz = nzL;
  S.operations = operations;
  S.perm = std::move(perm);
  S.iperm = std::move(iperm);
  S.parent = std::move(parent);
  S.rowcount = std::move(rowcount);
  S.colcount = std::move(colcount);
  S.rows = std::move(rowsL);
  S.ptr = std::move(ptrL);
  S.fsn = fsn;
  S.fsn_parent = std::move(fsn_parent);
  S.fsn_start = std::move(fsn_start);
  S.relind_cols = std::move(relind_cols);
  S.relind_clique = std::move(relind_clique);

  // clear the remaining data
  Clear();
}