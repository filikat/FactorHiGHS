#include "Analyze.h"

#include <iostream>
#include <stack>

Analyze::Analyze(const int* row_ind, const int* col_ptr, int size, int nonzeros,
                 bool is_upper) {
  // Input the symmetric matrix to be factorized in CSC format.
  // row_ind contains the row indices.
  // col_ptr contains the starting points of each column.
  // size is the number of rows/columns.
  // nonzeros is the number of nonzero entries.
  // is_upper is true if only upper triangular part is stored,
  //          is false if full matrix is stored

  // save pointers to original matrix
  original_ptr = col_ptr;
  original_rows = row_ind;
  original_upper = is_upper;
  original_nz = nonzeros;

  n = size;
  rows = std::vector<int>(row_ind, row_ind + nonzeros);
  ptr = std::vector<int>(col_ptr, col_ptr + n + 1);

  if (!original_upper) {
    // If original matrix is full, extract triangular part.
    // This can be done using Permute, with the identical permutation.

    std::vector<int> iperm(n);
    for (int i = 0; i < n; ++i) iperm[i] = i;

    Permute(iperm);
  }

  nz = ptr.back();
  ready = true;
}

void Analyze::GetPermutation() {
  // Use Metis to compute a nested dissection permutation of the original matrix

  perm.resize(n);
  iperm.resize(n);

  if (!original_upper) {
    // If original matrix is full, save a temporary local copy, because Metis
    // takes non-const pointers.

    std::vector<int> temp_ptr(original_ptr, original_ptr + n + 1);
    std::vector<int> temp_rows(original_rows, original_rows + original_nz);
    int status = METIS_NodeND(&n, temp_ptr.data(), temp_rows.data(), NULL, NULL,
                              perm.data(), iperm.data());
    assert(status == METIS_OK);

  } else {
    // If original matrix is upper triangular, build temporary full copy to be
    // used for Metis.

    std::vector<int> work(n, 0);

    // go through the columns to count nonzeros
    for (int j = 0; j < n; ++j) {
      for (int el = ptr[j]; el < ptr[j + 1]; ++el) {
        int i = rows[el];
        ++work[j];

        // if entry (i,j) is not on the diagonal, it is duplicated on the lower
        // part of column i
        if (i < j) {
          ++work[i];
        }
      }
    }

    // compute column pointers from column counts
    std::vector<int> temp_ptr(n + 1, 0);
    Counts2Ptr(temp_ptr, work);

    std::vector<int> temp_rows(temp_ptr.back(), 0);

    for (int j = 0; j < n; ++j) {
      for (int el = ptr[j]; el < ptr[j + 1]; ++el) {
        int i = rows[el];

        // insert row i in column j
        temp_rows[work[j]++] = i;

        if (i < j) {
          // insert row j in column i
          temp_rows[work[i]++] = j;
        }
      }
    }

    int status = METIS_NodeND(&n, temp_ptr.data(), temp_rows.data(), NULL, NULL,
                              perm.data(), iperm.data());
    assert(status == METIS_OK);
  }
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
    for (int el = ptr[j]; el < ptr[j + 1]; ++el) {
      int i = rows[el];

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
    for (int el = ptr[j]; el < ptr[j + 1]; ++el) {
      int i = rows[el];

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

  ptr = std::move(new_ptr);
  rows = std::move(new_rows);
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

    for (int el = ptr[j]; el < ptr[j + 1]; ++el) {
      for (int i = rows[el]; i != -1 && i < j; i = next) {
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
  // next[head[node]] is the second chilf,
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
  for (int i = 0; i < ptr.back(); ++i) {
    ++work[rows[i]];
  }

  // sum row sums to obtain pointers
  Counts2Ptr(ptrT, work);

  for (int j = 0; j < n; ++j) {
    for (int el = ptr[j]; el < ptr[j + 1]; ++el) {
      int i = rows[el];

      // entry (i,j) becomes entry (j,i)
      int pos = work[i]++;
      rowsT[pos] = j;
    }
  }
}

void Analyze::RowColCount() {
  // Compute the number of nonzero entries for each row and column of the
  // Cholesky factor.
  // It requires O(|L|) operations.
  // The method that uses O(|A|) operations does not work for now.
  // I will come back to it later.

  rowcount.resize(n);
  colcount.resize(n);

  // keep track of visited columns
  std::vector<int> mark(n, -1);

  // consider each row
  for (int i = 0; i < n; ++i) {
    // mark diagonal entry
    mark[i] = i;
    int rowc = 1;
    ++colcount[i];

    // for all entries in the row of lower triangle
    for (int el = ptr[i]; el < ptr[i + 1]; ++el) {
      int j = rows[el];
      if (j == i) continue;

      // while columns are not yet considered
      while (mark[j] != i) {
        mark[j] = i;
        ++rowc;
        ++colcount[j];

        // go up the elimination tree
        j = parent[j];
      }
    }

    rowcount[i] = rowc;
    nzL += rowc;
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
    for (int el = ptr[i]; el < ptr[i + 1]; ++el) {
      int j = rows[el];
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

  std::vector<int> prev_nonz(n, 0);

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

  // guarantee that node 0 is true (for some reason it does not happen)
  isSN[0] = true;

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
  fsn_ptr.resize(fsn + 1);
  int next = 0;
  for (int i = 0; i < n; ++i) {
    if (isSN[i]) {
      fsn_ptr[next] = i;
      ++next;
    }
  }
  fsn_ptr[next] = n;

  // build supernodal elimination tree
  fsn_parent.resize(fsn);
  for (int i = 0; i < fsn - 1; ++i) {
    int j = parent[fsn_ptr[i + 1] - 1];
    if (j != -1) {
      fsn_parent[i] = fsn_belong[j];
    } else {
      fsn_parent[i] = -1;
    }
  }
  fsn_parent.back() = -1;
}

void Analyze::RelativeInd() {
  // Find the relative indices of the child clique with respect to the parent
  // indexing.
  // E.g., supernode 0 has parent supernode 5.
  // Supernode 0 is made of 1 node and its column is {0,8,9}, i.e., it has
  // nonzero entries in its clique in rows {8,9}.
  // Parent supernode is made of nodes {7,8,9} and the full column has nonzeros
  // in rows {7,8,9,15}.
  // Then, the relative indices of 0 with respect to its parent 5 are {1,2},
  // because entry 8 in the clique of sn 0 corresponds to the entry number 1 in
  // the full column of the parent sn 5...

  relind.resize(fsn);

  for (int sn = 0; sn < fsn; ++sn) {
    // number of nodes in the supernode
    int sn_size = fsn_ptr[sn + 1] - fsn_ptr[sn];

    // column of the first node in the supernode
    int j = fsn_ptr[sn];

    // size of the first column of the supernode
    int sn_column_size = ptrL[j + 1] - ptrL[j];

    // size of the clique of the supernode
    int sn_clique_size = sn_column_size - sn_size;

    relind[sn].resize(sn_clique_size);

    // iterate through the clique of sn
    int ptr_current = ptrL[fsn_ptr[sn]] + sn_size;

    // iterate through the full column of parent sn
    int ptr_parent = ptrL[fsn_ptr[fsn_parent[sn]]];

    // keep track of start and end of parent sn column
    const int ptr_parent_start = ptr_parent;
    const int ptr_parent_end = ptrL[fsn_ptr[fsn_parent[sn]] + 1];

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
        relind[sn][index] = ptr_parent - ptr_parent_start;
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

void Analyze::clear() {
  original_rows = NULL;
  original_ptr = NULL;
  original_nz = 0;
  original_upper = false;
  ready = false;
  rows.clear();
  ptr.clear();
  rowsLower.clear();
  ptrLower.clear();
  n = 0;
  nz = 0;
  nzL = 0;
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
  fsn_ptr.clear();
  fsn_parent.clear();
  relind.clear();
}

void Analyze::Run(Symbolic& S, bool runSymbolic, bool runSupernodes) {
  // Perform analyze phase and store the result into the symbolic object S.
  // After Run returns, the current Analyze object is cleared.

  if (!ready) return;

  GetPermutation();
  Permute(iperm);
  ETree();
  Postorder();
  RowColCount();

  // potentially perform full symbolic factorization
  if (runSymbolic) {
    ColPattern();
  }

  // potentially look for supernodes
  if (runSupernodes) {
    FundamentalSupernodes();
    RelativeInd();
  }

  // move relevant stuff into S
  S.n = n;
  S.nz = nzL;
  S.perm = std::move(perm);
  S.iperm = std::move(iperm);
  S.parent = std::move(parent);
  S.rowcount = std::move(rowcount);
  S.colcount = std::move(colcount);
  S.rows = std::move(rowsL);
  S.ptr = std::move(ptrL);
  S.fsn = fsn;
  S.fsn_parent = std::move(fsn_parent);
  S.fsn_ptr = std::move(fsn_ptr);
  S.relind = std::move(relind);

  // clear the remaining data
  clear();
}