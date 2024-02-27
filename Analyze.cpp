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
    ColCount2Ptr(temp_ptr, work);

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
  ColCount2Ptr(new_ptr, work);

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
  Inverse_perm(postorder, ipost);
  std::vector<int> new_parent(n);
  for (int i = 0; i < n; ++i) {
    new_parent[ipost[i]] = ipost[parent[i]];
  }
  parent = std::move(new_parent);

  // Permute matrix based on postorder
  Permute(ipost);

  // Update perm and iperm
  Permute_vector(perm, postorder);
  Inverse_perm(perm, iperm);
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

