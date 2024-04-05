#include "Analyze.h"

#include <fstream>
#include <iostream>
#include <random>
#include <stack>

Analyze::Analyze(const std::vector<int>& rows_input,
                 const std::vector<int>& ptr_input) {
  // Input the symmetric matrix to be analyzed in CSC format.
  // row_ind contains the row indices.
  // col_ptr contains the starting points of each column.
  // size is the number of rows/columns.
  // nonzeros is the number of nonzero entries.
  // Only the lower triangular part is used.

  n = ptr_input.size() - 1;
  nz = rows_input.size();

  // Create upper triangular part
  rowsUpper.resize(nz);
  ptrUpper.resize(n + 1);
  Transpose(ptr_input, rows_input, ptrUpper, rowsUpper);

  // Permute the matrix with identical permutation, to extract upper triangular
  // part, if the input is not upper triangular.
  std::vector<int> iperm(n);
  for (int i = 0; i < n; ++i) iperm[i] = i;
  Permute(iperm);

  // actual number of nonzeros of only upper triangular part
  nz = ptrUpper.back();

  // number of nonzeros potentially changed after Permute.
  rowsUpper.resize(nz);

  // double transpose to sort columns
  ptrLower.resize(n + 1);
  rowsLower.resize(nz);
  Transpose(ptrUpper, rowsUpper, ptrLower, rowsLower);
  Transpose(ptrLower, rowsLower, ptrUpper, rowsUpper);

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
      const int i = rowsUpper[el];

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
      const int i = rowsUpper[el];

      if (i == j) continue;

      // insert row i in column j
      temp_rows[work[j]++] = i;

      // insert row j in column i
      temp_rows[work[i]++] = j;
    }
  }

  // call Metis
  int options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);
  int status = METIS_NodeND(&n, temp_ptr.data(), temp_rows.data(), NULL,
                            options, perm.data(), iperm.data());
  assert(status == METIS_OK);
}

void Analyze::Permute(const std::vector<int>& iperm) {
  // Symmetric permutation of the upper triangular matrix based on inverse
  // permutation iperm.
  // The resulting matrix is upper triangular, regardless of the input matrix.

  std::vector<int> work(n, 0);

  // go through the columns to count the nonzeros
  for (int j = 0; j < n; ++j) {
    // get new index of column
    const int col = iperm[j];

    // go through elements of column
    for (int el = ptrUpper[j]; el < ptrUpper[j + 1]; ++el) {
      const int i = rowsUpper[el];

      // ignore potential entries in lower triangular part
      if (i > j) continue;

      // get new index of row
      const int row = iperm[i];

      // since only upper triangular part is used, col is larger than row
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
    const int col = iperm[j];

    // go through elements of column
    for (int el = ptrUpper[j]; el < ptrUpper[j + 1]; ++el) {
      const int i = rowsUpper[el];

      // ignore potential entries in lower triangular part
      if (i > j) continue;

      // get new index of row
      const int row = iperm[i];

      // since only upper triangular part is used, column is larger than row
      const int actual_col = std::max(row, col);
      const int actual_row = std::min(row, col);

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

  // create linked list of children
  std::vector<int> head, next;
  ChildrenLinkedList(parent, head, next);

  // Execute depth first search only for root node(s)
  int start{};
  for (int node = 0; node < n; ++node) {
    if (parent[node] == -1) {
      DFS_post(node, start, head, next, postorder);
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

  // double transpose to sort columns and compute lower part
  Transpose(ptrUpper, rowsUpper, ptrLower, rowsLower);
  Transpose(ptrLower, rowsLower, ptrUpper, rowsUpper);

  // Update perm and iperm
  PermuteVector(perm, postorder);
  InversePerm(perm, iperm);
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

  // compute nonzeros of L
  for (int j = 0; j < n; ++j) {
    nzL += colcount[j];
  }
}

void Analyze::FundamentalSupernodes() {
  // Find fundamental supernodes.

  // isSN[i] is true if node i is the start of a fundamental supernode
  std::vector<bool> isSN(n, false);

  std::vector<int> prev_nonz(n, -1);

  // compute sizes of subtrees
  std::vector<int> subtreeSizes(n);
  SubtreeSize(parent, subtreeSizes);

  for (int j = 0; j < n; ++j) {
    for (int el = ptrLower[j]; el < ptrLower[j + 1]; ++el) {
      const int i = rowsLower[el];
      const int k = prev_nonz[i];

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
  sn_belong.resize(n);
  int sn_number = -1;
  for (int i = 0; i < n; ++i) {
    // if isSN[i] is true, then node i is the start of a new supernode
    if (isSN[i]) ++sn_number;

    // mark node i as belonging to the current supernode
    sn_belong[i] = sn_number;
  }

  // number of supernodes found
  sn_count = sn_belong.back() + 1;

  // fsn_ptr contains pointers to the starting node of each supernode
  sn_start.resize(sn_count + 1);
  int next = 0;
  for (int i = 0; i < n; ++i) {
    if (isSN[i]) {
      sn_start[next] = i;
      ++next;
    }
  }
  sn_start[next] = n;

  // build supernodal elimination tree
  sn_parent.resize(sn_count);
  for (int i = 0; i < sn_count - 1; ++i) {
    int j = parent[sn_start[i + 1] - 1];
    if (j != -1) {
      sn_parent[i] = sn_belong[j];
    } else {
      sn_parent[i] = -1;
    }
  }
  sn_parent.back() = -1;
}

void Analyze::RelaxSupernodes() {
  // =================================================
  // build information about supernodes
  // =================================================
  std::vector<int> sn_size(sn_count);
  std::vector<int> clique_size(sn_count);
  std::vector<int> fake_nonzeros(sn_count);
  for (int i = 0; i < sn_count; ++i) {
    sn_size[i] = sn_start[i + 1] - sn_start[i];
    clique_size[i] = colcount[sn_start[i]] - sn_size[i];
    fake_nonzeros[i] = 0;
  }

  // build linked lists of children
  std::vector<int> firstChild, nextChild;
  ChildrenLinkedList(sn_parent, firstChild, nextChild);

  // =================================================
  // Merge supernodes
  // =================================================
  std::vector<int> mergedInto(sn_count, -1);
  int merged_sn{};

  for (int sn = 0; sn < sn_count; ++sn) {
    // keep iterating through the children of the supernode, until there's no
    // more child to merge with

    std::vector<int> merged_sons{};

    while (true) {
      int child = firstChild[sn];

      int smallest_fake_nonzeros = INT_MAX;
      int size_smallest = 0;
      int child_smallest = -1;

      while (child != -1) {
        // how many zero rows would become nonzero
        int rowsFilled = sn_size[sn] + clique_size[sn] - clique_size[child];

        // how many zero entries would become nonzero
        int nzAdded = rowsFilled * sn_size[child];

        // how many artificial nonzeros would the merged supernode have
        int totalArtNz = nzAdded + fake_nonzeros[sn] + fake_nonzeros[child];

        // Save child with smallest number of artificial zeros created.
        // Ties are broken based on size of child.
        if (totalArtNz < smallest_fake_nonzeros ||
            (totalArtNz == smallest_fake_nonzeros &&
             size_smallest < sn_size[child])) {
          smallest_fake_nonzeros = totalArtNz;
          size_smallest = sn_size[child];
          child_smallest = child;
        }
        child = nextChild[child];
      }

      if (smallest_fake_nonzeros <= max_artificial_nz) {
        // if we found a child to merge with the parent

        // update information of parent
        sn_size[sn] += size_smallest;
        fake_nonzeros[sn] = smallest_fake_nonzeros;

        // count number of merges supernodes
        ++merged_sn;

        // save information about merging of supernodes
        mergedInto[child_smallest] = sn;
        merged_sons.push_back(child_smallest);

        // remove child from linked list of children
        child = firstChild[sn];
        if (child == child_smallest) {
          // child_smallest is the first child
          firstChild[sn] = nextChild[child_smallest];
        } else {
          while (nextChild[child] != child_smallest) {
            child = nextChild[child];
          }
          // now child is the previous child of child_smallest
          nextChild[child] = nextChild[child_smallest];
        }

      } else {
        // no more children can be merged with parent
        break;
      }
    }
  }

  // merging of supernodes completed

  // number of new supernodes
  int new_sn_count = sn_count - merged_sn;

  // keep track of number of row indices needed for each supernode
  sn_indices.resize(new_sn_count);

  // =================================================
  // Create supernodal permutation
  // =================================================

  // permutation of supernodes needed after merging
  std::vector<int> sn_perm(sn_count);

  // number of new sn that includes the old sn
  std::vector<int> new_id(sn_count);

  // new sn pointer vector
  std::vector<int> new_sn_start(new_sn_count + 1);

  {
    // keep track of the children merged into a given supernode
    std::vector<std::vector<int>> receivedFrom(sn_count, std::vector<int>());

    // index to write into sn_perm
    int start_perm{};

    // index to write into new_sn_start
    int sn_start_ind{};

    // next available number for new sn numbering
    int next_id{};

    for (int sn = 0; sn < sn_count; ++sn) {
      if (mergedInto[sn] > -1) {
        // Current sn was merged into its parent.
        // Save information about which supernode sn was merged into
        receivedFrom[mergedInto[sn]].push_back(sn);
      } else {
        // Current sn was not merged into its parent.
        // It is one of the new sn.

        // Add merged supernodes to the permutation, recursively.

        ++sn_start_ind;

        std::stack<int> toadd;
        toadd.push(sn);

        while (!toadd.empty()) {
          int current = toadd.top();

          if (!receivedFrom[current].empty()) {
            for (int i : receivedFrom[current]) toadd.push(i);
            receivedFrom[current].clear();
          } else {
            toadd.pop();
            sn_perm[start_perm++] = current;
            new_id[current] = next_id;

            // count number of nodes in each new supernode
            new_sn_start[sn_start_ind] +=
                sn_start[current + 1] - sn_start[current];

            // sum colcount for each column in the new sn
            for (int col = sn_start[current]; col < sn_start[current + 1];
                 ++col) {
              sn_indices[next_id] += colcount[col];
            }
          }
        }

        // keep track of total number of artificial nonzeros
        artificialNz += fake_nonzeros[sn];

        // add artificial nonzeros to the sum of colcounts
        sn_indices[next_id] += fake_nonzeros[sn];

        // Compute number of indices for new sn with the formula
        // [ \sum_i colcount_i + fake_nz + k*(k-1)/2 ] / k
        // where k is the number of columns in the new sn.
        int k = new_sn_start[sn_start_ind];
        sn_indices[next_id] += k * (k - 1) / 2;
        sn_indices[next_id] /= k;

        ++next_id;
      }
    }

    // new_sn_start contain the number of cols in each new sn.
    // sum them to obtain the sn pointers.
    for (int i = 0; i < new_sn_count; ++i) {
      new_sn_start[i + 1] += new_sn_start[i];
    }
  }

  // include artificial nonzeros in the nonzeros of the factor
  nzL += artificialNz;

  // compute number of flops needed for the factorization
  for (int sn = 0; sn < new_sn_count; ++sn) {
    double colccount_sn = (double)sn_indices[sn];
    for (int i = 0; i < new_sn_start[sn + 1] - new_sn_start[sn]; ++i) {
      operations += (colccount_sn - i) * (colccount_sn - i);
    }
  }

  // =================================================
  // Create nodal permutation
  // =================================================
  // Given the supernodal permutation, find the nodal permutation needed after
  // sn merging.

  // permutation to apply to the existing one
  std::vector<int> new_perm(n);

  // index to write into new_perm
  int start{};

  for (int i = 0; i < sn_count; ++i) {
    int sn = sn_perm[i];
    for (int j = sn_start[sn]; j < sn_start[sn + 1]; ++j) {
      new_perm[start++] = j;
    }
  }

  // obtain inverse permutation
  std::vector<int> new_iperm(n);
  InversePerm(new_perm, new_iperm);

  // =================================================
  // Create new sn elimination tree
  // =================================================
  std::vector<int> new_sn_parent(new_sn_count, -1);
  for (int i = 0; i < sn_count; ++i) {
    if (sn_parent[i] == -1) continue;

    int ii = new_id[i];
    int pp = new_id[sn_parent[i]];

    if (ii == pp) continue;

    new_sn_parent[ii] = pp;
  }

  // =================================================
  // Save new information
  // =================================================

  // build new sn_belong, i.e., the sn to which each columb belongs
  for (int sn = 0; sn < sn_count; ++sn) {
    for (int i = sn_start[sn]; i < sn_start[sn + 1]; ++i) {
      sn_belong[i] = new_id[sn];
    }
  }
  PermuteVector(sn_belong, new_perm);

  // Overwrite previous data
  sn_parent = std::move(new_sn_parent);
  sn_start = std::move(new_sn_start);
  sn_count = new_sn_count;

  // Permute matrix based on new permutation
  Permute(new_iperm);

  // double transpose to sort columns and compute lower part
  Transpose(ptrUpper, rowsUpper, ptrLower, rowsLower);
  Transpose(ptrLower, rowsLower, ptrUpper, rowsUpper);

  // Update perm and iperm
  PermuteVector(perm, new_perm);
  InversePerm(perm, iperm);
}

void Analyze::SnPattern() {
  // number of total indices needed
  int indices{};

  for (int i : sn_indices) indices += i;

  // allocate space for sn pattern
  rowsLsn.resize(indices);
  ptrLsn.resize(sn_count + 1);

  // keep track of visited supernodes
  std::vector<int> mark(sn_count, -1);

  // compute column pointers of L
  std::vector<int> work(sn_indices);
  Counts2Ptr(ptrLsn, work);

  // consider each row
  for (int i = 0; i < n; ++i) {
    // for all entries in the row of lower triangle
    for (int el = ptrUpper[i]; el < ptrUpper[i + 1]; ++el) {
      // there is nonzero (i,j)
      int j = rowsUpper[el];

      // supernode to which column j belongs to
      int snj = sn_belong[j];

      // while supernodes are not yet considered
      while (snj != -1 && mark[snj] != i) {
        // we may end up too far
        if (sn_start[snj] > i) break;

        // supernode snj is now considered for row i
        mark[snj] = i;

        // there is a nonzero entry in supernode snj at row i
        rowsLsn[work[snj]++] = i;

        // go up the elimination tree
        snj = sn_parent[snj];
      }
    }
  }
}

void Analyze::RelativeInd_cols() {
  // Find the relative indices of the original column wrt the frontal matrix of
  // the corresponding supernode

  relind_cols.resize(nz);

  // go through the supernodes
  for (int sn = 0; sn < sn_count; ++sn) {
    const int ptL_start = ptrLsn[sn];
    const int ptL_end = ptrLsn[sn + 1];

    // go through the columns of the supernode
    for (int col = sn_start[sn]; col < sn_start[sn + 1]; ++col) {
      // go through original column and supernodal column
      int ptA = ptrLower[col];
      int ptL = ptL_start;

      // offset wrt ptrLower[col]
      int index{};

      // size of the column of the original matrix
      int col_size = ptrLower[col + 1] - ptrLower[col];

      while (ptL < ptL_end) {
        // if found all the relative indices that are needed, stop
        if (index == col_size) {
          break;
        }

        // check if indices coincide
        if (rowsLsn[ptL] == rowsLower[ptA]) {
          // yes: save relative index and move pointers forward
          relind_cols[ptrLower[col] + index] = ptL - ptL_start;
          ++index;
          ++ptL;
          ++ptA;
        } else {
          // no: move pointer of L forward
          ++ptL;
        }
      }
    }
  }
}

void Analyze::RelativeInd_clique() {
  // Find the relative indices of the child clique wrt the frontal matrix of the
  // parent supernode

  relind_clique.resize(sn_count);
  consecutiveSums.resize(sn_count);

  for (int sn = 0; sn < sn_count; ++sn) {
    // if there is no parent, skip supernode
    if (sn_parent[sn] == -1) continue;

    // number of nodes in the supernode
    const int sn_size = sn_start[sn + 1] - sn_start[sn];

    // column of the first node in the supernode
    const int j = sn_start[sn];

    // size of the first column of the supernode
    const int sn_column_size = ptrLsn[sn + 1] - ptrLsn[sn];

    // size of the clique of the supernode
    const int sn_clique_size = sn_column_size - sn_size;

    relind_clique[sn].resize(sn_clique_size);

    // iterate through the clique of sn
    int ptr_current = ptrLsn[sn] + sn_size;

    // iterate through the full column of parent sn
    int ptr_parent = ptrLsn[sn_parent[sn]];

    // keep track of start and end of parent sn column
    const int ptr_parent_start = ptr_parent;
    const int ptr_parent_end = ptrLsn[sn_parent[sn] + 1];

    // where to write into relind
    int index{};

    // iterate though the column of the parent sn
    while (ptr_parent < ptr_parent_end) {
      // if found all the relative indices that are needed, stop
      if (index == sn_clique_size) {
        break;
      }

      // check if indices coincide
      if (rowsLsn[ptr_current] == rowsLsn[ptr_parent]) {
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

    // Difference between consecutive relative indices.
    // Useful to detect chains of consecutive indices.
    consecutiveSums[sn].resize(sn_clique_size);
    for (int i = 0; i < sn_clique_size - 1; ++i) {
      consecutiveSums[sn][i] = relind_clique[sn][i + 1] - relind_clique[sn][i];
    }

    // Number of consecutive sums that can be done in one blas call.
    consecutiveSums[sn].back() = 1;
    for (int i = sn_clique_size - 2; i >= 0; --i) {
      if (consecutiveSums[sn][i] > 1) {
        consecutiveSums[sn][i] = 1;
      } else if (consecutiveSums[sn][i] == 1) {
        consecutiveSums[sn][i] = consecutiveSums[sn][i + 1] + 1;
      } else {
        printf("Error in consecutiveSums\n");
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
  // rowsL.clear();
  // ptrL.clear();
  sn_count = 0;
  sn_belong.clear();
  sn_start.clear();
  sn_parent.clear();
  relind_cols.clear();
  relind_clique.clear();
  consecutiveSums.clear();
}

bool Analyze::Check() const {
  // Check that the symbolic factorization is correct, by using dense linear
  // algebra operations.
  // Return true if check is successful, or if matrix is too large.
  // To be used for debug.

  // Check symbolic factorization
  if (n > 5000) {
    printf("\n==> Matrix is too large for dense checking\n\n");
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
    printf("\n==> dpotrf failed\n\n");
    return false;
  }

  // assemble expected sparsity pattern into dense matrix
  std::vector<bool> L(n * n);
  for (int sn = 0; sn < sn_count; ++sn) {
    for (int col = sn_start[sn]; col < sn_start[sn + 1]; ++col) {
      for (int el = ptrLsn[sn]; el < ptrLsn[sn + 1]; ++el) {
        int row = rowsLsn[el];
        if (row < col) continue;
        L[row + n * col] = true;
      }
    }
  }

  int zeros_found{};
  int wrong_entries{};

  // check how many entries do not correspond
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < n; ++i) {
      double valM = M[j + n * i];
      bool valL = L[i + n * j];

      if (valL && valM == 0.0) {
        // count number of fake zeros found, to confront it with artificialNz
        ++zeros_found;
      } else if (!valL && valM != 0.0) {
        printf("==> (%d,%d) Found nonzero, expected zero\n", i, j);
        ++wrong_entries;
      }
    }
  }

  if (wrong_entries == 0 && zeros_found == artificialNz) {
    printf("\n==> Analyze check successful\n\n");
    return true;
  } else {
    printf("\n==> Analyze check failed\n\n");
    return false;
  }
}

void Analyze::PrintTimes() const {
  printf("\n----------------------------------------------------\n");
  printf("\t\tAnalyze\n");
  printf("----------------------------------------------------\n");
  printf("Analyze time            \t%f\n", time_total);
  printf("\tMetis:                  %f (%4.1f%%)\n", time_metis,
         time_metis / time_total * 100);
  printf("\tTree:                   %f (%4.1f%%)\n", time_tree,
         time_tree / time_total * 100);
  printf("\tCounts:                 %f (%4.1f%%)\n", time_count,
         time_count / time_total * 100);
  printf("\tSupernodes:             %f (%4.1f%%)\n", time_sn,
         time_sn / time_total * 100);
  printf("\tSn sparsity pattern:    %f (%4.1f%%)\n", time_pattern,
         time_pattern / time_total * 100);
  printf("\tRelative indices:       %f (%4.1f%%)\n", time_relind,
         time_relind / time_total * 100);
}

void Analyze::Run(Symbolic& S) {
  // Perform analyze phase and store the result into the symbolic object S.
  // After Run returns, the current Analyze object is cleared.

  if (!ready) return;

  Clock clock0{};
  clock0.start();

  Clock clock{};

  clock.start();
  GetPermutation();
  time_metis = clock.stop();

  clock.start();
  Permute(iperm);
  ETree();
  Postorder();
  time_tree = clock.stop();

  clock.start();
  RowColCount();
  time_count = clock.stop();

  clock.start();
  FundamentalSupernodes();
  RelaxSupernodes();
  time_sn = clock.stop();

  clock.start();
  SnPattern();
  time_pattern = clock.stop();

  clock.start();
  RelativeInd_cols();
  RelativeInd_clique();
  time_relind = clock.stop();

  time_total = clock0.stop();

  PrintTimes();

  Check();

  // move relevant stuff into S
  S.n = n;
  S.nz = nzL;
  S.sn = sn_count;
  S.artificialNz = artificialNz;
  S.largestFront = *std::max_element(sn_indices.begin(), sn_indices.end());
  std::vector<int> temp(sn_start);
  for (int i = sn_count; i > 0; --i) temp[i] -= temp[i - 1];
  S.largestSn = *std::max_element(temp.begin(), temp.end());

  S.operations = operations;
  S.perm = std::move(perm);
  S.iperm = std::move(iperm);
  S.rows = std::move(rowsLsn);
  S.ptr = std::move(ptrLsn);
  S.sn_parent = std::move(sn_parent);
  S.sn_start = std::move(sn_start);
  S.relind_cols = std::move(relind_cols);
  S.relind_clique = std::move(relind_clique);
  S.consecutiveSums = std::move(consecutiveSums);

  // clear the remaining data
  Clear();
}
