#include "Analyse.h"

#include <fstream>
#include <iostream>
#include <random>
#include <stack>

Analyse::Analyse(const std::vector<int>& rows_input,
                 const std::vector<int>& ptr_input, FactType type_input,
                 const std::vector<int>& order) {
  // Input the symmetric matrix to be analysed in CSC format.
  // row_ind contains the row indices.
  // col_ptr contains the starting points of each column.
  // size is the number of rows/columns.
  // nonzeros is the number of nonzero entries.
  // Only the lower triangular part is used.

  n = ptr_input.size() - 1;
  nz = rows_input.size();
  type = type_input;

  // Create upper triangular part
  rowsUpper.resize(nz);
  ptrUpper.resize(n + 1);
  Transpose(ptr_input, rows_input, ptrUpper, rowsUpper);

  // Permute the matrix with identical permutation, to extract upper triangular
  // part, if the input is not upper triangular.
  std::vector<int> id_perm(n);
  for (int i = 0; i < n; ++i) id_perm[i] = i;
  Permute(id_perm);

  // actual number of nonzeros of only upper triangular part
  nz = ptrUpper.back();

  // number of nonzeros potentially changed after Permute.
  rowsUpper.resize(nz);

  // double transpose to sort columns
  ptrLower.resize(n + 1);
  rowsLower.resize(nz);
  Transpose(ptrUpper, rowsUpper, ptrLower, rowsLower);
  Transpose(ptrLower, rowsLower, ptrUpper, rowsUpper);

  if (!order.empty()) {
    // inverse permutation provided by user
    iperm = order;
    perm.resize(n);
    InversePerm(iperm, perm);
  }

  ready = true;
}

void Analyse::GetPermutation() {
  // Use Metis to compute a nested dissection permutation of the original matrix

  if (!perm.empty()) {
    // permutation already provided by user
    return;
  }

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

  metis_order = iperm;
}

void Analyse::Permute(const std::vector<int>& iperm) {
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

void Analyse::ETree() {
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

void Analyse::Postorder() {
  // Find a postordering of the elimination tree using depth first search

  postorder.resize(n);

  // create linked list of children
  std::vector<int> head, next;
  ChildrenLinkedList(parent, head, next);

  // Execute depth first search only for root node(s)
  int start{};
  for (int node = 0; node < n; ++node) {
    if (parent[node] == -1) {
      Dfs_post(node, start, head, next, postorder);
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

void Analyse::ColCount() {
  // Columns count using skeleton matrix.
  // Taken from Tim Davis "Direct Methods for Sparse Linear Systems".

  std::vector<int> first(n, -1);
  std::vector<int> ancestor(n, -1);
  std::vector<int> max_first(n, -1);
  std::vector<int> prev_leaf(n, -1);

  colCount.resize(n);

  // find first descendant
  for (int k = 0; k < n; ++k) {
    int j = k;
    colCount[j] = (first[j] == -1) ? 1 : 0;
    while (j != -1 && first[j] == -1) {
      first[j] = k;
      j = parent[j];
    }
  }

  // each node belongs to a separate set
  for (int j = 0; j < n; j++) ancestor[j] = j;

  for (int k = 0; k < n; ++k) {
    int j = k;

    // if not a root, decrement
    if (parent[j] != -1) colCount[parent[j]]--;

    // process edges of matrix
    for (int el = ptrLower[j]; el < ptrLower[j + 1]; ++el) {
      ProcessEdge(j, rowsLower[el], first, max_first, colCount, prev_leaf,
                  ancestor);
    }

    if (parent[j] != -1) ancestor[j] = parent[j];
  }

  // sum contributions from each child
  for (int j = 0; j < n; ++j) {
    if (parent[j] != -1) {
      colCount[parent[j]] += colCount[j];
    }
  }

  // compute nonzeros of L
  operationsNorelax = 0.0;
  nzL = 0;
  for (int j = 0; j < n; ++j) {
    nzL += (double)colCount[j];
    operationsNorelax += (double)(colCount[j] - 1) * (colCount[j] - 1);
  }
}

void Analyse::FundamentalSupernodes() {
  // Find fundamental supernodes.

  // isSN[i] is true if node i is the start of a fundamental supernode
  std::vector<bool> is_sn(n, false);

  std::vector<int> prev_nonz(n, -1);

  // compute sizes of subtrees
  std::vector<int> subtree_sizes(n);
  SubtreeSize(parent, subtree_sizes);

  for (int j = 0; j < n; ++j) {
    for (int el = ptrLower[j]; el < ptrLower[j + 1]; ++el) {
      const int i = rowsLower[el];
      const int k = prev_nonz[i];

      // mark as fundamental sn, nodes which are leaf of subtrees
      if (k < j - subtree_sizes[j] + 1) {
        is_sn[j] = true;
      }

      // mark as fundamental sn, nodes which have more than one child
      if (parent[i] != -1 && subtree_sizes[i] + 1 != subtree_sizes[parent[i]]) {
        is_sn[parent[i]] = true;
      }

      prev_nonz[i] = j;
    }
  }

  // create information about fundamental supernodes
  snBelong.resize(n);
  int sn_number = -1;
  for (int i = 0; i < n; ++i) {
    // if isSN[i] is true, then node i is the start of a new supernode
    if (is_sn[i]) ++sn_number;

    // mark node i as belonging to the current supernode
    snBelong[i] = sn_number;
  }

  // number of supernodes found
  snCount = snBelong.back() + 1;

  // fsn_ptr contains pointers to the starting node of each supernode
  snStart.resize(snCount + 1);
  int next = 0;
  for (int i = 0; i < n; ++i) {
    if (is_sn[i]) {
      snStart[next] = i;
      ++next;
    }
  }
  snStart[next] = n;

  // build supernodal elimination tree
  snParent.resize(snCount);
  for (int i = 0; i < snCount - 1; ++i) {
    int j = parent[snStart[i + 1] - 1];
    if (j != -1) {
      snParent[i] = snBelong[j];
    } else {
      snParent[i] = -1;
    }
  }
  snParent.back() = -1;
}

void Analyse::RelaxSupernodes() {
  // Child which produces smallest number of fake nonzeros is merged if
  // resulting sn has fewer than max_artificial_nz fake nonzeros.
  // Multiple values of max_artificial_nz are tried, chosen with bisection
  // method, until the percentage of artificial nonzeros is in the range [1,2]%.

  int max_artificial_nz = k_start_thresh_relax;
  int largest_below = -1;
  int smallest_above = -1;

  for (int iter = 0; iter < k_max_iter_relax; ++iter) {
    // =================================================
    // build information about supernodes
    // =================================================
    std::vector<int> sn_size(snCount);
    std::vector<int> clique_size(snCount);
    fakeNonzeros.assign(snCount, 0);
    for (int i = 0; i < snCount; ++i) {
      sn_size[i] = snStart[i + 1] - snStart[i];
      clique_size[i] = colCount[snStart[i]] - sn_size[i];
      fakeNonzeros[i] = 0;
    }

    // build linked lists of children
    std::vector<int> first_child, next_child;
    ChildrenLinkedList(snParent, first_child, next_child);

    // =================================================
    // Merge supernodes
    // =================================================
    mergedInto.assign(snCount, -1);
    mergedSn = 0;

    for (int sn = 0; sn < snCount; ++sn) {
      // keep iterating through the children of the supernode, until there's no
      // more child to merge with

      while (true) {
        int child = first_child[sn];

        // info for first criterion
        int nz_fakenz = INT_MAX;
        int size_fakenz = 0;
        int child_fakenz = -1;

        while (child != -1) {
          // how many zero rows would become nonzero
          int rows_filled = sn_size[sn] + clique_size[sn] - clique_size[child];

          // how many zero entries would become nonzero
          int nz_added = rows_filled * sn_size[child];

          // how many artificial nonzeros would the merged supernode have
          int total_art_nz = nz_added + fakeNonzeros[sn] + fakeNonzeros[child];

          // Save child with smallest number of artificial zeros created.
          // Ties are broken based on size of child.
          if (total_art_nz < nz_fakenz ||
              (total_art_nz == nz_fakenz && size_fakenz < sn_size[child])) {
            nz_fakenz = total_art_nz;
            size_fakenz = sn_size[child];
            child_fakenz = child;
          }

          child = next_child[child];
        }

        if (nz_fakenz <= max_artificial_nz) {
          // merging creates fewer nonzeros than the maximum allowed

          // update information of parent
          sn_size[sn] += size_fakenz;
          fakeNonzeros[sn] = nz_fakenz;

          // count number of merged supernodes
          ++mergedSn;

          // save information about merging of supernodes
          mergedInto[child_fakenz] = sn;

          // remove child from linked list of children
          child = first_child[sn];
          if (child == child_fakenz) {
            // child_smallest is the first child
            first_child[sn] = next_child[child_fakenz];
          } else {
            while (next_child[child] != child_fakenz) {
              child = next_child[child];
            }
            // now child is the previous child of child_smallest
            next_child[child] = next_child[child_fakenz];
          }

        } else {
          // no more children can be merged with parent
          break;
        }
      }
    }

    // compute total number of artificial nonzeros and artificial ops for this
    // value of max_artificial_nz
    double temp_art_nz{};
    double temp_art_ops{};
    for (int sn = 0; sn < snCount; ++sn) {
      if (mergedInto[sn] == -1) {
        temp_art_nz += fakeNonzeros[sn];

        double nn = sn_size[sn];
        double cc = clique_size[sn];
        temp_art_ops += (nn + cc) * (nn + cc) * nn - (nn + cc) * nn * (nn + 1) +
                        nn * (nn + 1) * (2 * nn + 1) / 6;
      }
    }
    temp_art_ops -= operationsNorelax;

    // if enough fake nz or ops have been added, stop.
    // double ratio_fake = temp_art_nz / (nzL + temp_art_nz);
    double ratio_fake = temp_art_ops / (temp_art_ops + operationsNorelax);

    // try to find ratio in interval [0.01,0.02] using bisection
    if (ratio_fake < k_lower_ratio_relax) {
      // ratio too small
      largest_below = max_artificial_nz;
      if (smallest_above == -1) {
        max_artificial_nz *= 2;
      } else {
        max_artificial_nz = (largest_below + smallest_above) / 2;
      }
    } else if (ratio_fake > k_upper_ratio_relax) {
      // ratio too large
      smallest_above = max_artificial_nz;
      if (largest_below == -1) {
        max_artificial_nz /= 2;
      } else {
        max_artificial_nz = (largest_below + smallest_above) / 2;
      }
    } else {
      // good ratio
      return;
    }
  }
}

void Analyse::RelaxSupernodes_2() {
  // Smallest child is merged with parent, if child is small enough.

  // =================================================
  // build information about supernodes
  // =================================================
  std::vector<int> sn_size(snCount);
  std::vector<int> clique_size(snCount);
  fakeNonzeros.assign(snCount, 0);
  for (int i = 0; i < snCount; ++i) {
    sn_size[i] = snStart[i + 1] - snStart[i];
    clique_size[i] = colCount[snStart[i]] - sn_size[i];
    fakeNonzeros[i] = 0;
  }

  // build linked lists of children
  std::vector<int> first_child, next_child;
  ChildrenLinkedList(snParent, first_child, next_child);

  // =================================================
  // Merge supernodes
  // =================================================
  mergedInto.assign(snCount, -1);
  mergedSn = 0;

  for (int sn = 0; sn < snCount; ++sn) {
    // keep iterating through the children of the supernode, until there's no
    // more child to merge with

    while (true) {
      int child = first_child[sn];

      // info for first criterion
      int size_smallest = INT_MAX;
      int child_smallest = -1;
      int nz_smallest = 0;

      while (child != -1) {
        // how many zero rows would become nonzero
        int rows_filled = sn_size[sn] + clique_size[sn] - clique_size[child];

        // how many zero entries would become nonzero
        int nz_added = rows_filled * sn_size[child];

        // how many artificial nonzeros would the merged supernode have
        int total_art_nz = nz_added + fakeNonzeros[sn] + fakeNonzeros[child];

        if (sn_size[child] < size_smallest) {
          size_smallest = sn_size[child];
          child_smallest = child;
          nz_smallest = total_art_nz;
        }

        child = next_child[child];
      }

      if (size_smallest < 8 && sn_size[sn] < 8) {
        // smallest supernode is small enough to be merged with parent

        // update information of parent
        sn_size[sn] += size_smallest;
        fakeNonzeros[sn] = nz_smallest;

        // count number of merged supernodes
        ++mergedSn;

        // save information about merging of supernodes
        mergedInto[child_smallest] = sn;

        // remove child from linked list of children
        child = first_child[sn];
        if (child == child_smallest) {
          // child_smallest is the first child
          first_child[sn] = next_child[child_smallest];
        } else {
          while (next_child[child] != child_smallest) {
            child = next_child[child];
          }
          // now child is the previous child of child_smallest
          next_child[child] = next_child[child_smallest];
        }

      } else {
        // no more children can be merged with parent
        break;
      }
    }
  }
}

void Analyse::AfterRelaxSn() {
  // number of new supernodes
  int new_snCount = snCount - mergedSn;

  // keep track of number of row indices needed for each supernode
  snIndices.assign(new_snCount, 0);

  // =================================================
  // Create supernodal permutation
  // =================================================

  // permutation of supernodes needed after merging
  std::vector<int> sn_perm(snCount);

  // number of new sn that includes the old sn
  std::vector<int> new_id(snCount);

  // new sn pointer vector
  std::vector<int> new_snStart(new_snCount + 1);

  // keep track of the children merged into a given supernode
  std::vector<std::vector<int>> received_from(snCount, std::vector<int>());

  // index to write into sn_perm
  int start_perm{};

  // index to write into new_snStart
  int snStart_ind{};

  // next available number for new sn numbering
  int next_id{};

  for (int sn = 0; sn < snCount; ++sn) {
    if (mergedInto[sn] > -1) {
      // Current sn was merged into its parent.
      // Save information about which supernode sn was merged into
      received_from[mergedInto[sn]].push_back(sn);
    } else {
      // Current sn was not merged into its parent.
      // It is one of the new sn.

      // Add merged supernodes to the permutation, recursively.

      ++snStart_ind;

      std::stack<int> toadd;
      toadd.push(sn);

      while (!toadd.empty()) {
        int current = toadd.top();

        if (!received_from[current].empty()) {
          for (int i : received_from[current]) toadd.push(i);
          received_from[current].clear();
        } else {
          toadd.pop();
          sn_perm[start_perm++] = current;
          new_id[current] = next_id;

          // count number of nodes in each new supernode
          new_snStart[snStart_ind] += snStart[current + 1] - snStart[current];
        }
      }

      // keep track of total number of artificial nonzeros
      artificialNz += fakeNonzeros[sn];

      // Compure number of indices for new sn.
      // This is equal to the number of columns in the new sn plus the clique
      // size of the original supernode where the children where merged.
      snIndices[next_id] = new_snStart[snStart_ind] + colCount[snStart[sn]] -
                           snStart[sn + 1] + snStart[sn];

      ++next_id;
    }
  }

  // new_snStart contain the number of cols in each new sn.
  // sum them to obtain the sn pointers.
  for (int i = 0; i < new_snCount; ++i) {
    new_snStart[i + 1] += new_snStart[i];
  }

  // include artificial nonzeros in the nonzeros of the factor
  nzL += (double)artificialNz;

  // compute number of flops needed for the factorization
  operations = 0.0;
  for (int sn = 0; sn < new_snCount; ++sn) {
    double colcount_sn = (double)snIndices[sn];
    for (int i = 0; i < new_snStart[sn + 1] - new_snStart[sn]; ++i) {
      operations += (colcount_sn - i - 1) * (colcount_sn - i - 1);
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

  for (int i = 0; i < snCount; ++i) {
    int sn = sn_perm[i];
    for (int j = snStart[sn]; j < snStart[sn + 1]; ++j) {
      new_perm[start++] = j;
    }
  }

  // obtain inverse permutation
  std::vector<int> new_iperm(n);
  InversePerm(new_perm, new_iperm);

  // =================================================
  // Create new sn elimination tree
  // =================================================
  std::vector<int> new_snParent(new_snCount, -1);
  for (int i = 0; i < snCount; ++i) {
    if (snParent[i] == -1) continue;

    int ii = new_id[i];
    int pp = new_id[snParent[i]];

    if (ii == pp) continue;

    new_snParent[ii] = pp;
  }

  // =================================================
  // Save new information
  // =================================================

  // build new snBelong, i.e., the sn to which each columb belongs
  for (int sn = 0; sn < snCount; ++sn) {
    for (int i = snStart[sn]; i < snStart[sn + 1]; ++i) {
      snBelong[i] = new_id[sn];
    }
  }
  PermuteVector(snBelong, new_perm);

  // Overwrite previous data
  snParent = std::move(new_snParent);
  snStart = std::move(new_snStart);
  snCount = new_snCount;

  // Permute matrix based on new permutation
  Permute(new_iperm);

  // double transpose to sort columns and compute lower part
  Transpose(ptrUpper, rowsUpper, ptrLower, rowsLower);
  Transpose(ptrLower, rowsLower, ptrUpper, rowsUpper);

  // Update perm and iperm
  PermuteVector(perm, new_perm);
  InversePerm(perm, iperm);
}

void Analyse::SnPattern() {
  // number of total indices needed
  int indices{};

  for (int i : snIndices) indices += i;

  // allocate space for sn pattern
  rowsLsn.resize(indices);
  ptrLsn.resize(snCount + 1);

  // keep track of visited supernodes
  std::vector<int> mark(snCount, -1);

  // compute column pointers of L
  std::vector<int> work(snIndices);
  Counts2Ptr(ptrLsn, work);

  // consider each row
  for (int i = 0; i < n; ++i) {
    // for all entries in the row of lower triangle
    for (int el = ptrUpper[i]; el < ptrUpper[i + 1]; ++el) {
      // there is nonzero (i,j)
      int j = rowsUpper[el];

      // supernode to which column j belongs to
      int snj = snBelong[j];

      // while supernodes are not yet considered
      while (snj != -1 && mark[snj] != i) {
        // we may end up too far
        if (snStart[snj] > i) break;

        // supernode snj is now considered for row i
        mark[snj] = i;

        // there is a nonzero entry in supernode snj at row i
        rowsLsn[work[snj]++] = i;

        // go up the elimination tree
        snj = snParent[snj];
      }
    }
  }
}

void Analyse::RelativeIndCols() {
  // Find the relative indices of the original column wrt the frontal matrix of
  // the corresponding supernode

  relindCols.resize(nz);

  // go through the supernodes
  for (int sn = 0; sn < snCount; ++sn) {
    const int ptL_start = ptrLsn[sn];
    const int ptL_end = ptrLsn[sn + 1];

    // go through the columns of the supernode
    for (int col = snStart[sn]; col < snStart[sn + 1]; ++col) {
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
          relindCols[ptrLower[col] + index] = ptL - ptL_start;
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

void Analyse::RelativeIndClique() {
  // Find the relative indices of the child clique wrt the frontal matrix of the
  // parent supernode

  relindClique.resize(snCount);
  consecutiveSums.resize(snCount);

  for (int sn = 0; sn < snCount; ++sn) {
    // if there is no parent, skip supernode
    if (snParent[sn] == -1) continue;

    // number of nodes in the supernode
    const int sn_size = snStart[sn + 1] - snStart[sn];

    // column of the first node in the supernode
    const int j = snStart[sn];

    // size of the first column of the supernode
    const int sn_column_size = ptrLsn[sn + 1] - ptrLsn[sn];

    // size of the clique of the supernode
    const int sn_clique_size = sn_column_size - sn_size;

    // count number of assembly operations during factorize
    operationsAssembly += sn_clique_size * (sn_clique_size + 1) / 2;

    relindClique[sn].resize(sn_clique_size);

    // iterate through the clique of sn
    int ptr_current = ptrLsn[sn] + sn_size;

    // iterate through the full column of parent sn
    int ptr_parent = ptrLsn[snParent[sn]];

    // keep track of start and end of parent sn column
    const int ptr_parent_start = ptr_parent;
    const int ptr_parent_end = ptrLsn[snParent[sn] + 1];

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
        relindClique[sn][index] = ptr_parent - ptr_parent_start;
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
      consecutiveSums[sn][i] = relindClique[sn][i + 1] - relindClique[sn][i];
    }

    // Number of consecutive sums that can be done in one blas call.
    consecutiveSums[sn].back() = 1;
    for (int i = sn_clique_size - 2; i >= 0; --i) {
      if (consecutiveSums[sn][i] > 1) {
        consecutiveSums[sn][i] = 1;
      } else if (consecutiveSums[sn][i] == 1) {
        consecutiveSums[sn][i] = consecutiveSums[sn][i + 1] + 1;
      } else {
        printf("Error in consecutiveSums %d\n", consecutiveSums[sn][i]);
      }
    }
  }
}

bool Analyse::Check() const {
  // Check that the symbolic factorization is correct, by using dense linear
  // algebra operations.
  // Return true if check is successful, or if matrix is too large.
  // To be used for debug.

  // Check symbolic factorization
  if (n > 5000) {
    printf("\n==> Matrix is too large for dense check\n\n");
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
  for (int sn = 0; sn < snCount; ++sn) {
    for (int col = snStart[sn]; col < snStart[sn + 1]; ++col) {
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
      double val_M = M[j + n * i];
      bool val_L = L[i + n * j];

      if (val_L && val_M == 0.0) {
        // count number of fake zeros found, to confront it with artificialNz
        ++zeros_found;
      } else if (!val_L && val_M != 0.0) {
        printf("==> (%d,%d) Found nonzero, expected zero\n", i, j);
        ++wrong_entries;
      }
    }
  }

  if (wrong_entries == 0 && zeros_found == artificialNz) {
    printf("\n==> Analyse check successful\n\n");
    return true;
  } else {
    printf("\n==> Analyse check failed\n\n");
    return false;
  }
}

void Analyse::PrintTimes() const {
  printf("\n----------------------------------------------------\n");
  printf("\t\tAnalyse\n");
  printf("----------------------------------------------------\n");
  printf("\nAnalyse time            \t%8.4f\n", time_total);
  printf("\tMetis:                  %8.4f (%4.1f%%)\n", time_metis,
         time_metis / time_total * 100);
  printf("\tTree:                   %8.4f (%4.1f%%)\n", time_tree,
         time_tree / time_total * 100);
  printf("\tCounts:                 %8.4f (%4.1f%%)\n", time_count,
         time_count / time_total * 100);
  printf("\tSupernodes:             %8.4f (%4.1f%%)\n", time_sn,
         time_sn / time_total * 100);
  printf("\tSn sparsity pattern:    %8.4f (%4.1f%%)\n", time_pattern,
         time_pattern / time_total * 100);
  printf("\tRelative indices:       %8.4f (%4.1f%%)\n", time_relind,
         time_relind / time_total * 100);
}

void Analyse::Run(Symbolic& S) {
  // Perform analyse phase and store the result into the symbolic object S.
  // After Run returns, the Analyse object is not valid.

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
  ColCount();
  time_count = clock.stop();

  clock.start();
  FundamentalSupernodes();
  RelaxSupernodes();
  AfterRelaxSn();
  time_sn = clock.stop();

  clock.start();
  SnPattern();
  time_pattern = clock.stop();

  clock.start();
  RelativeIndCols();
  RelativeIndClique();
  time_relind = clock.stop();

  time_total = clock0.stop();

  PrintTimes();

  // move relevant stuff into S
  S.type = type;
  S.n = n;
  S.nz = nzL;
  S.fillin = (double)nzL / nz;
  S.sn = snCount;
  S.artificialNz = artificialNz;
  S.artificialOp = (double)operations - operationsNorelax;
  S.assemblyOp = operationsAssembly;
  S.largestFront = *std::max_element(snIndices.begin(), snIndices.end());

  std::vector<int> temp(snStart);
  for (int i = snCount; i > 0; --i) temp[i] -= temp[i - 1];
  S.largestSn = *std::max_element(temp.begin(), temp.end());

  S.operations = operations;
  S.perm = std::move(perm);
  S.iperm = std::move(iperm);
  S.rows = std::move(rowsLsn);
  S.ptr = std::move(ptrLsn);
  S.snParent = std::move(snParent);
  S.snStart = std::move(snStart);
  S.relindCols = std::move(relindCols);
  S.relindClique = std::move(relindClique);
  S.consecutiveSums = std::move(consecutiveSums);
}

void Analyse::GenerateLayer0() {
  // linked lists of children
  std::vector<int> head, next;
  ChildrenLinkedList(snParent, head, next);

  // compute number of operations for each supernode
  std::vector<double> sn_ops(snCount);
  for (int sn = 0; sn < snCount; ++sn) {
    // supernode size
    int sz = snStart[sn + 1] - snStart[sn];

    // frontal size
    int fr = ptrLsn[sn + 1] - ptrLsn[sn];

    // number of dense operations for this supernode
    for (int i = 0; i < sz; ++i) {
      sn_ops[sn] += (double)(fr - i - 1) * (fr - i - 1);
    }
  }

  // keep track of nodes in layer0
  std::vector<int> layer0{};

  // compute number of operations to process each subtree
  std::vector<double> subtree_ops(snCount, 0.0);
  for (int sn = 0; sn < snCount; ++sn) {
    subtree_ops[sn] += sn_ops[sn];
    if (snParent[sn] != -1) {
      subtree_ops[snParent[sn]] += subtree_ops[sn];
    } else {
      // add roots in layer0
      layer0.push_back(sn);
    }
  }

  for (int iter = 0; iter < 10; ++iter) {
    printf("\nlayer0: ");
    for (int i : layer0) printf("%d ", i + 1);
    printf("\n");
    // for (int i = 0; i < layer0.size(); ++i)
    // printf("%d %f\n", layer0[i], subtree_ops[layer0[i]]);

    std::vector<double> processors(3, 0.0);

    // sort nodes in layer0 according to cost in subtree_ops
    std::sort(layer0.begin(), layer0.end(),
              [&](int a, int b) { return subtree_ops[a] > subtree_ops[b]; });

    // allocate nodes in layer0 to processors
    for (int i = 0; i < layer0.size(); ++i) {
      // find processor with lowest load
      int proc_least_load =
          std::distance(processors.begin(),
                        std::min_element(processors.begin(), processors.end()));

      processors[proc_least_load] += subtree_ops[layer0[i]];
      // printf("Put %d in proc %d\n", layer0[i], proc_least_load);
    }

    // compute imbalance ratio
    double imbalance = *std::min_element(processors.begin(), processors.end()) /
                       *std::max_element(processors.begin(), processors.end());

    printf("Imbalance: %f\n", imbalance);
    if (imbalance > 0.5) {
      double left = operations;
      for (int i = 0; i < processors.size(); ++i) {
        printf("Proc %d: %e\n", i, processors[i]);
        left -= processors[i];
      }
      printf("Left  : %e\n", left);
      break;
    }

    // find most expensive node in layer0
    auto it_node_most_exp = std::max_element(
        layer0.begin(), layer0.end(),
        [&](int a, int b) { return subtree_ops[a] < subtree_ops[b]; });
    int node_most_exp = *it_node_most_exp;

    // printf("Node to remove: %d\n", node_most_exp);

    // remove it from layer0
    layer0.erase(it_node_most_exp);

    // and add its children
    int child = head[node_most_exp];
    while (child != -1) {
      layer0.push_back(child);
      child = next[child];
    }
  }
}