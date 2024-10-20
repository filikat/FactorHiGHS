#include "Analyse.h"

#include <fstream>
#include <iostream>
#include <random>
#include <stack>

Analyse::Analyse(const std::vector<int>& rows, const std::vector<int>& ptr,
                 Symbolic& S, const std::vector<int>& order,
                 int negative_pivots)
    : S_{S} {
  // Input the symmetric matrix to be analysed in CSC format.
  // row_ind contains the row indices.
  // col_ptr contains the starting points of each column.
  // size is the number of rows/columns.
  // nonzeros is the number of nonzero entries.
  // Only the lower triangular part is used.

  n_ = ptr.size() - 1;
  nz_ = rows.size();
  negative_pivots_ = negative_pivots;

  // Create upper triangular part
  rows_upper_.resize(nz_);
  ptr_upper_.resize(n_ + 1);
  transpose(ptr, rows, ptr_upper_, rows_upper_);

  // Permute the matrix with identical permutation, to extract upper triangular
  // part, if the input is not upper triangular.
  std::vector<int> id_perm(n_);
  for (int i = 0; i < n_; ++i) id_perm[i] = i;
  permute(id_perm);

  // actual number of nonzeros of only upper triangular part
  nz_ = ptr_upper_.back();

  // number of nonzeros potentially changed after Permute.
  rows_upper_.resize(nz_);

  // double transpose to sort columns
  ptr_lower_.resize(n_ + 1);
  rows_lower_.resize(nz_);
  transpose(ptr_upper_, rows_upper_, ptr_lower_, rows_lower_);
  transpose(ptr_lower_, rows_lower_, ptr_upper_, rows_upper_);

  if (!order.empty()) {
    // inverse permutation provided by user
    iperm_ = order;
    perm_.resize(n_);
    inversePerm(iperm_, perm_);
  }

  ready_ = true;
}

int Analyse::getPermutation() {
  // Use Metis to compute a nested dissection permutation of the original matrix

  if (!perm_.empty()) {
    // permutation already provided by user
    return kRetOk;
  }

  perm_.resize(n_);
  iperm_.resize(n_);

  // Build temporary full copy of the matrix, to be used for Metis.
  // NB: Metis adjacency list should not contain the vertex itself, so diagonal
  // element is skipped.

  std::vector<int> work(n_, 0);

  // go through the columns to count nonzeros
  for (int j = 0; j < n_; ++j) {
    for (int el = ptr_upper_[j]; el < ptr_upper_[j + 1]; ++el) {
      const int i = rows_upper_[el];

      // skip diagonal entries
      if (i == j) continue;

      // nonzero in column j
      ++work[j];

      // duplicated on the lower part of column i
      ++work[i];
    }
  }

  // compute column pointers from column counts
  std::vector<int> temp_ptr(n_ + 1, 0);
  counts2Ptr(temp_ptr, work);

  std::vector<int> temp_rows(temp_ptr.back(), 0);

  for (int j = 0; j < n_; ++j) {
    for (int el = ptr_upper_[j]; el < ptr_upper_[j + 1]; ++el) {
      const int i = rows_upper_[el];

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
  // fix seed of rng inside Metis, to make it deterministic (?)
  options[METIS_OPTION_SEED] = 42;

  int status = METIS_NodeND(&n_, temp_ptr.data(), temp_rows.data(), NULL,
                            options, perm_.data(), iperm_.data());
  if (status != METIS_OK) {
    printf("Error with Metis\n");
    return kRetMetisError;
  }

  metis_order_ = iperm_;
  return kRetOk;
}

void Analyse::permute(const std::vector<int>& iperm) {
  // Symmetric permutation of the upper triangular matrix based on inverse
  // permutation iperm.
  // The resulting matrix is upper triangular, regardless of the input matrix.

  std::vector<int> work(n_, 0);

  // go through the columns to count the nonzeros
  for (int j = 0; j < n_; ++j) {
    // get new index of column
    const int col = iperm[j];

    // go through elements of column
    for (int el = ptr_upper_[j]; el < ptr_upper_[j + 1]; ++el) {
      const int i = rows_upper_[el];

      // ignore potential entries in lower triangular part
      if (i > j) continue;

      // get new index of row
      const int row = iperm[i];

      // since only upper triangular part is used, col is larger than row
      int actual_col = std::max(row, col);
      ++work[actual_col];
    }
  }

  std::vector<int> new_ptr(n_ + 1);

  // get column pointers by summing the count of nonzeros in each column.
  // copy column pointers into work
  counts2Ptr(new_ptr, work);

  std::vector<int> new_rows(new_ptr.back());

  // go through the columns to assign row indices
  for (int j = 0; j < n_; ++j) {
    // get new index of column
    const int col = iperm[j];

    // go through elements of column
    for (int el = ptr_upper_[j]; el < ptr_upper_[j + 1]; ++el) {
      const int i = rows_upper_[el];

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

  ptr_upper_ = std::move(new_ptr);
  rows_upper_ = std::move(new_rows);
}

void Analyse::eTree() {
  // Find elimination tree.
  // It works only for upper triangular matrices.
  // The tree is stored in the vector parent:
  //  parent[i] = j
  // means that j is the parent of i in the tree.
  // For the root(s) of the tree, parent[root] = -1.

  parent_.resize(n_);
  std::vector<int> ancestor(n_);
  int next{};

  for (int j = 0; j < n_; ++j) {
    // initialize parent and ancestor, which are still unknown
    parent_[j] = -1;
    ancestor[j] = -1;

    for (int el = ptr_upper_[j]; el < ptr_upper_[j + 1]; ++el) {
      for (int i = rows_upper_[el]; i != -1 && i < j; i = next) {
        // next is used to move up the tree
        next = ancestor[i];

        // ancestor keeps track of the known part of the tree, to avoid
        // repeating (aka path compression): from j there is a known path to i
        ancestor[i] = j;

        if (next == -1) parent_[i] = j;
      }
    }
  }
}

void Analyse::postorder() {
  // Find a postordering of the elimination tree using depth first search

  postorder_.resize(n_);

  // create linked list of children
  std::vector<int> head, next;
  childrenLinkedList(parent_, head, next);

  // Execute depth first search only for root node(s)
  int start{};
  for (int node = 0; node < n_; ++node) {
    if (parent_[node] == -1) {
      dfsPostorder(node, start, head, next, postorder_);
    }
  }

  // Permute elimination tree based on postorder
  std::vector<int> ipost(n_);
  inversePerm(postorder_, ipost);
  std::vector<int> new_parent(n_);
  for (int i = 0; i < n_; ++i) {
    if (parent_[i] != -1) {
      new_parent[ipost[i]] = ipost[parent_[i]];
    } else {
      new_parent[ipost[i]] = -1;
    }
  }
  parent_ = std::move(new_parent);

  // Permute matrix based on postorder
  permute(ipost);

  // double transpose to sort columns and compute lower part
  transpose(ptr_upper_, rows_upper_, ptr_lower_, rows_lower_);
  transpose(ptr_lower_, rows_lower_, ptr_upper_, rows_upper_);

  // Update perm and iperm
  permuteVector(perm_, postorder_);
  inversePerm(perm_, iperm_);
}

void Analyse::colCount() {
  // Columns count using skeleton matrix.
  // Taken from Tim Davis "Direct Methods for Sparse Linear Systems".

  std::vector<int> first(n_, -1);
  std::vector<int> ancestor(n_, -1);
  std::vector<int> max_first(n_, -1);
  std::vector<int> prev_leaf(n_, -1);

  col_count_.resize(n_);

  // find first descendant
  for (int k = 0; k < n_; ++k) {
    int j = k;
    col_count_[j] = (first[j] == -1) ? 1 : 0;
    while (j != -1 && first[j] == -1) {
      first[j] = k;
      j = parent_[j];
    }
  }

  // each node belongs to a separate set
  for (int j = 0; j < n_; j++) ancestor[j] = j;

  for (int k = 0; k < n_; ++k) {
    const int j = k;

    // if not a root, decrement
    if (parent_[j] != -1) col_count_[parent_[j]]--;

    // process edges of matrix
    for (int el = ptr_lower_[j]; el < ptr_lower_[j + 1]; ++el) {
      processEdge(j, rows_lower_[el], first, max_first, col_count_, prev_leaf,
                  ancestor);
    }

    if (parent_[j] != -1) ancestor[j] = parent_[j];
  }

  // sum contributions from each child
  for (int j = 0; j < n_; ++j) {
    if (parent_[j] != -1) {
      col_count_[parent_[j]] += col_count_[j];
    }
  }

  // compute nonzeros of L
  operations_no_relax_ = 0.0;
  nz_factor_ = 0;
  for (int j = 0; j < n_; ++j) {
    nz_factor_ += (double)col_count_[j];
    operations_no_relax_ += (double)(col_count_[j] - 1) * (col_count_[j] - 1);
  }
}

void Analyse::fundamentalSupernodes() {
  // Find fundamental supernodes.

  // isSN[i] is true if node i is the start of a fundamental supernode
  std::vector<bool> is_sn(n_, false);

  std::vector<int> prev_nonz(n_, -1);

  // compute sizes of subtrees
  std::vector<int> subtree_sizes(n_);
  subtreeSize(parent_, subtree_sizes);

  for (int j = 0; j < n_; ++j) {
    for (int el = ptr_lower_[j]; el < ptr_lower_[j + 1]; ++el) {
      const int i = rows_lower_[el];
      const int k = prev_nonz[i];

      // mark as fundamental sn, nodes which are leaf of subtrees
      if (k < j - subtree_sizes[j] + 1) {
        is_sn[j] = true;
      }

      // mark as fundamental sn, nodes which have more than one child
      if (parent_[i] != -1 &&
          subtree_sizes[i] + 1 != subtree_sizes[parent_[i]]) {
        is_sn[parent_[i]] = true;
      }

      prev_nonz[i] = j;
    }
  }

  // create information about fundamental supernodes
  sn_belong_.resize(n_);
  int sn_number = -1;
  for (int i = 0; i < n_; ++i) {
    // if isSN[i] is true, then node i is the start of a new supernode
    if (is_sn[i]) ++sn_number;

    // mark node i as belonging to the current supernode
    sn_belong_[i] = sn_number;
  }

  // number of supernodes found
  sn_count_ = sn_belong_.back() + 1;

  // fsn_ptr contains pointers to the starting node of each supernode
  sn_start_.resize(sn_count_ + 1);
  int next = 0;
  for (int i = 0; i < n_; ++i) {
    if (is_sn[i]) {
      sn_start_[next] = i;
      ++next;
    }
  }
  sn_start_[next] = n_;

  // build supernodal elimination tree
  sn_parent_.resize(sn_count_);
  for (int i = 0; i < sn_count_ - 1; ++i) {
    int j = parent_[sn_start_[i + 1] - 1];
    if (j != -1) {
      sn_parent_[i] = sn_belong_[j];
    } else {
      sn_parent_[i] = -1;
    }
  }
  sn_parent_.back() = -1;
}

void Analyse::relaxSupernodes() {
  // Child which produces smallest number of fake nonzeros is merged if
  // resulting sn has fewer than max_artificial_nz fake nonzeros.
  // Multiple values of max_artificial_nz are tried, chosen with bisection
  // method, until the percentage of artificial nonzeros is in the range [1,2]%.

  int max_artificial_nz = kStartThreshRelax;
  int largest_below = -1;
  int smallest_above = -1;

  for (int iter = 0; iter < kMaxIterRelax; ++iter) {
    // =================================================
    // Build information about supernodes
    // =================================================
    std::vector<int> sn_size(sn_count_);
    std::vector<int> clique_size(sn_count_);
    fake_nz_.assign(sn_count_, 0);
    for (int i = 0; i < sn_count_; ++i) {
      sn_size[i] = sn_start_[i + 1] - sn_start_[i];
      clique_size[i] = col_count_[sn_start_[i]] - sn_size[i];
      fake_nz_[i] = 0;
    }

    // build linked lists of children
    std::vector<int> first_child, next_child;
    childrenLinkedList(sn_parent_, first_child, next_child);

    // =================================================
    // Merge supernodes
    // =================================================
    merged_into_.assign(sn_count_, -1);
    merged_sn_ = 0;

    for (int sn = 0; sn < sn_count_; ++sn) {
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
          const int rows_filled =
              sn_size[sn] + clique_size[sn] - clique_size[child];

          // how many zero entries would become nonzero
          const int nz_added = rows_filled * sn_size[child];

          // how many artificial nonzeros would the merged supernode have
          const int total_art_nz = nz_added + fake_nz_[sn] + fake_nz_[child];

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
          fake_nz_[sn] = nz_fakenz;

          // count number of merged supernodes
          ++merged_sn_;

          // save information about merging of supernodes
          merged_into_[child_fakenz] = sn;

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
    for (int sn = 0; sn < sn_count_; ++sn) {
      if (merged_into_[sn] == -1) {
        temp_art_nz += fake_nz_[sn];

        const double nn = sn_size[sn];
        const double cc = clique_size[sn];
        temp_art_ops += (nn + cc) * (nn + cc) * nn - (nn + cc) * nn * (nn + 1) +
                        nn * (nn + 1) * (2 * nn + 1) / 6;
      }
    }
    temp_art_ops -= operations_no_relax_;

    // if enough fake nz or ops have been added, stop.
    // double ratio_fake = temp_art_nz / (nzL + temp_art_nz);
    const double ratio_fake =
        temp_art_ops / (temp_art_ops + operations_no_relax_);

    // try to find ratio in interval [0.01,0.02] using bisection
    if (ratio_fake < kLowerRatioRelax) {
      // ratio too small
      largest_below = max_artificial_nz;
      if (smallest_above == -1) {
        max_artificial_nz *= 2;
      } else {
        max_artificial_nz = (largest_below + smallest_above) / 2;
      }
    } else if (ratio_fake > kUpperRatioRelax) {
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

void Analyse::relaxSupernodes2() {
  // Smallest child is merged with parent, if child is small enough.

  // =================================================
  // build information about supernodes
  // =================================================
  std::vector<int> sn_size(sn_count_);
  std::vector<int> clique_size(sn_count_);
  fake_nz_.assign(sn_count_, 0);
  for (int i = 0; i < sn_count_; ++i) {
    sn_size[i] = sn_start_[i + 1] - sn_start_[i];
    clique_size[i] = col_count_[sn_start_[i]] - sn_size[i];
    fake_nz_[i] = 0;
  }

  // build linked lists of children
  std::vector<int> first_child, next_child;
  childrenLinkedList(sn_parent_, first_child, next_child);

  // =================================================
  // Merge supernodes
  // =================================================
  merged_into_.assign(sn_count_, -1);
  merged_sn_ = 0;

  for (int sn = 0; sn < sn_count_; ++sn) {
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
        const int rows_filled =
            sn_size[sn] + clique_size[sn] - clique_size[child];

        // how many zero entries would become nonzero
        const int nz_added = rows_filled * sn_size[child];

        // how many artificial nonzeros would the merged supernode have
        const int total_art_nz = nz_added + fake_nz_[sn] + fake_nz_[child];

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
        fake_nz_[sn] = nz_smallest;

        // count number of merged supernodes
        ++merged_sn_;

        // save information about merging of supernodes
        merged_into_[child_smallest] = sn;

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

void Analyse::afterRelaxSn() {
  // number of new supernodes
  const int new_snCount = sn_count_ - merged_sn_;

  // keep track of number of row indices needed for each supernode
  sn_indices_.assign(new_snCount, 0);

  // =================================================
  // Create supernodal permutation
  // =================================================

  // permutation of supernodes needed after merging
  std::vector<int> sn_perm(sn_count_);

  // number of new sn that includes the old sn
  std::vector<int> new_id(sn_count_);

  // new sn pointer vector
  std::vector<int> new_snStart(new_snCount + 1);

  // keep track of the children merged into a given supernode
  std::vector<std::vector<int>> received_from(sn_count_, std::vector<int>());

  // index to write into sn_perm
  int start_perm{};

  // index to write into new_snStart
  int snStart_ind{};

  // next available number for new sn numbering
  int next_id{};

  for (int sn = 0; sn < sn_count_; ++sn) {
    if (merged_into_[sn] > -1) {
      // Current sn was merged into its parent.
      // Save information about which supernode sn was merged into
      received_from[merged_into_[sn]].push_back(sn);
    } else {
      // Current sn was not merged into its parent.
      // It is one of the new sn.

      // Add merged supernodes to the permutation, recursively.

      ++snStart_ind;

      std::stack<int> toadd;
      toadd.push(sn);

      while (!toadd.empty()) {
        const int current = toadd.top();

        if (!received_from[current].empty()) {
          for (int i : received_from[current]) toadd.push(i);
          received_from[current].clear();
        } else {
          toadd.pop();
          sn_perm[start_perm++] = current;
          new_id[current] = next_id;

          // count number of nodes in each new supernode
          new_snStart[snStart_ind] +=
              sn_start_[current + 1] - sn_start_[current];
        }
      }

      // keep track of total number of artificial nonzeros
      artificial_nz_ += fake_nz_[sn];

      // Compute number of indices for new sn.
      // This is equal to the number of columns in the new sn plus the clique
      // size of the original supernode where the children where merged.
      sn_indices_[next_id] = new_snStart[snStart_ind] +
                             col_count_[sn_start_[sn]] - sn_start_[sn + 1] +
                             sn_start_[sn];

      ++next_id;
    }
  }

  // new_snStart contain the number of cols in each new sn.
  // sum them to obtain the sn pointers.
  for (int i = 0; i < new_snCount; ++i) {
    new_snStart[i + 1] += new_snStart[i];
  }

  // include artificial nonzeros in the nonzeros of the factor
  nz_factor_ += (double)artificial_nz_;

  // compute number of flops needed for the factorization
  operations_ = 0.0;
  for (int sn = 0; sn < new_snCount; ++sn) {
    const double colcount_sn = (double)sn_indices_[sn];
    for (int i = 0; i < new_snStart[sn + 1] - new_snStart[sn]; ++i) {
      operations_ += (colcount_sn - i - 1) * (colcount_sn - i - 1);
    }
  }

  // =================================================
  // Create nodal permutation
  // =================================================
  // Given the supernodal permutation, find the nodal permutation needed after
  // sn merging.

  // permutation to apply to the existing one
  std::vector<int> new_perm(n_);

  // index to write into new_perm
  int start{};

  for (int i = 0; i < sn_count_; ++i) {
    const int sn = sn_perm[i];
    for (int j = sn_start_[sn]; j < sn_start_[sn + 1]; ++j) {
      new_perm[start++] = j;
    }
  }

  // obtain inverse permutation
  std::vector<int> new_iperm(n_);
  inversePerm(new_perm, new_iperm);

  // =================================================
  // Create new sn elimination tree
  // =================================================
  std::vector<int> new_snParent(new_snCount, -1);
  for (int i = 0; i < sn_count_; ++i) {
    if (sn_parent_[i] == -1) continue;

    const int ii = new_id[i];
    const int pp = new_id[sn_parent_[i]];

    if (ii == pp) continue;

    new_snParent[ii] = pp;
  }

  // =================================================
  // Save new information
  // =================================================

  // build new snBelong, i.e., the sn to which each column belongs
  for (int sn = 0; sn < sn_count_; ++sn) {
    for (int i = sn_start_[sn]; i < sn_start_[sn + 1]; ++i) {
      sn_belong_[i] = new_id[sn];
    }
  }
  permuteVector(sn_belong_, new_perm);

  permuteVector(col_count_, new_perm);

  // Overwrite previous data
  sn_parent_ = std::move(new_snParent);
  sn_start_ = std::move(new_snStart);
  sn_count_ = new_snCount;

  // Permute matrix based on new permutation
  permute(new_iperm);

  // double transpose to sort columns and compute lower part
  transpose(ptr_upper_, rows_upper_, ptr_lower_, rows_lower_);
  transpose(ptr_lower_, rows_lower_, ptr_upper_, rows_upper_);

  // Update perm and iperm
  permuteVector(perm_, new_perm);
  inversePerm(perm_, iperm_);
}

void Analyse::snPattern() {
  // number of total indices needed
  int indices{};

  for (int i : sn_indices_) indices += i;

  // allocate space for sn pattern
  rows_sn_.resize(indices);
  ptr_sn_.resize(sn_count_ + 1);

  // keep track of visited supernodes
  std::vector<int> mark(sn_count_, -1);

  // compute column pointers of L
  std::vector<int> work(sn_indices_);
  counts2Ptr(ptr_sn_, work);

  // consider each row
  for (int i = 0; i < n_; ++i) {
    // for all entries in the row of lower triangle
    for (int el = ptr_upper_[i]; el < ptr_upper_[i + 1]; ++el) {
      // there is nonzero (i,j)
      const int j = rows_upper_[el];

      // supernode to which column j belongs to
      int snj = sn_belong_[j];

      // while supernodes are not yet considered
      while (snj != -1 && mark[snj] != i) {
        // we may end up too far
        if (sn_start_[snj] > i) break;

        // supernode snj is now considered for row i
        mark[snj] = i;

        // there is a nonzero entry in supernode snj at row i
        rows_sn_[work[snj]++] = i;

        // go up the elimination tree
        snj = sn_parent_[snj];
      }
    }
  }
}

void Analyse::relativeIndCols() {
  // Find the relative indices of the original column wrt the frontal matrix of
  // the corresponding supernode

  relind_cols_.resize(nz_);

  // go through the supernodes
  for (int sn = 0; sn < sn_count_; ++sn) {
    const int ptL_start = ptr_sn_[sn];
    const int ptL_end = ptr_sn_[sn + 1];

    // go through the columns of the supernode
    for (int col = sn_start_[sn]; col < sn_start_[sn + 1]; ++col) {
      // go through original column and supernodal column
      int ptA = ptr_lower_[col];
      int ptL = ptL_start;

      // offset wrt ptrLower[col]
      int index{};

      // size of the column of the original matrix
      int col_size = ptr_lower_[col + 1] - ptr_lower_[col];

      while (ptL < ptL_end) {
        // if found all the relative indices that are needed, stop
        if (index == col_size) {
          break;
        }

        // check if indices coincide
        if (rows_sn_[ptL] == rows_lower_[ptA]) {
          // yes: save relative index and move pointers forward
          relind_cols_[ptr_lower_[col] + index] = ptL - ptL_start;
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

void Analyse::relativeIndClique() {
  // Find the relative indices of the child clique wrt the frontal matrix of the
  // parent supernode

  relind_clique_.resize(sn_count_);
  consecutive_sums_.resize(sn_count_);

  for (int sn = 0; sn < sn_count_; ++sn) {
    // if there is no parent, skip supernode
    if (sn_parent_[sn] == -1) continue;

    // number of nodes in the supernode
    const int sn_size = sn_start_[sn + 1] - sn_start_[sn];

    // column of the first node in the supernode
    const int j = sn_start_[sn];

    // size of the first column of the supernode
    const int sn_column_size = ptr_sn_[sn + 1] - ptr_sn_[sn];

    // size of the clique of the supernode
    const int sn_clique_size = sn_column_size - sn_size;

    // count number of assembly operations during factorize
    operations_assembly_ += sn_clique_size * (sn_clique_size + 1) / 2;

    relind_clique_[sn].resize(sn_clique_size);

    // iterate through the clique of sn
    int ptr_current = ptr_sn_[sn] + sn_size;

    // iterate through the full column of parent sn
    int ptr_parent = ptr_sn_[sn_parent_[sn]];

    // keep track of start and end of parent sn column
    const int ptr_parent_start = ptr_parent;
    const int ptr_parent_end = ptr_sn_[sn_parent_[sn] + 1];

    // where to write into relind
    int index{};

    // iterate though the column of the parent sn
    while (ptr_parent < ptr_parent_end) {
      // if found all the relative indices that are needed, stop
      if (index == sn_clique_size) {
        break;
      }

      // check if indices coincide
      if (rows_sn_[ptr_current] == rows_sn_[ptr_parent]) {
        // yes: save relative index and move pointers forward
        relind_clique_[sn][index] = ptr_parent - ptr_parent_start;
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
    consecutive_sums_[sn].resize(sn_clique_size);
    for (int i = 0; i < sn_clique_size - 1; ++i) {
      consecutive_sums_[sn][i] =
          relind_clique_[sn][i + 1] - relind_clique_[sn][i];
    }

    // Number of consecutive sums that can be done in one blas call.
    consecutive_sums_[sn].back() = 1;
    for (int i = sn_clique_size - 2; i >= 0; --i) {
      if (consecutive_sums_[sn][i] > 1) {
        consecutive_sums_[sn][i] = 1;
      } else if (consecutive_sums_[sn][i] == 1) {
        consecutive_sums_[sn][i] = consecutive_sums_[sn][i + 1] + 1;
      } else {
        printf("Error in consecutiveSums %d\n", consecutive_sums_[sn][i]);
      }
    }
  }
}

bool Analyse::check() const {
  // Check that the symbolic factorization is correct, by using dense linear
  // algebra operations.
  // Return true if check is successful, or if matrix is too large.
  // To be used for debug.

  // Check symbolic factorization
  if (n_ > 5000) {
    printf("\n==> Matrix is too large for dense check\n\n");
    return true;
  }

  // initialize random number generator (to avoid numerical cancellation)
  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_real_distribution<double> distr(0.1, 10.0);

  // assemble sparse matrix into dense matrix
  std::vector<double> M(n_ * n_);
  for (int col = 0; col < n_; ++col) {
    for (int el = ptr_upper_[col]; el < ptr_upper_[col + 1]; ++el) {
      int row = rows_upper_[el];

      // insert random element in position (row,col)
      M[row + col * n_] = distr(rng);

      // guarantee matrix is diagonally dominant (thus positive definite)
      if (row == col) {
        M[row + col * n_] += n_ * 10;
      }
    }
  }

  // use Lapack to factorize the dense matrix
  char uplo = 'U';
  int N = n_;
  int info;
  dpotrf_(&uplo, &N, M.data(), &N, &info);
  if (info != 0) {
    printf("\n==> dpotrf failed\n\n");
    return false;
  }

  // assemble expected sparsity pattern into dense matrix
  std::vector<bool> L(n_ * n_);
  for (int sn = 0; sn < sn_count_; ++sn) {
    for (int col = sn_start_[sn]; col < sn_start_[sn + 1]; ++col) {
      for (int el = ptr_sn_[sn]; el < ptr_sn_[sn + 1]; ++el) {
        int row = rows_sn_[el];
        if (row < col) continue;
        L[row + n_ * col] = true;
      }
    }
  }

  int zeros_found{};
  int wrong_entries{};

  // check how many entries do not correspond
  for (int j = 0; j < n_; ++j) {
    for (int i = 0; i < n_; ++i) {
      double val_M = M[j + n_ * i];
      bool val_L = L[i + n_ * j];

      if (val_L && val_M == 0.0) {
        // count number of fake zeros found, to confront it with artificialNz
        ++zeros_found;
      } else if (!val_L && val_M != 0.0) {
        printf("==> (%d,%d) Found nonzero, expected zero\n", i, j);
        ++wrong_entries;
      }
    }
  }

  if (wrong_entries == 0 && zeros_found == artificial_nz_) {
    printf("\n==> Analyse check successful\n\n");
    return true;
  } else {
    printf("\n==> Analyse check failed\n\n");
    return false;
  }
}

void Analyse::generateLayer0(double imbalance_ratio) {
  if (S_.n_threads_ <= 1) return;

  // linked lists of children
  std::vector<int> head, next;
  childrenLinkedList(sn_parent_, head, next);

  double total_ops = operations_;

  // compute number of operations for each supernode
  std::vector<double> sn_ops(sn_count_);
  for (int sn = 0; sn < sn_count_; ++sn) {
    // supernode size
    const int sz = sn_start_[sn + 1] - sn_start_[sn];

    // frontal size
    const int fr = ptr_sn_[sn + 1] - ptr_sn_[sn];

    // number of operations for this supernode
    for (int i = 0; i < sz; ++i) {
      sn_ops[sn] += (double)(fr - i - 1) * (fr - i - 1);
    }

    // add assembly operations times 100 to the parent
    if (sn_parent_[sn] != -1) {
      const int ldc = fr - sz;
      sn_ops[sn_parent_[sn]] += ldc * (ldc + 1) / 2 * 100;
      total_ops += ldc * (ldc + 1) / 2 * 100;
    }
  }

  // keep track of nodes in layer0
  std::vector<int> layer0{};

  // compute number of operations to process each subtree
  std::vector<double> subtree_ops(sn_count_, 0.0);
  for (int sn = 0; sn < sn_count_; ++sn) {
    subtree_ops[sn] += sn_ops[sn];
    if (sn_parent_[sn] != -1) {
      subtree_ops[sn_parent_[sn]] += subtree_ops[sn];
    } else {
      // add roots in layer0
      layer0.push_back(sn);
    }
  }

  for (int iter = 0; iter < 10; ++iter) {
    ops_per_thread_.assign(S_.n_threads_, 0.0);
    subtrees_per_thread_.assign(S_.n_threads_, {});

    // sort nodes in layer0 according to cost in subtree_ops.
    // expensive nodes are ordered last
    std::sort(layer0.begin(), layer0.end(),
              [&](int a, int b) { return subtree_ops[a] < subtree_ops[b]; });

    // allocate nodes in layer0 to processors:
    // least loaded processor receives the next most expensive node
    for (int i = layer0.size() - 1; i >= 0; --i) {
      // find processor with lowest load
      const int proc_least_load = std::distance(
          ops_per_thread_.begin(),
          std::min_element(ops_per_thread_.begin(), ops_per_thread_.end()));

      ops_per_thread_[proc_least_load] += subtree_ops[layer0[i]];
      subtrees_per_thread_[proc_least_load].push_back(layer0[i]);
    }

    // compute imbalance ratio
    const double imbalance =
        *std::min_element(ops_per_thread_.begin(), ops_per_thread_.end()) /
        *std::max_element(ops_per_thread_.begin(), ops_per_thread_.end());

    if (imbalance > imbalance_ratio) break;

    // most expensive node in layer0 is the last one, because layer0 is sorted
    const int node_most_exp = layer0.back();

    // if node to be removed does not have children, stop
    if (head[node_most_exp] == -1) break;

    // remove it from layer0
    layer0.pop_back();
    ++sn_above_layer0_;

    // and add its children
    int child = head[node_most_exp];
    while (child != -1) {
      layer0.push_back(child);
      child = next[child];
    }
  }
}

void Analyse::computeStorage(int fr, int sz, int& fr_entries,
                             int& cl_entries) const {
  // compute storage required by frontal and clique, based on the format used

  const int cl = fr - sz;
  switch (S_.formatType()) {
    case FormatType::Full:
      // full format stores a rectangle for each
      fr_entries = fr * sz;
      cl_entries = cl * cl;
      break;

    case FormatType::HybridPacked:
    case FormatType::HybridHybrid:
      // frontal is stored as a trapezoid
      fr_entries = fr * sz - sz * (sz - 1) / 2;

      // clique is stored as a collection of rectangles
      const int nb = S_.blockSize();
      const int n_blocks = (cl - 1) / nb + 1;
      int schur_size{};
      for (int j = 0; j < n_blocks; ++j) {
        const int jb = std::min(nb, cl - j * nb);
        schur_size += (cl - j * nb) * jb;
      }
      cl_entries = schur_size;
      break;
  }
}

void Analyse::computeStorage() {
  std::vector<int> clique_entries(sn_count_);
  std::vector<int> frontal_entries(sn_count_);
  std::vector<int> storage(sn_count_);
  std::vector<int> storage_factors(sn_count_);
  std::vector<int> stack_size(sn_count_);

  // initialize data of supernodes
  for (int sn = 0; sn < sn_count_; ++sn) {
    // supernode size
    const int sz = sn_start_[sn + 1] - sn_start_[sn];

    // frontal size
    const int fr = ptr_sn_[sn + 1] - ptr_sn_[sn];

    // compute storage based on format used
    computeStorage(fr, sz, frontal_entries[sn], clique_entries[sn]);

    max_clique_entries_ = std::max(max_clique_entries_, clique_entries[sn]);

    // compute number of entries in factors within the subtree
    storage_factors[sn] += frontal_entries[sn];
    if (sn_parent_[sn] != -1)
      storage_factors[sn_parent_[sn]] += storage_factors[sn];
  }

  // linked lists of children
  std::vector<int> head, next;
  childrenLinkedList(sn_parent_, head, next);

  // go through the supernodes
  for (int sn = 0; sn < sn_count_; ++sn) {
    // leaf node
    if (head[sn] == -1) {
      storage[sn] = frontal_entries[sn] + clique_entries[sn];
      // stack_size[sn] = clique_entries[sn];
      stack_size[sn] = 0;
      continue;
    }

    int clique_total_entries{};
    int factors_total_entries{};
    int child = head[sn];
    while (child != -1) {
      int value =
          storage[child] - clique_entries[child] - storage_factors[child];
      clique_total_entries += clique_entries[child];
      factors_total_entries += storage_factors[child];
      child = next[child];
    }

    // Compute storage, based on order just found.
    // storage is found as max(storage_1,storage_2), where
    // storage_1 = max_j storage[j] + \sum_{k up to j-1} clique_entries[k] +
    //                                                   storage_factors[k]
    // storage_2 = frontal_entries + clique_entries + clique_total_entries +
    //             factors_total_entries
    const int storage_2 = frontal_entries[sn] + clique_entries[sn] +
                          clique_total_entries + factors_total_entries;
    const int storage_2_stack = clique_total_entries;

    int clique_partial_entries{};
    int factors_partial_entries{};
    int storage_1{};
    int storage_1_stack{};

    child = head[sn];
    while (child != -1) {
      int current =
          storage[child] + clique_partial_entries + factors_partial_entries;

      storage_1_stack =
          std::max(storage_1_stack, stack_size[child] + clique_partial_entries);

      clique_partial_entries += clique_entries[child];
      factors_partial_entries += storage_factors[child];
      storage_1 = std::max(storage_1, current);

      child = next[child];
    }
    storage[sn] = std::max(storage_1, storage_2);
    stack_size[sn] = std::max(storage_1_stack, storage_2_stack);

    // save max storage needed, multiply by 8 because double needs 8 bytes
    max_storage_ = std::max(max_storage_, 8 * storage[sn]);
    max_stack_entries_ = std::max(max_stack_entries_, stack_size[sn]);
  }
}

void Analyse::reorderChildren() {
  std::vector<int> clique_entries(sn_count_);
  std::vector<int> frontal_entries(sn_count_);
  std::vector<int> storage(sn_count_);
  std::vector<int> storage_factors(sn_count_);

  // initialize data of supernodes
  for (int sn = 0; sn < sn_count_; ++sn) {
    // supernode size
    const int sz = sn_start_[sn + 1] - sn_start_[sn];

    // frontal size
    const int fr = col_count_[sn_start_[sn]];

    // compute storage based on format used
    computeStorage(fr, sz, frontal_entries[sn], clique_entries[sn]);

    // compute number of entries in factors within the subtree
    storage_factors[sn] += frontal_entries[sn];
    if (sn_parent_[sn] != -1)
      storage_factors[sn_parent_[sn]] += storage_factors[sn];
  }

  // linked lists of children
  std::vector<int> head, next;
  childrenLinkedList(sn_parent_, head, next);

  // go through the supernodes
  for (int sn = 0; sn < sn_count_; ++sn) {
    // leaf node
    if (head[sn] == -1) {
      storage[sn] = frontal_entries[sn] + clique_entries[sn];
      continue;
    }

    // save children and values to sort
    std::vector<std::pair<int, int>> children{};
    int child = head[sn];
    while (child != -1) {
      int value =
          storage[child] - clique_entries[child] - storage_factors[child];
      children.push_back({child, value});
      child = next[child];
    }

    // sort children in decreasing order of the values
    std::sort(children.begin(), children.end(),
              [&](std::pair<int, int>& a, std::pair<int, int>& b) {
                return a.second > b.second;
              });

    // modify linked lists with new order of children
    head[sn] = children.front().first;
    for (int i = 0; i < children.size() - 1; ++i) {
      next[children[i].first] = children[i + 1].first;
    }
    next[children.back().first] = -1;
  }

  // =================================================
  // Create supernodal permutation
  // =================================================
  // build supernodal permutation with dfs
  std::vector<int> sn_perm(sn_count_);
  int start{};
  for (int sn = 0; sn < sn_count_; ++sn) {
    if (sn_parent_[sn] == -1) dfsPostorder(sn, start, head, next, sn_perm);
  }

  // =================================================
  // Create nodal permutation
  // =================================================
  // Given the supernodal permutation, find the nodal permutation

  // permutation to apply to the existing one
  std::vector<int> new_perm(n_);

  // index to write into new_perm
  start = 0;

  for (int i = 0; i < sn_count_; ++i) {
    const int sn = sn_perm[i];
    for (int j = sn_start_[sn]; j < sn_start_[sn + 1]; ++j) {
      new_perm[start++] = j;
    }
  }

  // obtain inverse permutation
  std::vector<int> new_iperm(n_);
  inversePerm(new_perm, new_iperm);

  // =================================================
  // Create new sn elimination tree
  // =================================================
  std::vector<int> isn_perm(sn_count_);
  inversePerm(sn_perm, isn_perm);
  std::vector<int> new_sn_parent(sn_count_);
  for (int i = 0; i < sn_count_; ++i) {
    if (sn_parent_[i] != -1) {
      new_sn_parent[isn_perm[i]] = isn_perm[sn_parent_[i]];
    } else {
      new_sn_parent[isn_perm[i]] = -1;
    }
  }

  // =================================================
  // Create new snBelong
  // =================================================

  // build new snBelong, i.e., the sn to which each colum belongs
  for (int sn = 0; sn < sn_count_; ++sn) {
    for (int i = sn_start_[sn]; i < sn_start_[sn + 1]; ++i) {
      sn_belong_[i] = isn_perm[sn];
    }
  }
  permuteVector(sn_belong_, new_perm);

  // permute other vectors that may be needed
  permuteVector(col_count_, new_perm);
  permuteVector(sn_indices_, sn_perm);

  // =================================================
  // Create new snStart
  // =================================================
  std::vector<int> cols_per_sn(sn_count_);

  // compute size of each supernode
  for (int sn = 0; sn < sn_count_; ++sn) {
    cols_per_sn[sn] = sn_start_[sn + 1] - sn_start_[sn];
  }

  // permute according to new order of supernodes
  permuteVector(cols_per_sn, sn_perm);

  // sum number of columns to obtain pointers
  for (int i = 0; i < sn_count_ - 1; ++i) {
    cols_per_sn[i + 1] += cols_per_sn[i];
  }

  for (int i = 0; i < sn_count_; ++i) {
    sn_start_[i + 1] = cols_per_sn[i];
  }

  // =================================================
  // Save new data
  // =================================================

  // Overwrite previous data
  sn_parent_ = std::move(new_sn_parent);

  // Permute matrix based on new permutation
  permute(new_iperm);

  // double transpose to sort columns and compute lower part
  transpose(ptr_upper_, rows_upper_, ptr_lower_, rows_lower_);
  transpose(ptr_lower_, rows_lower_, ptr_upper_, rows_upper_);

  // Update perm and iperm
  permuteVector(perm_, new_perm);
  inversePerm(perm_, iperm_);
}

int Analyse::run() {
  // Perform analyse phase and store the result into the symbolic object S.
  // After Run returns, the Analyse object is not valid.

  if (!ready_) return kRetGeneric;

  S_.times().resize(kTimeSize);

  Clock clock_total{};
  Clock clock_items{};

#ifdef COARSE_TIMING
  clock_total.start();
#endif

#ifdef FINE_TIMING
  clock_items.start();
#endif
  int metis_status = getPermutation();
  if (metis_status) return kRetMetisError;
#ifdef FINE_TIMING
  S_.times(kTimeAnalyseMetis) += clock_items.stop();
#endif

#ifdef FINE_TIMING
  clock_items.start();
#endif
  permute(iperm_);
  eTree();
  postorder();
#ifdef FINE_TIMING
  S_.times(kTimeAnalyseTree) += clock_items.stop();
#endif

#ifdef FINE_TIMING
  clock_items.start();
#endif
  colCount();
#ifdef FINE_TIMING
  S_.times(kTimeAnalyseCount) += clock_items.stop();
#endif

#ifdef FINE_TIMING
  clock_items.start();
#endif
  fundamentalSupernodes();
  relaxSupernodes();
  afterRelaxSn();
#ifdef FINE_TIMING
  S_.times(kTimeAnalyseSn) += clock_items.stop();
#endif

#ifdef FINE_TIMING
  clock_items.start();
#endif
  reorderChildren();
#ifdef FINE_TIMING
  S_.times(kTimeAnalyseReorder) += clock_items.stop();
#endif

#ifdef FINE_TIMING
  clock_items.start();
#endif
  snPattern();
#ifdef FINE_TIMING
  S_.times(kTimeAnalysePattern) += clock_items.stop();
#endif

#ifdef FINE_TIMING
  clock_items.start();
#endif
  relativeIndCols();
  relativeIndClique();
  computeStorage();
#ifdef FINE_TIMING
  S_.times(kTimeAnalyseRelInd) += clock_items.stop();
#endif

#ifdef FINE_TIMING
  clock_items.start();
#endif
  generateLayer0(0.7);
#ifdef FINE_TIMING
  S_.times(kTimeAnalyseLayer0) += clock_items.stop();
#endif

  // move relevant stuff into S
  S_.n_ = n_;
  S_.nz_ = nz_factor_;
  S_.fillin_ = (double)nz_factor_ / nz_;
  S_.sn_ = sn_count_;
  S_.artificial_nz_ = artificial_nz_;
  S_.artificial_ops_ = (double)operations_ - operations_no_relax_;
  S_.assembly_ops_ = operations_assembly_;
  S_.largest_front_ = *std::max_element(sn_indices_.begin(), sn_indices_.end());
  S_.max_storage_ = max_storage_;
  S_.max_stack_entries_ = max_stack_entries_;
  S_.max_clique_entries_ = max_clique_entries_;
  S_.sn_above_layer0_ = sn_above_layer0_;

  // compute largest supernode
  std::vector<int> temp(sn_start_);
  for (int i = sn_count_; i > 0; --i) temp[i] -= temp[i - 1];
  S_.largest_sn_ = *std::max_element(temp.begin(), temp.end());

  // initialize sign of pivots and permute them
  S_.pivot_sign_.insert(S_.pivot_sign_.end(), negative_pivots_, -1);
  S_.pivot_sign_.insert(S_.pivot_sign_.end(), n_ - negative_pivots_, 1);
  permuteVector(S_.pivot_sign_, perm_);

  S_.dense_ops_ = operations_;
  S_.iperm_ = std::move(iperm_);
  S_.rows_ = std::move(rows_sn_);
  S_.ptr_ = std::move(ptr_sn_);
  S_.sn_parent_ = std::move(sn_parent_);
  S_.sn_start_ = std::move(sn_start_);
  S_.relind_cols_ = std::move(relind_cols_);
  S_.relind_clique_ = std::move(relind_clique_);
  S_.consecutive_sums_ = std::move(consecutive_sums_);
  S_.subtrees_per_thread_ = std::move(subtrees_per_thread_);
  S_.ops_per_thread_ = std::move(ops_per_thread_);

#ifdef COARSE_TIMING
  S_.times(kTimeAnalyse) += clock_total.stop();
#endif

  return kRetOk;
}