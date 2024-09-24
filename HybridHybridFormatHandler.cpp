#include "HybridHybridFormatHandler.h"

void HybridHybridFormatHandler::initFrontal() {
  // frontal is initialized to zero
  frontal_->resize(ldf_ * sn_size_ - sn_size_ * (sn_size_ - 1) / 2);
}

void HybridHybridFormatHandler::initClique() {
  // clique is not initialized to zero.
  // it works provided that assembly is done in the right order

  const int n_blocks = (ldc_ - 1) / nb_ + 1;
  (*clique_block_start_)[sn_].resize(n_blocks + 1);
  int schur_size{};
  for (int j = 0; j < n_blocks; ++j) {
    (*clique_block_start_)[sn_][j] = schur_size;
    const int jb = std::min(nb_, ldc_ - j * nb_);
    schur_size += (ldc_ - j * nb_) * jb;
  }
  (*clique_block_start_)[sn_].back() = schur_size;
  *clique_ = new double[schur_size];
}

void HybridHybridFormatHandler::assembleFrontal(int i, int j, double val) {
  (*frontal_)[i + j * ldf_ - j * (j + 1) / 2] = val;
}

void HybridHybridFormatHandler::assembleFrontalMultiple(int num, double* child,
                                                        int nc, int child_sn,
                                                        int row, int col, int i,
                                                        int j) {
  const int jblock = col / nb_;
  const int jb = std::min(nb_, nc - nb_ * jblock);
  const int row_ = row - jblock * nb_;
  const int col_ = col - jblock * nb_;
  const int start_block = (*clique_block_start_)[child_sn][jblock];
  daxpy_(&num, &d_one, &child[start_block + col_ + jb * row_], &jb,
         &(*frontal_)[i + ldf_ * j - j * (j + 1) / 2], &i_one);
}

int HybridHybridFormatHandler::denseFactorise(
    double reg_thresh, std::vector<double>& regularization,
    std::vector<double>& times) {
  int status;

  status = dense_fact_l2h(frontal_->data(), ldf_, sn_size_, nb_, times.data());
  if (status) return status;

  if (S_->factType() == FactType::Chol) {
    // find the position within regularization corresponding to this supernode
    int sn_start = S_->snStart(sn_);
    double* regul = &regularization[sn_start];

    status = dense_fact_pdbs(ldf_, sn_size_, nb_, frontal_->data(), *clique_,
                             reg_thresh, regul, times.data());
  } else {
    // find the position within pivot_sign corresponding to this supernode
    int sn_start = S_->snStart(sn_);
    const int* pivot_sign = &S_->pivotSign().data()[sn_start];
    double* regul = &regularization[sn_start];

    status =
        dense_fact_pibs(ldf_, sn_size_, S_->blockSize(), frontal_->data(),
                        *clique_, pivot_sign, reg_thresh, regul, times.data());
  }

  return status;
}

void HybridHybridFormatHandler::assembleClique(double* child, int nc,
                                               int child_sn) {
  // assemble the child clique into the current clique by blocks of columns.
  // within a block, assemble by rows.

  const int n_blocks = (nc - 1) / nb_ + 1;

  int row_start{};

  // go through the blocks of columns of the child sn
  for (int b = 0; b < n_blocks; ++b) {
    const int b_start = (*clique_block_start_)[child_sn][b];

    const int col_start = row_start;
    const int col_end = std::min(col_start + nb_, nc);

    // go through the rows within this block
    for (int row = row_start; row < nc; ++row) {
      const int i = S_->relindClique(child_sn, row) - sn_size_;

      // already assembled into frontal
      if (i < 0) continue;

      // go through the columns of the block
      int col = col_start;
      while (col < col_end) {
        int j = S_->relindClique(child_sn, col);
        if (j < sn_size_) {
          ++col;
          continue;
        }
        j -= sn_size_;

        // information and sizes of child sn
        const int jblock_c = b;
        const int jb_c = std::min(nb_, nc - nb_ * jblock_c);
        const int row_ = row - jblock_c * nb_;
        const int col_ = col - jblock_c * nb_;
        const int start_block_c = b_start;

        // sun consecutive entries in a row.
        // consecutive need to be reduced, to account for edge of the block
        const int zeros_stored_row = std::max(0, jb_c - (row - row_start) - 1);
        int consecutive = S_->consecutiveSums(child_sn, col);
        const int left_in_child = col_end - col - zeros_stored_row;
        consecutive = std::min(consecutive, left_in_child);

        // consecutive need to account also for edge of block in parent
        const int block_in_parent = j / nb_;
        const int col_end_parent = std::min((block_in_parent + 1) * nb_, ldc_);
        const int left_in_parent = col_end_parent - j;
        consecutive = std::min(consecutive, left_in_parent);

        // needed to deal with zeros stored in upper right part of block
        if (consecutive == 0) break;

        // information and sizes of current sn
        const int jblock = block_in_parent;
        const int jb = std::min(nb_, ldc_ - nb_ * jblock);
        const int i_ = i - jblock * nb_;
        const int j_ = j - jblock * nb_;
        const int start_block = (*clique_block_start_)[sn_][jblock];

        const double d_one = 1.0;
        const int i_one = 1;
        daxpy_(&consecutive, &d_one, &child[start_block_c + col_ + jb_c * row_],
               &i_one, &(*clique_)[start_block + j_ + jb * i_], &i_one);

        col += consecutive;
      }
    }

    row_start += nb_;
  }
}