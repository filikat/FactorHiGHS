#include "PackedPackedFormatHandler.h"

PackedPackedFormatHandler::PackedPackedFormatHandler(const Symbolic& S,
                                                     DataCollector& DC, int sn)
    : FormatHandler(S, DC, sn) {
  // initialize frontal and clique
  initFrontal();
  initClique();
}

void PackedPackedFormatHandler::initFrontal() {
  const int n_blocks = (sn_size_ - 1) / nb_ + 1;
  diag_start_.resize(n_blocks);
  int frontal_size = getDiagStart(ldf_, sn_size_, nb_, n_blocks, diag_start_);
  frontal_.resize(frontal_size + 10);
  // NB: the plus 10 is not needed, but it avoids weird problems later on.
}

void PackedPackedFormatHandler::initClique() {
  clique_.resize(S_->cliqueSize(sn_));
}

void PackedPackedFormatHandler::assembleFrontal(int i, int j, double val) {
  int block = j / nb_;
  int ldb = ldf_ - block * nb_;
  int ii = i - block * nb_;
  int jj = j - block * nb_;
  frontal_[diag_start_[block] + ii + ldb * jj] = val;
}

void PackedPackedFormatHandler::assembleFrontalMultiple(
    int num, const std::vector<double>& child, int nc, int child_sn, int row,
    int col, int i, int j) {
  const int jblock = col / nb_;
  row -= jblock * nb_;
  col -= jblock * nb_;
  const int start_block = S_->cliqueBlockStart(child_sn, jblock);
  const int ld = nc - nb_ * jblock;

  int block = j / nb_;
  int ldb = ldf_ - block * nb_;
  int ii = i - block * nb_;
  int jj = j - block * nb_;

  daxpy_(&num, &d_one, &child[start_block + row + ld * col], &i_one,
         &frontal_[diag_start_[block] + ii + ldb * jj], &i_one);
}

int PackedPackedFormatHandler::denseFactorise(double reg_thresh) {
  // find the position within pivot_sign corresponding to this supernode
  int sn_start = S_->snStart(sn_);
  const int* pivot_sign = &S_->pivotSign().data()[sn_start];

  int status = denseFactFP(ldf_, sn_size_, nb_, frontal_.data(), clique_.data(),
                           pivot_sign, reg_thresh, local_reg_.data(), DC_);

  return status;
}

void PackedPackedFormatHandler::assembleClique(const std::vector<double>& child,
                                               int nc, int child_sn) {
  //   go through the columns of the contribution of the child
  for (int col = 0; col < nc; ++col) {
    // relative index of column in the frontal matrix
    int j = S_->relindClique(child_sn, col);

    if (j >= sn_size_) {
      // assemble into clique

      // adjust relative index to access clique
      j -= sn_size_;

      // go through the rows of the contribution of the child
      int row = col;
      while (row < nc) {
        // relative index of the entry in the matrix clique
        const int i = S_->relindClique(child_sn, row) - sn_size_;

        // how many entries to sum
        const int consecutive = S_->consecutiveSums(child_sn, row);

        // use daxpy_ for summing consecutive entries

        const int jblock_c = col / nb_;
        const int jb_c = std::min(nb_, nc - nb_ * jblock_c);
        const int row_c = row - jblock_c * nb_;
        const int col_c = col - jblock_c * nb_;
        const int start_block_c = S_->cliqueBlockStart(child_sn, jblock_c);
        const int ld_c = nc - nb_ * jblock_c;

        const int jblock = j / nb_;
        const int jb = std::min(nb_, ldc_ - nb_ * jblock);
        const int ii = i - jblock * nb_;
        const int jj = j - jblock * nb_;
        const int start_block = S_->cliqueBlockStart(sn_, jblock);
        const int ld = ldc_ - nb_ * jblock;

        daxpy_(&consecutive, &d_one,
               &child[start_block_c + row_c + ld_c * col_c], &i_one,
               &clique_[start_block + ii + ld * jj], &i_one);

        row += consecutive;
      }
    }
  }
}

void PackedPackedFormatHandler::extremeEntries() {
  double minD = std::numeric_limits<double>::max();
  double maxD = 0.0;
  double minoffD = std::numeric_limits<double>::max();
  double maxoffD = 0.0;

  // number of blocks of columns
  const int n_blocks = (sn_size_ - 1) / nb_ + 1;

  // index to access frontal
  int index{};

  // go through blocks of columns for this supernode
  for (int j = 0; j < n_blocks; ++j) {
    // number of columns in the block
    const int jb = std::min(nb_, sn_size_ - nb_ * j);

    for (int k = 0; k < jb; ++k) {
      // off diagonal entries
      for (int i = 0; i < k; ++i) {
        if (frontal_[index] != 0.0) {
          minoffD = std::min(minoffD, std::abs(frontal_[index]));
          maxoffD = std::max(maxoffD, std::abs(frontal_[index]));
        }
        index++;
      }

      // diagonal entry
      minD = std::min(minD, std::abs(frontal_[index]));
      maxD = std::max(maxD, std::abs(frontal_[index]));
      index++;
    }

    const int entries_left = (ldf_ - nb_ * j - jb) * jb;

    for (int i = 0; i < entries_left; ++i) {
      if (frontal_[index] != 0.0) {
        minoffD = std::min(minoffD, std::abs(frontal_[index]));
        maxoffD = std::max(maxoffD, std::abs(frontal_[index]));
      }
      index++;
    }
  }

  DC_.extremeEntries(minD, maxD, minoffD, maxoffD);
}