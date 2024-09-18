#include "FullFormatHandler.h"

void FullFormatHandler::initFrontal() {
  // frontal is initialized to zero
  frontal_->resize(ldf_ * sn_size_);
}

void FullFormatHandler::initClique() {
  // clique is not initialized to zero.
  // it works provided that assembly is done in the right order
  *clique_ = new double[ldc_ * ldc_];
}

void FullFormatHandler::assembleFrontal(int i, int j, double val) {
  (*frontal_)[i + j * ldf_] = val;
}

void FullFormatHandler::assembleFrontalMultiple(int num, double* child, int nc,
                                                int child_sn, int row, int col,
                                                int i, int j) {
  daxpy_(&num, &d_one, &child[row + nc * col], &i_one,
         &(*frontal_)[i + ldf_ * j], &i_one);
}

int FullFormatHandler::denseFactorise(double reg_thresh,
                                      std::vector<double>& times) {
  int status;
  if (S_->type() == FactType::NormEq) {
    status = dense_fact_pdbf(ldf_, sn_size_, nb_, frontal_->data(), ldf_,
                             *clique_, ldc_, times.data(),reg_thresh);
  } else {
    status =
        dense_fact_pibf(ldf_, sn_size_, nb_, frontal_->data(), ldf_, *clique_,
                        ldc_, S_->pivotSign().data(), times.data());
  }
  return status;
}

void FullFormatHandler::assembleClique(double* child, int nc, int child_sn) {
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
        daxpy_(&consecutive, &d_one, &child[row + nc * col], &i_one,
               &(*clique_)[i + ldc_ * j], &i_one);

        row += consecutive;
      }
    }

    // j < sn_size was already done before, because it was needed before the
    // partial factorisation. Assembling into the clique instead can be done
    // after.
  }
}