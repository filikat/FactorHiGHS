#include "FormatHandler.h"

FormatHandler::FormatHandler(const Symbolic& S, int sn)
    : S_{&S},
      sn_{sn},
      nb_{S_->blockSize()},
      sn_size_{S_->snStart(sn_ + 1) - S_->snStart(sn_)},
      ldf_{S_->ptr(sn_ + 1) - S_->ptr(sn_)},
      ldc_{ldf_ - sn_size_} {}

void FormatHandler::terminate(std::vector<double>& frontal,
                              std::vector<double>& clique) {
  // Move local copies of data into their final position.
  // In this way, the shared objects sn_columns_ and schur_contribution_ are
  // accessed only here, while a local copy is used for the assembly and dense
  // factorization. This should avoid the problem of false sharing.

  frontal = std::move(frontal_);
  clique = std::move(clique_);
}
