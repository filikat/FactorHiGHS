#include "FormatHandler.h"

void FormatHandler::init(const Symbolic* S) {
  S_ = S;
}

void FormatHandler::attach(int sn) {
  sn_ = sn;
  const int sn_begin = S_->snStart(sn_);
  const int sn_end = S_->snStart(sn_ + 1);
  sn_size_ = sn_end - sn_begin;
  ldf_ = S_->ptr(sn_ + 1) - S_->ptr(sn_);
  ldc_ = ldf_ - sn_size_;
  nb_ = S_->blockSize();

  // initialize frontal and clique
  initFrontal();
  initClique();
}

void FormatHandler::detach(std::vector<double>& frontal,
                           std::vector<double>& clique) {
  // move local copies of data into their final position
  frontal = std::move(frontal_);
  clique = std::move(clique_);
  sn_ = -1;
}
