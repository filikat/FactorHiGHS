#include "FormatHandler.h"

void FormatHandler::init(const Symbolic* S) {
  S_ = S;
  clique_block_start_.resize(S_->sn());
}

void FormatHandler::attach(std::vector<double>* frontal, int sn) {
  frontal_ = frontal;

  sn_ = sn;
  const int sn_begin = S_->snStart(sn_);
  const int sn_end = S_->snStart(sn_ + 1);
  sn_size_ = sn_end - sn_begin;
  ldf_ = S_->ptr(sn_ + 1) - S_->ptr(sn_);
  ldc_ = ldf_ - sn_size_;
  nb_ = S_->blockSize();
}

void FormatHandler::attachClique(double* clique) { clique_ = clique; }

void FormatHandler::detach() {
  frontal_ = nullptr;
  clique_ = nullptr;
}
