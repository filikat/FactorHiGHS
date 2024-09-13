#include "FormatHandler.h"

void FormatHandler::attach(std::vector<double>* frontal, double** clique,
                           std::vector<std::vector<int>>* clique_block_start,
                           const Symbolic* S, int sn) {
  frontal_ = frontal;
  clique_ = clique;
  clique_block_start_ = clique_block_start;
  S_ = S;

  sn_ = sn;
  const int sn_begin = S_->snStart(sn_);
  const int sn_end = S_->snStart(sn_ + 1);
  sn_size_ = sn_end - sn_begin;
  ldf_ = S_->ptr(sn_ + 1) - S_->ptr(sn_);
  ldc_ = ldf_ - sn_size_;
  nb_ = S_->blockSize();
}

void FormatHandler::detach() {
  frontal_ = nullptr;
  clique_ = nullptr;
  clique_block_start_ = nullptr;
  S_ = nullptr;
}
