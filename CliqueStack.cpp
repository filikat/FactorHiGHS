#include "CliqueStack.h"

void CliqueStack::init(int stack_size, int max_clique_size) {
  stack_.reserve(stack_size);
  top_ = 0;
  work_.reserve(max_clique_size);
}

double* CliqueStack::setup(int clique_size) {
  // Clear workspace

  // This should not trigger reallocation, because the reserve in init is done
  // with the maximum possible size of a clique.
  if (clique_size > work_.capacity())
    printf("=== Warning, reallocation of workspace\n");
  work_.assign(clique_size, 0.0);

  worksize_ = clique_size;
  return work_.data();
}

const double* CliqueStack::getChild(int& child_sn) const {
  // Get the top of the stack, in terms of supernode ID of the child and pointer
  // to its data.

  child_sn = sn_pushed_.top().first;
  int child_size = sn_pushed_.top().second;
  const double* child = &stack_[top_ - child_size];

  return child;
}

void CliqueStack::popChild() {
  // Remove top child from the stack

  int child_size = sn_pushed_.top().second;
  sn_pushed_.pop();

  top_ -= child_size;
}

void CliqueStack::pushWork(int sn) {
  // Put the content of the workspace at the top of the stack

  // This should not trigger reallocation, because the reserve in init is done
  // with the maximum possible size of the stack.
  if (top_ + worksize_ > stack_.capacity())
    printf("=== Warning, reallocation of stack\n");
  stack_.resize(top_ + worksize_);

  // copy workspace into top
  std::memcpy(&stack_[top_], work_.data(), worksize_ * sizeof(double));

  // update top
  top_ += worksize_;

  // keep track of supernodes pushed
  sn_pushed_.push({sn, worksize_});

  worksize_ = 0;
}