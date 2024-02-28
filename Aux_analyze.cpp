#include "Aux_analyze.h"

void Counts2Ptr(std::vector<int>& ptr, std::vector<int>& w) {
  // Given the column counts in the vector w (of size n),
  // compute the column pointers in the vector ptr (of size n+1),
  // and copy the first n pointers back into w.

  int temp_nz{};
  int n = w.size();
  for (int j = 0; j < n; ++j) {
    ptr[j] = temp_nz;
    temp_nz += w[j];
    w[j] = ptr[j];
  }
  ptr[n] = temp_nz;
}

void InversePerm(const std::vector<int>& perm, std::vector<int>& iperm) {
  // Given the permutation perm, produce the inverse permutation iperm.
  // perm[i] : i-th entry to use in the new order.
  // iperm[i]: where entry i is located in the new order.

  for (int i = 0; i < perm.size(); ++i) {
    iperm[perm[i]] = i;
  }
}

void PermuteVector(std::vector<int>& v, const std::vector<int>& perm) {
  // Permute vector v according to permutation perm.

  std::vector<int> new_v(v.size());
  for (int i = 0; i < v.size(); ++i) {
    new_v[i] = v[perm[i]];
  }
  v = std::move(new_v);
}