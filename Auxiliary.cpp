#include "Auxiliary.h"

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

void SubtreeSize(const std::vector<int>& parent, std::vector<int>& sizes) {
  // Compute sizes of subtrees of the tree given by parent

  int n = parent.size();
  sizes.assign(n, 1);

  for (int i = 0; i < n; ++i) {
    int k = parent[i];
    if (k != -1) sizes[k] += sizes[i];
  }
}

void Transpose(const std::vector<int>& ptr, const std::vector<int>& rows,
               std::vector<int>& ptrT, std::vector<int>& rowsT) {
  // Compute the transpose of the matrix and return it in rowsT and ptrT

  int n = ptr.size() - 1;

  std::vector<int> work(n);

  // count the entries in each row into work
  for (int i = 0; i < ptr.back(); ++i) {
    ++work[rows[i]];
  }

  // sum row sums to obtain pointers
  Counts2Ptr(ptrT, work);

  for (int j = 0; j < n; ++j) {
    for (int el = ptr[j]; el < ptr[j + 1]; ++el) {
      int i = rows[el];

      // entry (i,j) becomes entry (j,i)
      int pos = work[i]++;
      rowsT[pos] = j;
    }
  }
}

void Transpose(const std::vector<int>& ptr, const std::vector<int>& rows,
               const std::vector<double>& val, std::vector<int>& ptrT,
               std::vector<int>& rowsT, std::vector<double>& valT) {
  // Compute the transpose of the matrix and return it in rowsT, ptrT and valT

  int n = ptr.size() - 1;

  std::vector<int> work(n);

  // count the entries in each row into work
  for (int i = 0; i < ptr.back(); ++i) {
    ++work[rows[i]];
  }

  // sum row sums to obtain pointers
  Counts2Ptr(ptrT, work);

  for (int j = 0; j < n; ++j) {
    for (int el = ptr[j]; el < ptr[j + 1]; ++el) {
      int i = rows[el];

      // entry (i,j) becomes entry (j,i)
      int pos = work[i]++;
      rowsT[pos] = j;
      valT[pos] = val[el];
    }
  }
}

void ChildrenLinkedList(const std::vector<int>& parent, std::vector<int>& head,
                        std::vector<int>& next) {
  // Create linked lists of children in elimination tree.
  // parent gives the dependencies of the tree,
  // head[node] is the first child of node,
  // next[head[node]] is the second child,
  // next[next[head[node]]] is the third child...
  // until -1 is reached.

  int n = parent.size();
  head.resize(n, -1);
  next.resize(n, -1);
  for (int node = n - 1; node >= 0; --node) {
    if (parent[node] == -1) continue;
    next[node] = head[parent[node]];
    head[parent[node]] = node;
  }
}

void Clock::start() { t0 = std::chrono::high_resolution_clock::now(); }
double Clock::stop() {
  auto t1 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> d = t1 - t0;
  return d.count();
}