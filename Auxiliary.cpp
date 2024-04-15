#include "Auxiliary.h"

#include <stack>

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
  head.assign(n, -1);
  next.assign(n, -1);
  for (int node = n - 1; node >= 0; --node) {
    if (parent[node] == -1) continue;
    next[node] = head[parent[node]];
    head[parent[node]] = node;
  }
}

void DFS_post(int node, int& start, std::vector<int>& head,
              const std::vector<int>& next, std::vector<int>& order) {
  // Perform depth first search starting from root node and order the nodes
  // starting from the value start. head and next contain the linked list of
  // children.

  std::stack<int> stack;
  stack.push(node);

  while (!stack.empty()) {
    const int current = stack.top();
    const int child = head[current];

    if (child == -1) {
      // no children left to order,
      // remove from the stack and order
      stack.pop();
      order[start++] = current;
    } else {
      // at least one child left to order,
      // add it to the stack and remove it from the list of children
      stack.push(child);
      head[current] = next[child];
    }
  }
}

void ProcessEdge(int j, int i, const std::vector<int>& first,
                 std::vector<int>& maxfirst, std::vector<int>& delta,
                 std::vector<int>& prevleaf, std::vector<int>& ancestor) {
  // Process edge of skeleton matrix.
  // Taken from Tim Davis "Direct Methods for Sparse Linear Systems".

  // j not a leaf of ith row subtree
  if (i <= j || first[j] <= maxfirst[i]) {
    return;
  }

  // max first[j] so far
  maxfirst[i] = first[j];

  // previous leaf of ith row subtree
  int jprev = prevleaf[i];

  // A(i,j) is in the skeleton matrix
  delta[j]++;

  if (jprev != -1) {
    // find least common ancestor of jprev and j
    int q = jprev;
    while (q != ancestor[q]) {
      q = ancestor[q];
    }

    // path compression
    int sparent;
    for (int s = jprev; s != q; s = sparent) {
      sparent = ancestor[s];
      ancestor[s] = q;
    }

    // consider overlap
    delta[q]--;
  }

  // previous leaf of ith subtree set to j
  prevleaf[i] = j;
}

void Clock::start() { t0 = std::chrono::high_resolution_clock::now(); }
double Clock::stop() {
  auto t1 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> d = t1 - t0;
  return d.count();
}
