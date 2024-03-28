#ifndef AUXILIARY_H
#define AUXILIARY_H

#include <chrono>
#include <fstream>
#include <string>
#include <vector>

void Counts2Ptr(std::vector<int>& ptr, std::vector<int>& w);
void InversePerm(const std::vector<int>& perm, std::vector<int>& iperm);
void PermuteVector(std::vector<int>& v, const std::vector<int>& perm);
void SubtreeSize(const std::vector<int>& parent, std::vector<int>& sizes);
void Transpose(const std::vector<int>& ptr, const std::vector<int>& rows,
               std::vector<int>& ptrT, std::vector<int>& rowsT);
void Transpose(const std::vector<int>& ptr, const std::vector<int>& rows,
               const std::vector<double>& val, std::vector<int>& ptrT,
               std::vector<int>& rowsT, std::vector<double>& valT);
void ChildrenLinkedList(const std::vector<int>& parent, std::vector<int>& head,
                        std::vector<int>& next);

template <typename T>
void print(std::ofstream& out_file, const std::vector<T>& v,
           const std::string s) {
  char name[80];
  snprintf(name, 80, "matlab/%s.txt", s.c_str());
  out_file.open(name);
  for (T i : v) {
    out_file << i << '\n';
  }
  out_file.close();
}

class Clock {
  std::chrono::high_resolution_clock::time_point t0;

 public:
  void start();
  double stop();
};

#endif
