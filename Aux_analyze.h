#ifndef AUX_ANALYZE_H
#define AUX_ANALYZE_H

#include <vector>
#include <fstream>
#include <string>

void Counts2Ptr(std::vector<int>& ptr, std::vector<int>& w);
void InversePerm(const std::vector<int>& perm, std::vector<int>& iperm);
void PermuteVector(std::vector<int>& v, const std::vector<int>& perm);
void SubtreeSize(const std::vector<int>& parent, std::vector<int>& sizes);

void print(std::ofstream& out_file, const std::vector<int>& v,
           const std::string s);

#endif
