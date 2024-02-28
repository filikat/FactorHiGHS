#ifndef AUX_ANALYZE_H
#define AUX_ANALYZE_H

#include <vector>

void Counts2Ptr(std::vector<int>& ptr, std::vector<int>& w);
void InversePerm(const std::vector<int>& perm, std::vector<int>& iperm);
void PermuteVector(std::vector<int>& v, const std::vector<int>& perm);

#endif
