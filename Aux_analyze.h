#ifndef AUX_ANALYZE_H
#define AUX_ANALYZE_H

#include <vector>

void ColCount2Ptr(std::vector<int>& ptr, std::vector<int>& w);
void Inverse_perm(const std::vector<int>& perm, std::vector<int>& iperm);
void Permute_vector(std::vector<int>& v, const std::vector<int>& perm);

#endif
