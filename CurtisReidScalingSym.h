#ifndef CURTIS_REID_SCALING_SYM_H
#define CURTIS_REID_SCALING_SYM_H

#include <cmath>
#include <vector>

void CurtisReidScalingSym(const std::vector<int>& ptr,
                          const std::vector<int>& rows,
                          const std::vector<double>& val,
                          std::vector<int>& colexp);

#endif