#include "SolveHandler.h"

SolveHandler::SolveHandler(const Symbolic& S, DataCollector& DC,
                           const std::vector<std::vector<double>>& sn_columns)
    : S_{S}, DC_{DC}, sn_columns_{sn_columns} {}