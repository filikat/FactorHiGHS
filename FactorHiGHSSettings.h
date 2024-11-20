#ifndef FACTOR_HIGHS_SETTINGS_H
#define FACTOR_HIGHS_SETTINGS_H

// switch for tree parallelism
#define PARALLEL_TREE

// switch for node parallelism (not active now)
// #define PARALLEN_NODE

// switch for pivoting
#define PIVOTING

// parameters for supernode amalgamation
const int kStartThreshRelax = 256;
const double kUpperRatioRelax = 0.02;
const double kLowerRatioRelax = 0.01;
const int kMaxIterRelax = 10;
const int kSnSizeRelax = 8;

// parameters for dense factorization
const int kBlockSize = 128;

#endif