#ifndef FACTOR_HIGHS_SETTINGS_H
#define FACTOR_HIGHS_SETTINGS_H

// ===========================================================================
// SWITCHES
// ===========================================================================

#define PARALLEL_TREE
// #define PARALLEN_NODE    // not active now
#define PIVOTING
#define DATA_COLLECTION
// #define PRINT_ITER_REF
// #define PRINT_REGULARIZATION

// ===========================================================================
// PARAMETERS
// ===========================================================================

// supernode amalgamation
const int kStartThreshRelax = 256;
const double kUpperRatioRelax = 0.02;
const double kLowerRatioRelax = 0.01;
const int kMaxIterRelax = 10;
const int kSnSizeRelax = 8;

// dense factorization
const int kBlockSize = 128;
const double kAlphaBK = (sqrt(17.0) + 1.0) / 8.0;

// regularization
const double kPrimalStaticRegularization = 1e-12;
const double kDualStaticRegularization = 1e-10;
const double kDynamicDiagCoeff = 1e-16;

// refinement
const int kMaxRefinementIter = 5;
const double kRefinementTolerance = 1e-8;

#endif