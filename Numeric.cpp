#include "Numeric.h"

#include "Auxiliary.h"
#include "CallAndTimeBlas.h"
#include "FormatHandler.h"
#include "FullSolveHandler.h"
#include "HybridSolveHandler.h"
#include "PackedSolveHandler.h"
#include "Timing.h"
#include "DataCollector.h"

Numeric::Numeric(const Symbolic& S) : S_{S} {
  // initialize solve handler
  switch (S_.formatType()) {
    case FormatType::Full:
      SH_.reset(new FullSolveHandler(S_, sn_columns_));
      break;
    case FormatType::HybridPacked:
    case FormatType::HybridHybrid:
      SH_.reset(new HybridSolveHandler(S_, sn_columns_, swaps_, pivot_2x2_));
      break;
    case FormatType::PackedPacked:
      SH_.reset(new PackedSolveHandler(S_, sn_columns_));
      break;
  }
}

void Numeric::solve(std::vector<double>& x) const {
#ifdef COARSE_TIMING
  Clock clock{};
  clock.start();
#endif

  // permute rhs
  permuteVectorInverse(x, S_.iperm());

  // solve
  SH_->forwardSolve(x);
  SH_->diagSolve(x);
  SH_->backwardSolve(x);

  // unpermute solution
  permuteVector(x, S_.iperm());

#ifdef COARSE_TIMING
  DataCollector::get()->sumTime(kTimeSolve, clock.stop());
#endif
}
