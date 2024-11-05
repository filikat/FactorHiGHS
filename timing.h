#ifndef TIMING_H
#define TIMING_H

// choose level of timing:
// - TIMING_0: no timing
// - TIMING_1: basic timing
// - TIMING_2: advanced timing
// - TIMING_3: extreme timing (timing of each BLAS call, considerably slower)
#define TIMING_2

// define for timing
#if (defined(TIMING_1) || defined(TIMING_2) || defined(TIMING_3))
#define COARSE_TIMING
#endif

#if (defined(TIMING_2) || defined(TIMING_3))
#define FINE_TIMING
#endif

#if (defined(TIMING_3))
#define BLAS_TIMING
#endif

enum TimeItems {
  // Analyse timer
  kTimeAnalyse,         // TIMING_1
  kTimeAnalyseMetis,    // TIMING_2
  kTimeAnalyseTree,     // TIMING_2
  kTimeAnalyseCount,    // TIMING_2
  kTimeAnalysePattern,  // TIMING_2
  kTimeAnalyseSn,       // TIMING_2
  kTimeAnalyseReorder,  // TIMING_2
  kTimeAnalyseRelInd,   // TIMING_2
  // Factorise timer
  kTimeFactorise,                         // TIMING_1
  kTimeFactorisePrepare,                  // TIMING_2
  kTimeFactoriseAssembleOriginal,         // TIMING_2
  kTimeFactoriseAssembleChildrenFrontal,  // TIMING_2
  kTimeFactoriseAssembleChildrenClique,   // TIMING_2
  kTimeFactoriseDenseFact,                // TIMING_2
  kTimeDenseFact_fact,                    // TIMING_2
  kTimeDenseFact_convert,                 // TIMING_2
  kTimeFactoriseTerminate,                // TIMING_2
  // Solve timer
  kTimeSolve,  // TIMING_1
  // BLAS times
  kTimeBlasStart,
  kTimeBlas_copy = kTimeBlasStart,  // TIMING_3
  kTimeBlas_axpy,                   // TIMING_3
  kTimeBlas_scal,                   // TIMING_3
  kTimeBlas_gemv,                   // TIMING_3
  kTimeBlas_trsv,                   // TIMING_3
  kTimeBlas_tpsv,                   // TIMING_3
  kTimeBlas_trsm,                   // TIMING_3
  kTimeBlas_syrk,                   // TIMING_3
  kTimeBlas_gemm,                   // TIMING_3
  kTimeBlasEnd = kTimeBlas_gemm,
  // enum size
  kTimeSize
};

#endif