#ifndef TIMING_H
#define TIMING_H

// choose level of timing: 1, 2, 3
#define TIMING_2

// define for timing
#if (defined(TIMING_1) || defined(TIMING_2) || defined(TIMING_3))
#define COARSE_TIMING
#endif

#if (defined(TIMING_2) || defined(TIMING_3))
#define FINE_TIMING
#endif

#if (defined(TIMING_3))
#define FINEST_TIMING
#endif

enum TimeItems {
  // Analyse timer
  kTimeAnalyse,  // TIMING_1
  // Analyse items
  kTimeAnalyseMetis,    // TIMING_2
  kTimeAnalyseTree,     // TIMING_2
  kTimeAnalyseCount,    // TIMING_2
  kTimeAnalysePattern,  // TIMING_2
  kTimeAnalyseSn,       // TIMING_2
  kTimeAnalyseReorder,  // TIMING_2
  kTimeAnalyseRelInd,   // TIMING_2
  kTimeAnalyseLayer0,   // TIMING_2
  // Factorise timer
  kTimeFactorise,  // TIMING_1
  // Factorise items
  kTimeFactorisePrepare,            // TIMING_2
  kTimeFactoriseAssembleOriginal,   // TIMING_2
  kTimeFactoriseAssembleChildrenF,  // TIMING_2
  kTimeFactoriseAssembleChildrenC,  // TIMING_2
  kTimeFactoriseDenseFact,          // TIMING_2
  // DenseFact items
  kTimeDenseFact_trsm,       // TIMING_3
  kTimeDenseFact_syrk,       // TIMING_3
  kTimeDenseFact_gemm,       // TIMING_3
  kTimeDenseFact_fact,       // TIMING_3
  kTimeDenseFact_copy,       // TIMING_3
  kTimeDenseFact_copyschur,  // TIMING_3
  kTimeDenseFact_scal,       // TIMING_3
  kTimeDenseFact_convert,    // TIMING_3
  // enum size
  kTimeSize
};

double GetTime();

#endif