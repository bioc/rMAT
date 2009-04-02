#include <R.h>
#include <Rmath.h>
#include <float.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void MAT (double *C, double *I, int *nProbes, int *nArraysC, int *nArraysI, int *position, double *dMax, double *dMerge, int *nProbesMin, double *threshold, double *MATScores, double *pValues, int *method, int *regions, int *verbose, int *seqNum, int *numRegions);
void NormalizeProbes (char **seq, double *y, double *yNormalized, int *nProbes, int *nArrays, double *copyNumber, int *method, int *robust, double *adjRSquare, double *RSquare, double *BIC, double *outBeta, int *permVector, int *betaLength, int *all,  int *MATScaling,int *isVerbose);
SEXP Parser (SEXP fileName);
SEXP convertSeqToChNo (SEXP seqNum, SEXP cTableSeq, SEXP cTableChNo);
SEXP BPMAPCelMerger (SEXP list, SEXP celList);
SEXP WriteBAR (SEXP fileName, SEXP barNameString, SEXP MATScores, SEXP pValue, SEXP position, SEXP chNos, SEXP numSeq);
SEXP readBPMAPSeq (SEXP fileName, SEXP rangelist,  SEXP readPM, SEXP readMM, SEXP readProbeLength, SEXP readPMProbeSeq, SEXP readMatchScore, SEXP readPosition, SEXP readTopStrand, SEXP verbose, SEXP maxRange, SEXP minRange) ;
SEXP readBPMAPFileHeader (SEXP fileName);
SEXP readBPMAPSeqHeader (SEXP fileName, SEXP seqToRead);
SEXP readBPMAPAllSeqHeader (SEXP fileName);



R_CallMethodDef CallEntries[]  = {
      {"convertSeqToChNo", (DL_FUNC)&convertSeqToChNo, 4},
      {"Parser", (DL_FUNC)&Parser, 2},
      {"WriteBAR", (DL_FUNC)&WriteBAR,8},       
      {"readBPMAPFileHeader", (DL_FUNC)&readBPMAPFileHeader, 2},
      {"readBPMAPSeqHeader", (DL_FUNC)&readBPMAPAllSeqHeader, 3},
      {"readBPMAPAllSeqHeader", (DL_FUNC)&readBPMAPAllSeqHeader, 2},
      {"BPMAPCelMerger", (DL_FUNC)&BPMAPCelMerger, 3},
      {"readBPMAPSeq", (DL_FUNC)&readBPMAPSeq, 13},
      {NULL, NULL, 0}
    };

static const R_CMethodDef CEntries[] = {
  {"MAT", (DL_FUNC)&MAT, 18},
  {"NormalizeProbes", (DL_FUNC)&NormalizeProbes, 18},  
  {NULL, NULL, 0}
};

void R_init_MATR(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    // R_useDynamicSymbols(dll, FALSE);
}
