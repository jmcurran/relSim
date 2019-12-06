// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// breed
IntegerVector breed(IntegerVector Parents, int ns, int Ns, int nLoci);
RcppExport SEXP _relSim_breed(SEXP ParentsSEXP, SEXP nsSEXP, SEXP NsSEXP, SEXP nLociSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type Parents(ParentsSEXP);
    Rcpp::traits::input_parameter< int >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< int >::type Ns(NsSEXP);
    Rcpp::traits::input_parameter< int >::type nLoci(nLociSEXP);
    rcpp_result_gen = Rcpp::wrap(breed(Parents, ns, Ns, nLoci));
    return rcpp_result_gen;
END_RCPP
}
// calcFst
NumericVector calcFst(const IntegerVector& Pop, IntegerVector SubPopIdx, int N, int ns, int nLoci, IntegerVector NumLocusAlleles);
RcppExport SEXP _relSim_calcFst(SEXP PopSEXP, SEXP SubPopIdxSEXP, SEXP NSEXP, SEXP nsSEXP, SEXP nLociSEXP, SEXP NumLocusAllelesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type Pop(PopSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type SubPopIdx(SubPopIdxSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< int >::type nLoci(nLociSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type NumLocusAlleles(NumLocusAllelesSEXP);
    rcpp_result_gen = Rcpp::wrap(calcFst(Pop, SubPopIdx, N, ns, nLoci, NumLocusAlleles));
    return rcpp_result_gen;
END_RCPP
}
// calcFStatistics
List calcFStatistics(const IntegerVector& Pop, IntegerVector SubPopIdx, int N, int ns, int nLoci, IntegerVector NumLocusAlleles);
RcppExport SEXP _relSim_calcFStatistics(SEXP PopSEXP, SEXP SubPopIdxSEXP, SEXP NSEXP, SEXP nsSEXP, SEXP nLociSEXP, SEXP NumLocusAllelesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type Pop(PopSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type SubPopIdx(SubPopIdxSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type ns(nsSEXP);
    Rcpp::traits::input_parameter< int >::type nLoci(nLociSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type NumLocusAlleles(NumLocusAllelesSEXP);
    rcpp_result_gen = Rcpp::wrap(calcFStatistics(Pop, SubPopIdx, N, ns, nLoci, NumLocusAlleles));
    return rcpp_result_gen;
END_RCPP
}
// IS
List IS(List freqs, int N, int numContributors, int maxAllelesShowing, List Perms, bool bTail);
RcppExport SEXP _relSim_IS(SEXP freqsSEXP, SEXP NSEXP, SEXP numContributorsSEXP, SEXP maxAllelesShowingSEXP, SEXP PermsSEXP, SEXP bTailSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type freqs(freqsSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type numContributors(numContributorsSEXP);
    Rcpp::traits::input_parameter< int >::type maxAllelesShowing(maxAllelesShowingSEXP);
    Rcpp::traits::input_parameter< List >::type Perms(PermsSEXP);
    Rcpp::traits::input_parameter< bool >::type bTail(bTailSEXP);
    rcpp_result_gen = Rcpp::wrap(IS(freqs, N, numContributors, maxAllelesShowing, Perms, bTail));
    return rcpp_result_gen;
END_RCPP
}
// locusLRmix_Caller
double locusLRmix_Caller(IntegerVector ProfVic, IntegerVector ProfSus, NumericVector Freq);
RcppExport SEXP _relSim_locusLRmix_Caller(SEXP ProfVicSEXP, SEXP ProfSusSEXP, SEXP FreqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type ProfVic(ProfVicSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ProfSus(ProfSusSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Freq(FreqSEXP);
    rcpp_result_gen = Rcpp::wrap(locusLRmix_Caller(ProfVic, ProfSus, Freq));
    return rcpp_result_gen;
END_RCPP
}
// LRmix
NumericVector LRmix(IntegerVector ProfVic, IntegerVector ProfSus, List listFreqs);
RcppExport SEXP _relSim_LRmix(SEXP ProfVicSEXP, SEXP ProfSusSEXP, SEXP listFreqsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type ProfVic(ProfVicSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ProfSus(ProfSusSEXP);
    Rcpp::traits::input_parameter< List >::type listFreqs(listFreqsSEXP);
    rcpp_result_gen = Rcpp::wrap(LRmix(ProfVic, ProfSus, listFreqs));
    return rcpp_result_gen;
END_RCPP
}
// locusIBS
IntegerVector locusIBS(IntegerVector ProfMat, int N);
RcppExport SEXP _relSim_locusIBS(SEXP ProfMatSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type ProfMat(ProfMatSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(locusIBS(ProfMat, N));
    return rcpp_result_gen;
END_RCPP
}
// IBS_Caller
int IBS_Caller(IntegerVector Prof1, IntegerVector Prof2, int nLoci);
RcppExport SEXP _relSim_IBS_Caller(SEXP Prof1SEXP, SEXP Prof2SEXP, SEXP nLociSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type Prof1(Prof1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Prof2(Prof2SEXP);
    Rcpp::traits::input_parameter< int >::type nLoci(nLociSEXP);
    rcpp_result_gen = Rcpp::wrap(IBS_Caller(Prof1, Prof2, nLoci));
    return rcpp_result_gen;
END_RCPP
}
// randomProfiles
IntegerVector randomProfiles(List listFreqs, int nBlockSize);
RcppExport SEXP _relSim_randomProfiles(SEXP listFreqsSEXP, SEXP nBlockSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type listFreqs(listFreqsSEXP);
    Rcpp::traits::input_parameter< int >::type nBlockSize(nBlockSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(randomProfiles(listFreqs, nBlockSize));
    return rcpp_result_gen;
END_RCPP
}
// randomSibs
IntegerVector randomSibs(IntegerVector ProfSib1, List listFreqs, int nBlockSize);
RcppExport SEXP _relSim_randomSibs(SEXP ProfSib1SEXP, SEXP listFreqsSEXP, SEXP nBlockSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type ProfSib1(ProfSib1SEXP);
    Rcpp::traits::input_parameter< List >::type listFreqs(listFreqsSEXP);
    Rcpp::traits::input_parameter< int >::type nBlockSize(nBlockSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(randomSibs(ProfSib1, listFreqs, nBlockSize));
    return rcpp_result_gen;
END_RCPP
}
// randomChildren
IntegerVector randomChildren(IntegerVector ProfParent, List listFreqs, int nBlockSize);
RcppExport SEXP _relSim_randomChildren(SEXP ProfParentSEXP, SEXP listFreqsSEXP, SEXP nBlockSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type ProfParent(ProfParentSEXP);
    Rcpp::traits::input_parameter< List >::type listFreqs(listFreqsSEXP);
    Rcpp::traits::input_parameter< int >::type nBlockSize(nBlockSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(randomChildren(ProfParent, listFreqs, nBlockSize));
    return rcpp_result_gen;
END_RCPP
}
// lrPC_Caller
double lrPC_Caller(IntegerVector ProfParent, IntegerVector ProfChild, List listFreqs);
RcppExport SEXP _relSim_lrPC_Caller(SEXP ProfParentSEXP, SEXP ProfChildSEXP, SEXP listFreqsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type ProfParent(ProfParentSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ProfChild(ProfChildSEXP);
    Rcpp::traits::input_parameter< List >::type listFreqs(listFreqsSEXP);
    rcpp_result_gen = Rcpp::wrap(lrPC_Caller(ProfParent, ProfChild, listFreqs));
    return rcpp_result_gen;
END_RCPP
}
// maximizeLRPC
List maximizeLRPC(List listFreqs, int nBlockSize);
RcppExport SEXP _relSim_maximizeLRPC(SEXP listFreqsSEXP, SEXP nBlockSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type listFreqs(listFreqsSEXP);
    Rcpp::traits::input_parameter< int >::type nBlockSize(nBlockSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(maximizeLRPC(listFreqs, nBlockSize));
    return rcpp_result_gen;
END_RCPP
}
// lrSib_Caller
double lrSib_Caller(IntegerVector ProfSib1, IntegerVector ProfSib2, List listFreqs);
RcppExport SEXP _relSim_lrSib_Caller(SEXP ProfSib1SEXP, SEXP ProfSib2SEXP, SEXP listFreqsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type ProfSib1(ProfSib1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ProfSib2(ProfSib2SEXP);
    Rcpp::traits::input_parameter< List >::type listFreqs(listFreqsSEXP);
    rcpp_result_gen = Rcpp::wrap(lrSib_Caller(ProfSib1, ProfSib2, listFreqs));
    return rcpp_result_gen;
END_RCPP
}
// blockStatCounts
IntegerVector blockStatCounts(IntegerVector Prof1, IntegerVector Prof2, int nProf, List listFreqs, int nCode, bool bFalseNeg, IntegerVector IBSthresh, NumericVector LRthresh, int nNumResults);
RcppExport SEXP _relSim_blockStatCounts(SEXP Prof1SEXP, SEXP Prof2SEXP, SEXP nProfSEXP, SEXP listFreqsSEXP, SEXP nCodeSEXP, SEXP bFalseNegSEXP, SEXP IBSthreshSEXP, SEXP LRthreshSEXP, SEXP nNumResultsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type Prof1(Prof1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Prof2(Prof2SEXP);
    Rcpp::traits::input_parameter< int >::type nProf(nProfSEXP);
    Rcpp::traits::input_parameter< List >::type listFreqs(listFreqsSEXP);
    Rcpp::traits::input_parameter< int >::type nCode(nCodeSEXP);
    Rcpp::traits::input_parameter< bool >::type bFalseNeg(bFalseNegSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type IBSthresh(IBSthreshSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type LRthresh(LRthreshSEXP);
    Rcpp::traits::input_parameter< int >::type nNumResults(nNumResultsSEXP);
    rcpp_result_gen = Rcpp::wrap(blockStatCounts(Prof1, Prof2, nProf, listFreqs, nCode, bFalseNeg, IBSthresh, LRthresh, nNumResults));
    return rcpp_result_gen;
END_RCPP
}
// blockStats
List blockStats(IntegerVector Prof1, IntegerVector Prof2, int nProf, List listFreqs, int nCode);
RcppExport SEXP _relSim_blockStats(SEXP Prof1SEXP, SEXP Prof2SEXP, SEXP nProfSEXP, SEXP listFreqsSEXP, SEXP nCodeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type Prof1(Prof1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Prof2(Prof2SEXP);
    Rcpp::traits::input_parameter< int >::type nProf(nProfSEXP);
    Rcpp::traits::input_parameter< List >::type listFreqs(listFreqsSEXP);
    Rcpp::traits::input_parameter< int >::type nCode(nCodeSEXP);
    rcpp_result_gen = Rcpp::wrap(blockStats(Prof1, Prof2, nProf, listFreqs, nCode));
    return rcpp_result_gen;
END_RCPP
}
// simNpersonMixture
NumericMatrix simNpersonMixture(List listFreqs, int numContributors, int numIterations);
RcppExport SEXP _relSim_simNpersonMixture(SEXP listFreqsSEXP, SEXP numContributorsSEXP, SEXP numIterationsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type listFreqs(listFreqsSEXP);
    Rcpp::traits::input_parameter< int >::type numContributors(numContributorsSEXP);
    Rcpp::traits::input_parameter< int >::type numIterations(numIterationsSEXP);
    rcpp_result_gen = Rcpp::wrap(simNpersonMixture(listFreqs, numContributors, numIterations));
    return rcpp_result_gen;
END_RCPP
}
// famSearch
List famSearch(IntegerVector& profiles, IntegerVector& siblings, IntegerVector& children, List& listFreqs, int step);
RcppExport SEXP _relSim_famSearch(SEXP profilesSEXP, SEXP siblingsSEXP, SEXP childrenSEXP, SEXP listFreqsSEXP, SEXP stepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector& >::type profiles(profilesSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type siblings(siblingsSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type children(childrenSEXP);
    Rcpp::traits::input_parameter< List& >::type listFreqs(listFreqsSEXP);
    Rcpp::traits::input_parameter< int >::type step(stepSEXP);
    rcpp_result_gen = Rcpp::wrap(famSearch(profiles, siblings, children, listFreqs, step));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_relSim_breed", (DL_FUNC) &_relSim_breed, 4},
    {"_relSim_calcFst", (DL_FUNC) &_relSim_calcFst, 6},
    {"_relSim_calcFStatistics", (DL_FUNC) &_relSim_calcFStatistics, 6},
    {"_relSim_IS", (DL_FUNC) &_relSim_IS, 6},
    {"_relSim_locusLRmix_Caller", (DL_FUNC) &_relSim_locusLRmix_Caller, 3},
    {"_relSim_LRmix", (DL_FUNC) &_relSim_LRmix, 3},
    {"_relSim_locusIBS", (DL_FUNC) &_relSim_locusIBS, 2},
    {"_relSim_IBS_Caller", (DL_FUNC) &_relSim_IBS_Caller, 3},
    {"_relSim_randomProfiles", (DL_FUNC) &_relSim_randomProfiles, 2},
    {"_relSim_randomSibs", (DL_FUNC) &_relSim_randomSibs, 3},
    {"_relSim_randomChildren", (DL_FUNC) &_relSim_randomChildren, 3},
    {"_relSim_lrPC_Caller", (DL_FUNC) &_relSim_lrPC_Caller, 3},
    {"_relSim_maximizeLRPC", (DL_FUNC) &_relSim_maximizeLRPC, 2},
    {"_relSim_lrSib_Caller", (DL_FUNC) &_relSim_lrSib_Caller, 3},
    {"_relSim_blockStatCounts", (DL_FUNC) &_relSim_blockStatCounts, 9},
    {"_relSim_blockStats", (DL_FUNC) &_relSim_blockStats, 5},
    {"_relSim_simNpersonMixture", (DL_FUNC) &_relSim_simNpersonMixture, 3},
    {"_relSim_famSearch", (DL_FUNC) &_relSim_famSearch, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_relSim(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
