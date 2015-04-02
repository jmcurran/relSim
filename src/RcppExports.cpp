// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP relSim_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    __result = Rcpp::wrap(rcpp_hello_world());
    return __result;
END_RCPP
}
// randomProfilesC
IntegerVector randomProfilesC(List listFreqs, int nBlockSize);
RcppExport SEXP relSim_randomProfilesC(SEXP listFreqsSEXP, SEXP nBlockSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< List >::type listFreqs(listFreqsSEXP);
    Rcpp::traits::input_parameter< int >::type nBlockSize(nBlockSizeSEXP);
    __result = Rcpp::wrap(randomProfilesC(listFreqs, nBlockSize));
    return __result;
END_RCPP
}
// randomSibsC
IntegerVector randomSibsC(IntegerVector ProfSib1, List listFreqs, int nBlockSize);
RcppExport SEXP relSim_randomSibsC(SEXP ProfSib1SEXP, SEXP listFreqsSEXP, SEXP nBlockSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type ProfSib1(ProfSib1SEXP);
    Rcpp::traits::input_parameter< List >::type listFreqs(listFreqsSEXP);
    Rcpp::traits::input_parameter< int >::type nBlockSize(nBlockSizeSEXP);
    __result = Rcpp::wrap(randomSibsC(ProfSib1, listFreqs, nBlockSize));
    return __result;
END_RCPP
}
// randomChildrenC
IntegerVector randomChildrenC(IntegerVector ProfParent, List listFreqs, int nBlockSize);
RcppExport SEXP relSim_randomChildrenC(SEXP ProfParentSEXP, SEXP listFreqsSEXP, SEXP nBlockSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type ProfParent(ProfParentSEXP);
    Rcpp::traits::input_parameter< List >::type listFreqs(listFreqsSEXP);
    Rcpp::traits::input_parameter< int >::type nBlockSize(nBlockSizeSEXP);
    __result = Rcpp::wrap(randomChildrenC(ProfParent, listFreqs, nBlockSize));
    return __result;
END_RCPP
}
// maximizeLRPC
double maximizeLRPC(IntegerVector ProfParent, IntegerVector ProfChild, List listFreqs, int nBlockSize);
RcppExport SEXP relSim_maximizeLRPC(SEXP ProfParentSEXP, SEXP ProfChildSEXP, SEXP listFreqsSEXP, SEXP nBlockSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type ProfParent(ProfParentSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ProfChild(ProfChildSEXP);
    Rcpp::traits::input_parameter< List >::type listFreqs(listFreqsSEXP);
    Rcpp::traits::input_parameter< int >::type nBlockSize(nBlockSizeSEXP);
    __result = Rcpp::wrap(maximizeLRPC(ProfParent, ProfChild, listFreqs, nBlockSize));
    return __result;
END_RCPP
}
