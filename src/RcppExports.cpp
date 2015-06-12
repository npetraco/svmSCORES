// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// BootstrapIterKNMScores
NumericVector BootstrapIterKNMScores(NumericMatrix iterProbMat);
RcppExport SEXP svmSCORES_BootstrapIterKNMScores(SEXP iterProbMatSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type iterProbMat(iterProbMatSEXP);
    __result = Rcpp::wrap(BootstrapIterKNMScores(iterProbMat));
    return __result;
END_RCPP
}
// empiricalPvalues
SEXP empiricalPvalues(NumericVector nullvec, NumericVector vec);
RcppExport SEXP svmSCORES_empiricalPvalues(SEXP nullvecSEXP, SEXP vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type nullvec(nullvecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vec(vecSEXP);
    __result = Rcpp::wrap(empiricalPvalues(nullvec, vec));
    return __result;
END_RCPP
}
// compPvals2
SEXP compPvals2(NumericVector nullvec, NumericVector vec);
RcppExport SEXP svmSCORES_compPvals2(SEXP nullvecSEXP, SEXP vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type nullvec(nullvecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vec(vecSEXP);
    __result = Rcpp::wrap(compPvals2(nullvec, vec));
    return __result;
END_RCPP
}
// compPvals3
SEXP compPvals3(NumericVector nullvec, NumericVector vec);
RcppExport SEXP svmSCORES_compPvals3(SEXP nullvecSEXP, SEXP vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type nullvec(nullvecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vec(vecSEXP);
    __result = Rcpp::wrap(compPvals3(nullvec, vec));
    return __result;
END_RCPP
}
