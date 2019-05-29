// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// AIPW_coxscore_cpp
Rcpp::NumericMatrix AIPW_coxscore_cpp(NumericVector time, NumericVector status, NumericMatrix covar, IntegerVector eventid, IntegerVector id, NumericVector score, IntegerVector strata, IntegerMatrix marker, IntegerVector R, NumericMatrix pR, IntegerMatrix total_R, List marker_r, IntegerVector whereX, IntegerVector whereW, NumericVector gamma, NumericMatrix comb_y, int nvar, int n_marker, int nR, int ngamma, int nalp, bool second_cont_bl, bool second_cont_rr);
RcppExport SEXP _AIPWsubtype_AIPW_coxscore_cpp(SEXP timeSEXP, SEXP statusSEXP, SEXP covarSEXP, SEXP eventidSEXP, SEXP idSEXP, SEXP scoreSEXP, SEXP strataSEXP, SEXP markerSEXP, SEXP RSEXP, SEXP pRSEXP, SEXP total_RSEXP, SEXP marker_rSEXP, SEXP whereXSEXP, SEXP whereWSEXP, SEXP gammaSEXP, SEXP comb_ySEXP, SEXP nvarSEXP, SEXP n_markerSEXP, SEXP nRSEXP, SEXP ngammaSEXP, SEXP nalpSEXP, SEXP second_cont_blSEXP, SEXP second_cont_rrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type time(timeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type status(statusSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type covar(covarSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type eventid(eventidSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type id(idSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type score(scoreSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type strata(strataSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type marker(markerSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type R(RSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type pR(pRSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type total_R(total_RSEXP);
    Rcpp::traits::input_parameter< List >::type marker_r(marker_rSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type whereX(whereXSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type whereW(whereWSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type comb_y(comb_ySEXP);
    Rcpp::traits::input_parameter< int >::type nvar(nvarSEXP);
    Rcpp::traits::input_parameter< int >::type n_marker(n_markerSEXP);
    Rcpp::traits::input_parameter< int >::type nR(nRSEXP);
    Rcpp::traits::input_parameter< int >::type ngamma(ngammaSEXP);
    Rcpp::traits::input_parameter< int >::type nalp(nalpSEXP);
    Rcpp::traits::input_parameter< bool >::type second_cont_bl(second_cont_blSEXP);
    Rcpp::traits::input_parameter< bool >::type second_cont_rr(second_cont_rrSEXP);
    rcpp_result_gen = Rcpp::wrap(AIPW_coxscore_cpp(time, status, covar, eventid, id, score, strata, marker, R, pR, total_R, marker_r, whereX, whereW, gamma, comb_y, nvar, n_marker, nR, ngamma, nalp, second_cont_bl, second_cont_rr));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_AIPWsubtype_AIPW_coxscore_cpp", (DL_FUNC) &_AIPWsubtype_AIPW_coxscore_cpp, 23},
    {NULL, NULL, 0}
};

RcppExport void R_init_AIPWsubtype(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
