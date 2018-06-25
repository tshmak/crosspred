// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/pthresh.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

// multiBed4
arma::mat multiBed4(const std::string fileName, int N, int P, const arma::vec weights, arma::Col<int> pbin, int nbin, const arma::Col<int> col_skip_pos, const arma::Col<int> col_skip, const arma::Col<int> keepbytes, const arma::Col<int> keepoffset);
static SEXP _pthresh_multiBed4_try(SEXP fileNameSEXP, SEXP NSEXP, SEXP PSEXP, SEXP weightsSEXP, SEXP pbinSEXP, SEXP nbinSEXP, SEXP col_skip_posSEXP, SEXP col_skipSEXP, SEXP keepbytesSEXP, SEXP keepoffsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const std::string >::type fileName(fileNameSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type P(PSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< arma::Col<int> >::type pbin(pbinSEXP);
    Rcpp::traits::input_parameter< int >::type nbin(nbinSEXP);
    Rcpp::traits::input_parameter< const arma::Col<int> >::type col_skip_pos(col_skip_posSEXP);
    Rcpp::traits::input_parameter< const arma::Col<int> >::type col_skip(col_skipSEXP);
    Rcpp::traits::input_parameter< const arma::Col<int> >::type keepbytes(keepbytesSEXP);
    Rcpp::traits::input_parameter< const arma::Col<int> >::type keepoffset(keepoffsetSEXP);
    rcpp_result_gen = Rcpp::wrap(multiBed4(fileName, N, P, weights, pbin, nbin, col_skip_pos, col_skip, keepbytes, keepoffset));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _pthresh_multiBed4(SEXP fileNameSEXP, SEXP NSEXP, SEXP PSEXP, SEXP weightsSEXP, SEXP pbinSEXP, SEXP nbinSEXP, SEXP col_skip_posSEXP, SEXP col_skipSEXP, SEXP keepbytesSEXP, SEXP keepoffsetSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_pthresh_multiBed4_try(fileNameSEXP, NSEXP, PSEXP, weightsSEXP, pbinSEXP, nbinSEXP, col_skip_posSEXP, col_skipSEXP, keepbytesSEXP, keepoffsetSEXP));
    }
    Rboolean rcpp_isInterrupt_gen = Rf_inherits(rcpp_result_gen, "interrupted-error");
    if (rcpp_isInterrupt_gen) {
        UNPROTECT(1);
        Rf_onintr();
    }
    Rboolean rcpp_isError_gen = Rf_inherits(rcpp_result_gen, "try-error");
    if (rcpp_isError_gen) {
        SEXP rcpp_msgSEXP_gen = Rf_asChar(rcpp_result_gen);
        UNPROTECT(1);
        Rf_error(CHAR(rcpp_msgSEXP_gen));
    }
    UNPROTECT(1);
    return rcpp_result_gen;
}

// validate (ensure exported C++ functions exist before calling them)
static int _pthresh_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("arma::mat(*multiBed4)(const std::string,int,int,const arma::vec,arma::Col<int>,int,const arma::Col<int>,const arma::Col<int>,const arma::Col<int>,const arma::Col<int>)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _pthresh_RcppExport_registerCCallable() { 
    R_RegisterCCallable("pthresh", "_pthresh_multiBed4", (DL_FUNC)_pthresh_multiBed4_try);
    R_RegisterCCallable("pthresh", "_pthresh_RcppExport_validate", (DL_FUNC)_pthresh_RcppExport_validate);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_pthresh_multiBed4", (DL_FUNC) &_pthresh_multiBed4, 10},
    {"_pthresh_RcppExport_registerCCallable", (DL_FUNC) &_pthresh_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_pthresh(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
