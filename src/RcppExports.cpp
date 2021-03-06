// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/crosspred.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

// multiBed4
arma::mat multiBed4(const std::string fileName, int N, int P, const arma::vec weights, arma::Col<int> pbin, int nbin, const arma::Col<int> col_skip_pos, const arma::Col<int> col_skip, const arma::Col<int> keepbytes, const arma::Col<int> keepoffset, const int trace);
static SEXP _crosspred_multiBed4_try(SEXP fileNameSEXP, SEXP NSEXP, SEXP PSEXP, SEXP weightsSEXP, SEXP pbinSEXP, SEXP nbinSEXP, SEXP col_skip_posSEXP, SEXP col_skipSEXP, SEXP keepbytesSEXP, SEXP keepoffsetSEXP, SEXP traceSEXP) {
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
    Rcpp::traits::input_parameter< const int >::type trace(traceSEXP);
    rcpp_result_gen = Rcpp::wrap(multiBed4(fileName, N, P, weights, pbin, nbin, col_skip_pos, col_skip, keepbytes, keepoffset, trace));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _crosspred_multiBed4(SEXP fileNameSEXP, SEXP NSEXP, SEXP PSEXP, SEXP weightsSEXP, SEXP pbinSEXP, SEXP nbinSEXP, SEXP col_skip_posSEXP, SEXP col_skipSEXP, SEXP keepbytesSEXP, SEXP keepoffsetSEXP, SEXP traceSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_crosspred_multiBed4_try(fileNameSEXP, NSEXP, PSEXP, weightsSEXP, pbinSEXP, nbinSEXP, col_skip_posSEXP, col_skipSEXP, keepbytesSEXP, keepoffsetSEXP, traceSEXP));
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
// overallbeta
arma::vec overallbeta(const std::string fileName, int N, int P, arma::Col<int> col_skip_pos, arma::Col<int> col_skip, arma::Col<int> keepbytes, arma::Col<int> keepoffset, arma::vec pred, arma::vec meanbeta, const std::string save, const std::string load);
static SEXP _crosspred_overallbeta_try(SEXP fileNameSEXP, SEXP NSEXP, SEXP PSEXP, SEXP col_skip_posSEXP, SEXP col_skipSEXP, SEXP keepbytesSEXP, SEXP keepoffsetSEXP, SEXP predSEXP, SEXP meanbetaSEXP, SEXP saveSEXP, SEXP loadSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< const std::string >::type fileName(fileNameSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type P(PSEXP);
    Rcpp::traits::input_parameter< arma::Col<int> >::type col_skip_pos(col_skip_posSEXP);
    Rcpp::traits::input_parameter< arma::Col<int> >::type col_skip(col_skipSEXP);
    Rcpp::traits::input_parameter< arma::Col<int> >::type keepbytes(keepbytesSEXP);
    Rcpp::traits::input_parameter< arma::Col<int> >::type keepoffset(keepoffsetSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pred(predSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type meanbeta(meanbetaSEXP);
    Rcpp::traits::input_parameter< const std::string >::type save(saveSEXP);
    Rcpp::traits::input_parameter< const std::string >::type load(loadSEXP);
    rcpp_result_gen = Rcpp::wrap(overallbeta(fileName, N, P, col_skip_pos, col_skip, keepbytes, keepoffset, pred, meanbeta, save, load));
    return rcpp_result_gen;
END_RCPP_RETURN_ERROR
}
RcppExport SEXP _crosspred_overallbeta(SEXP fileNameSEXP, SEXP NSEXP, SEXP PSEXP, SEXP col_skip_posSEXP, SEXP col_skipSEXP, SEXP keepbytesSEXP, SEXP keepoffsetSEXP, SEXP predSEXP, SEXP meanbetaSEXP, SEXP saveSEXP, SEXP loadSEXP) {
    SEXP rcpp_result_gen;
    {
        Rcpp::RNGScope rcpp_rngScope_gen;
        rcpp_result_gen = PROTECT(_crosspred_overallbeta_try(fileNameSEXP, NSEXP, PSEXP, col_skip_posSEXP, col_skipSEXP, keepbytesSEXP, keepoffsetSEXP, predSEXP, meanbetaSEXP, saveSEXP, loadSEXP));
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
static int _crosspred_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
        signatures.insert("arma::mat(*multiBed4)(const std::string,int,int,const arma::vec,arma::Col<int>,int,const arma::Col<int>,const arma::Col<int>,const arma::Col<int>,const arma::Col<int>,const int)");
        signatures.insert("arma::vec(*overallbeta)(const std::string,int,int,arma::Col<int>,arma::Col<int>,arma::Col<int>,arma::Col<int>,arma::vec,arma::vec,const std::string,const std::string)");
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _crosspred_RcppExport_registerCCallable() { 
    R_RegisterCCallable("crosspred", "_crosspred_multiBed4", (DL_FUNC)_crosspred_multiBed4_try);
    R_RegisterCCallable("crosspred", "_crosspred_overallbeta", (DL_FUNC)_crosspred_overallbeta_try);
    R_RegisterCCallable("crosspred", "_crosspred_RcppExport_validate", (DL_FUNC)_crosspred_RcppExport_validate);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_crosspred_multiBed4", (DL_FUNC) &_crosspred_multiBed4, 11},
    {"_crosspred_overallbeta", (DL_FUNC) &_crosspred_overallbeta, 11},
    {"_crosspred_RcppExport_registerCCallable", (DL_FUNC) &_crosspred_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_crosspred(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
