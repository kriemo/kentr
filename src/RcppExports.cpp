// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/kentr.h"
#include <Rcpp.h>

using namespace Rcpp;

// read_bam
DataFrame read_bam(std::string bampath);
RcppExport SEXP _kentr_read_bam(SEXP bampathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type bampath(bampathSEXP);
    rcpp_result_gen = Rcpp::wrap(read_bam(bampath));
    return rcpp_result_gen;
END_RCPP
}
// getSeq
DataFrame getSeq(DataFrame df, std::string fapath);
RcppExport SEXP _kentr_getSeq(SEXP dfSEXP, SEXP fapathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type df(dfSEXP);
    Rcpp::traits::input_parameter< std::string >::type fapath(fapathSEXP);
    rcpp_result_gen = Rcpp::wrap(getSeq(df, fapath));
    return rcpp_result_gen;
END_RCPP
}
// get_hamming
std::vector<int> get_hamming(std::string bc_to_test, std::vector< std::string > all_bcs);
RcppExport SEXP _kentr_get_hamming(SEXP bc_to_testSEXP, SEXP all_bcsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type bc_to_test(bc_to_testSEXP);
    Rcpp::traits::input_parameter< std::vector< std::string > >::type all_bcs(all_bcsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_hamming(bc_to_test, all_bcs));
    return rcpp_result_gen;
END_RCPP
}
// find_mismatch
std::vector<int> find_mismatch(std::string bc_to_test, std::vector< std::string > all_bcs);
RcppExport SEXP _kentr_find_mismatch(SEXP bc_to_testSEXP, SEXP all_bcsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type bc_to_test(bc_to_testSEXP);
    Rcpp::traits::input_parameter< std::vector< std::string > >::type all_bcs(all_bcsSEXP);
    rcpp_result_gen = Rcpp::wrap(find_mismatch(bc_to_test, all_bcs));
    return rcpp_result_gen;
END_RCPP
}
