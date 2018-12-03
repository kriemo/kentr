// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/kentr.h"
#include <Rcpp.h>

using namespace Rcpp;

// read_bam_tags
DataFrame read_bam_tags(std::string bampath, std::vector<std::string> tag_ids, std::vector<std::string> tag_types, std::string region);
RcppExport SEXP _kentr_read_bam_tags(SEXP bampathSEXP, SEXP tag_idsSEXP, SEXP tag_typesSEXP, SEXP regionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type bampath(bampathSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type tag_ids(tag_idsSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type tag_types(tag_typesSEXP);
    Rcpp::traits::input_parameter< std::string >::type region(regionSEXP);
    rcpp_result_gen = Rcpp::wrap(read_bam_tags(bampath, tag_ids, tag_types, region));
    return rcpp_result_gen;
END_RCPP
}
// read_bam
DataFrame read_bam(std::string bampath, std::string region, std::vector<std::string> tag_ids, std::vector<std::string> tag_types);
RcppExport SEXP _kentr_read_bam(SEXP bampathSEXP, SEXP regionSEXP, SEXP tag_idsSEXP, SEXP tag_typesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type bampath(bampathSEXP);
    Rcpp::traits::input_parameter< std::string >::type region(regionSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type tag_ids(tag_idsSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type tag_types(tag_typesSEXP);
    rcpp_result_gen = Rcpp::wrap(read_bam(bampath, region, tag_ids, tag_types));
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
// revComp
CharacterVector revComp(CharacterVector vec);
RcppExport SEXP _kentr_revComp(SEXP vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type vec(vecSEXP);
    rcpp_result_gen = Rcpp::wrap(revComp(vec));
    return rcpp_result_gen;
END_RCPP
}
// get_hamming_pairs
IntegerVector get_hamming_pairs(std::vector< std::string > bcs_to_test, std::vector< std::string > all_bcs);
RcppExport SEXP _kentr_get_hamming_pairs(SEXP bcs_to_testSEXP, SEXP all_bcsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector< std::string > >::type bcs_to_test(bcs_to_testSEXP);
    Rcpp::traits::input_parameter< std::vector< std::string > >::type all_bcs(all_bcsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_hamming_pairs(bcs_to_test, all_bcs));
    return rcpp_result_gen;
END_RCPP
}
// get_hamming
IntegerVector get_hamming(std::string bc_to_test, std::vector< std::string > all_bcs);
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
// get_kmers
List get_kmers(CharacterVector seqs, int n);
RcppExport SEXP _kentr_get_kmers(SEXP seqsSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type seqs(seqsSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(get_kmers(seqs, n));
    return rcpp_result_gen;
END_RCPP
}
// mw_test_impl
DataFrame mw_test_impl(DataFrame df, IntegerVector idx1, IntegerVector idx2, std::string alternative, bool correct, int mu);
RcppExport SEXP _kentr_mw_test_impl(SEXP dfSEXP, SEXP idx1SEXP, SEXP idx2SEXP, SEXP alternativeSEXP, SEXP correctSEXP, SEXP muSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type df(dfSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type idx1(idx1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type idx2(idx2SEXP);
    Rcpp::traits::input_parameter< std::string >::type alternative(alternativeSEXP);
    Rcpp::traits::input_parameter< bool >::type correct(correctSEXP);
    Rcpp::traits::input_parameter< int >::type mu(muSEXP);
    rcpp_result_gen = Rcpp::wrap(mw_test_impl(df, idx1, idx2, alternative, correct, mu));
    return rcpp_result_gen;
END_RCPP
}
// get_sw
DataFrame get_sw(std::string query_seq, std::vector<std::string> ref_seqs);
RcppExport SEXP _kentr_get_sw(SEXP query_seqSEXP, SEXP ref_seqsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type query_seq(query_seqSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type ref_seqs(ref_seqsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_sw(query_seq, ref_seqs));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_kentr_read_bam_tags", (DL_FUNC) &_kentr_read_bam_tags, 4},
    {"_kentr_read_bam", (DL_FUNC) &_kentr_read_bam, 4},
    {"_kentr_getSeq", (DL_FUNC) &_kentr_getSeq, 2},
    {"_kentr_revComp", (DL_FUNC) &_kentr_revComp, 1},
    {"_kentr_get_hamming_pairs", (DL_FUNC) &_kentr_get_hamming_pairs, 2},
    {"_kentr_get_hamming", (DL_FUNC) &_kentr_get_hamming, 2},
    {"_kentr_find_mismatch", (DL_FUNC) &_kentr_find_mismatch, 2},
    {"_kentr_get_kmers", (DL_FUNC) &_kentr_get_kmers, 2},
    {"_kentr_mw_test_impl", (DL_FUNC) &_kentr_mw_test_impl, 6},
    {"_kentr_get_sw", (DL_FUNC) &_kentr_get_sw, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_kentr(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
