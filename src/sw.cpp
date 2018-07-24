// sw.cpp
//
// Copyright (C) 2018 Kent Riemondy
//
// This file is part of kentr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include "kentr.h"

DataFrame FormatAlignment(const StripedSmithWaterman::Alignment& alignment){

  return DataFrame::create(_("sw_score") = alignment.sw_score,
                    _("secondary_sw_score") = alignment.sw_score_next_best,
                    _("rstart") = alignment.ref_begin,
                    _("rend") = alignment.ref_end,
                    _("qstart") = alignment.query_begin,
                    _("qend") = alignment.query_end,
                    _("secondary_rend") = alignment.ref_end_next_best,
                    _("nmismatches") = alignment.mismatches,
                    _("cigar") = alignment.cigar_string,
                    _("stringsAsFactors") = false) ;
}

//' Compute complete striped smith waterman alignment between two dna strings (ATCG)
//' @param query_seq query
//' @param ref_seq ref
//' @export
// [[Rcpp::export]]
DataFrame get_sw(std::string query_seq,
                 std::string ref_seq) {

  int32_t maskLen = strlen(query_seq.c_str())/2;
  maskLen = maskLen < 15 ? 15 : maskLen;

  // Declares a default Aligner
  StripedSmithWaterman::Aligner aligner;
  // Declares a default filter
  StripedSmithWaterman::Filter filter;
  // Declares an alignment that stores the result
  StripedSmithWaterman::Alignment alignment;
  // Aligns the query to the ref
  aligner.Align(query_seq.c_str(),
                ref_seq.c_str(),
                ref_seq.size(),
                filter,
                &alignment,
                maskLen);

  return FormatAlignment(alignment) ;
}

/*** R

query <- c("ATCGATCG")
ref <- c("ATCGATCATGCTGATATGG")

get_sw(query, ref)

*/

