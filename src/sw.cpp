// sw.cpp
//
// Copyright (C) 2018 Kent Riemondy
//
// This file is part of kentr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include "kentr.h"


//' Compute complete striped smith waterman alignment between two dna strings (ATCG)
//' @description Alignments will be computed for both forward and reverse complement. Forward
//' alignments will be indicated with strand "+", whereas reverse complement alignments will
//' be reported as "-". The best alignment is reported
//' @param query_seq query sequence as character vector of length one
//' @param ref_seqs reference sequences as a character vector
//' @export
// [[Rcpp::export(rng = false)]]
DataFrame get_sw(std::string query_seq,
                 std::vector<std::string> ref_seqs) {

  std::transform(query_seq.begin(),
                 query_seq.end(),
                 query_seq.begin(),
                 ::toupper) ;

  int nx = ref_seqs.size();
  auto query_seq_rev = revComp(query_seq) ;

  NumericVector sw_score(nx) ;
  NumericVector sec_sw_score(nx) ;
  IntegerVector rstart(nx) ;
  IntegerVector rend(nx) ;
  IntegerVector qstart(nx) ;
  IntegerVector qend(nx) ;
  IntegerVector sec_rend(nx) ;
  IntegerVector nmis(nx) ;
  CharacterVector cigar(nx) ;
  CharacterVector strand(nx) ;

  for(int i = 0; i < nx; ++i){
    auto ref_seq = ref_seqs[i] ;
    std::transform(ref_seq.begin(),
                   ref_seq.end(),
                   ref_seq.begin(),
                   ::toupper) ;

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

    StripedSmithWaterman::Alignment rv_alignment;

    // Aligns the revcomp query to the ref
    aligner.Align(query_seq_rev.c_str(),
                  ref_seq.c_str(),
                  ref_seq.size(),
                  filter,
                  &rv_alignment,
                  maskLen);
   if (alignment.sw_score >=  rv_alignment.sw_score) {

     sw_score[i] = alignment.sw_score ;
     sec_sw_score[i] = alignment.sw_score_next_best ;
     rstart[i] = alignment.ref_begin;
     rend[i] = alignment.ref_end;
     qstart[i] = alignment.query_begin;
     qend[i] = alignment.query_end;
     sec_rend[i]  = alignment.ref_end_next_best;
     nmis[i] = alignment.mismatches ;
     cigar[i] = alignment.cigar_string ;
     strand[i] = "+" ;
   } else {
     sw_score[i] = rv_alignment.sw_score ;
     sec_sw_score[i] = rv_alignment.sw_score_next_best ;
     rstart[i] = rv_alignment.ref_begin;
     rend[i] = rv_alignment.ref_end;
     qstart[i] = rv_alignment.query_begin;
     qend[i] = rv_alignment.query_end;
     sec_rend[i]  = rv_alignment.ref_end_next_best;
     nmis[i] = rv_alignment.mismatches ;
     cigar[i] = rv_alignment.cigar_string ;
     strand[i] = "-" ;
   }
  }
  return DataFrame::create(_("sw_score") = sw_score,
                          _("secondary_sw_score") = sec_sw_score,
                          _("rstart") = rstart,
                          _("rend") = rend,
                          _("qstart") = qstart,
                          _("qend") = qend,
                          _("secondary_rend") = sec_rend,
                          _("nmismatches") = nmis,
                          _("cigar") = cigar,
                          _("strand") = strand,
                          _("stringsAsFactors") = false) ;
}

/*** R

query <- c("ATCGATCG")
ref <- c("ATCGATCATGCTGATATGG",
         "AATCGATTGATGATGCTGATGCATGCTGATATGG")

get_sw(query, ref)

*/

