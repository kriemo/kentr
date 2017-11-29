// hamming_dist.cpp
//
// Copyright (C) 2017 Kent Riemondy
//
// This file is part of kentr.
//
// This software may be modified and distributed under the terms
// of the MIT license. See the LICENSE file for details.

#include "kentr.h"

std::vector<int> mismatch_position(const std::string& fs, const std::string& ss){
  std::vector<int> mismatch_positions ;

  if((fs.length() == ss.length())){

    for(size_t i = 0; i < fs.length(); i++){
      if(!(fs[i] == ss[i])){
        mismatch_positions.push_back(1) ;
      } else {
        mismatch_positions.push_back(0) ;
      }
    }
  }

  return mismatch_positions ;
}

int hamming_distance(const std::string& fs, const std::string& ss){
  int hm_distance = 0;

  if((fs.length() == ss.length())){

    for(size_t i = 0; i < fs.length(); i++){
      if(!(fs[i] == ss[i])){
        hm_distance++;
      }
    }
  }

  return hm_distance ;
}

//' Compare a string to a vector of strings and report hamming distance
//' @param bc_to_test String to get hamming distance
//' @param all_bcs character vector to compare
//' @export
// [[Rcpp::export]]
std::vector<int> get_hamming(std::string bc_to_test,
                                           std::vector< std::string > all_bcs) {

  int n = all_bcs.size() ;
  std::vector<int> hdists ;
  hdists.reserve(n) ;
  for (int i = 0; i<n; i++){
    int hdist = hamming_distance(bc_to_test, all_bcs[i]) ;
    hdists.push_back(hdist) ;
  }

  return hdists ;
  }

//' Compare a string to a vector of strings report mismatch position
//' @param bc_to_test String to check
//' @param all_bcs character vector to compare
//' @export
// [[Rcpp::export]]
std::vector<int> find_mismatch(std::string bc_to_test,
                               std::vector< std::string > all_bcs) {
  int n = all_bcs.size() ;
  int bc_length = bc_to_test.size() ;
  std::vector<int> mismatch_pos(bc_length) ;
  for (int i = 0; i<n; i++){
    std::vector<int> positions = mismatch_position(bc_to_test, all_bcs[i]) ;
    for(int j = 0; j < bc_length; j++) {
      mismatch_pos[j] += positions[j] ;
    }
  }

  return mismatch_pos ;
}

/*** R

dna_seq <- c("ATCGATCG")
dna_seqs <- c("ATCGATCG", "ATCGATGG")

get_hamming(dna_seq, dna_seqs)
find_mismatch(dna_seq, dna_seqs)

*/

