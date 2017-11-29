#include "kentr.h"
#include <string>
#include <algorithm>
#include <unordered_map>

//' generate overlapping kmers from a string
//' @param seqs string containing sequences to split into kmers
//' @param k kmer length
//' @export
// [[Rcpp::export]]
CharacterVector getKmers(std::string seqs, int k = 2){
  int lim = seqs.length() - k + 1;
  CharacterVector result( lim );
  for ( int j = 0; j < lim; j++ )
  {
    result[j] = seqs.substr( j, k );
  }
  return result;
}

//' count kmers
//' @param kmers character vector of kmers to count
//' @export
// [[Rcpp::export]]
DataFrame countKmers(CharacterVector kmers){

  std::vector<std::string> kmer_name ;
  std::vector<int> counts ;

  // sort kmers
  std::sort(kmers.begin(), kmers.end()) ;

  // initalize groups
  std::string current_kmer("") ;
  int kmer_count(0) ;

  for(auto i:kmers){
    if(current_kmer == "") {
      current_kmer = i ;
      kmer_count += 1 ;
      continue ;
    }
    if(i == current_kmer) {
      kmer_count += 1 ;
    } else {
      kmer_name.push_back(current_kmer) ;
      counts.push_back(kmer_count) ;
      current_kmer = i ;
      kmer_count = 1 ;
    }
  }
  // add final kmer to list
  kmer_name.push_back(current_kmer) ;
  counts.push_back(kmer_count) ;

  return DataFrame::create(_["kmer"] = kmer_name,
                           _["counts"] = counts) ;
}
/*** R
getKmers("ATGCTAGCTAGCTGATATATATCGATGTAGCTG", 4)
*/
