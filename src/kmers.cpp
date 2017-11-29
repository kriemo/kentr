#include "kentr.h"
#include <string>
#include <algorithm>
#include <unordered_map>

// generate overlapping kmers from a string
CharacterVector getKmers(std::string seqs, int k = 2){

  // make all uppercase
  std::transform(seqs.begin(), seqs.end(), seqs.begin(), toupper) ;

  int lim = seqs.length() - k + 1;
  CharacterVector result( lim );
  for ( int j = 0; j < lim; j++ )
  {
    result[j] = seqs.substr( j, k );
  }
  return result;
}

// count kmers
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
                           _["counts"] = counts,
                           _("stringsAsFactors") = false) ;
}


//' get kmer-counts for character vector of sequences
//' @param n kmer size
//' @export
// [[Rcpp::export]]
List get_kmers(CharacterVector seqs, int n = 2){

  int n_seqs = seqs.size() ;
  List res(n_seqs) ;

  for(int i = 0; i < n_seqs; i++){
    auto seq = as<std::string>(seqs[i]) ;
    auto kmers = getKmers(seq, n) ;
    auto kmer_count = countKmers(kmers) ;
    res[i] = kmer_count ;
  }

  return res ;
}

/*** R
get_kmers(c("ATCGATGCTGATCGT", "ATGCTAGCTAGCTGATATATATCGATGTAGCTG"), 4)
*/
