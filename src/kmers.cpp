#include "kentr.h"
#include <string>
#include <algorithm>
#include <unordered_map>

std::string canonicalKmer(std::string kmer){
  std::string rc_kmer = revComp(kmer) ;
  if(kmer < rc_kmer){
    return kmer ;
  } else {
    return rc_kmer ;
  }
}

// generate overlapping kmers from a string
CharacterVector getKmers(std::string seqs,
                         int k = 2,
                         bool report_canonical = false){

  // make all uppercase
  std::transform(seqs.begin(), seqs.end(), seqs.begin(), toupper) ;

  int lim = seqs.length() - k + 1;
  CharacterVector result( lim ) ;

  if(report_canonical) {
    for ( int j = 0; j < lim; j++ )
    {
      auto canonical_kmer = canonicalKmer(seqs.substr( j, k )) ;
      result[j] = canonical_kmer ;
    }
  }
  else {
    for ( int j = 0; j < lim; j++ )
    {
      result[j] = seqs.substr( j, k );
    }
  }
  return result;
}

// count kmers
DataFrame countKmers(CharacterVector kmers){

  std::vector<std::string> kmer_name ;
  std::vector<int> counts ;

  // sort kmers, don't use std::sort, see Rcpp known issues
  kmers.sort() ;

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
//' @param seqs character vector of sequences
//' @param n kmer size
//' @param both_strands if TRUE, each kmer and it's reverse complement
//' will be counted as the same kmer. The kmer reported will be
//' the kmer first in lexicographic ordering. Set this if analyzing
//' strand non-specific sequences i.e. some DNA motifs,
//' or unstranded sequencing reads. Default = FALSE
//'
//' @examples
//' get_kmers(c("ATCGATGCTGATCGT",
//'             "ATGCTAGCTAGCTGATATATATCGATGTAGCTG"),
//'             4)
//'
//' get_kmers(c("TTTTAAAA"),
//'          n = 4)
//'
//' get_kmers(c("TTTTAAAA"),
//'           n = 4,
//'           both_strands = TRUE)
//'
//'
//' @export
// [[Rcpp::export]]
List get_kmers(CharacterVector seqs,
               int n = 2,
               bool both_strands = false){

  int n_seqs = seqs.size() ;
  List res(n_seqs) ;

  for(int i = 0; i < n_seqs; i++){
    auto seq = as<std::string>(seqs[i]) ;
    if (seq.size() < n){
      //int zero = 0;
      Rcpp::String chr_s = NA_STRING ;
      Rcpp::String int_s = NA_INTEGER ;
      res[i] = DataFrame::create(_["kmer"] = chr_s,
                                 _["counts"] = int_s,
                                 _("stringsAsFactors") = false) ;
      continue ;
    }
    auto kmers = getKmers(seq, n, both_strands) ;
    auto kmer_count = countKmers(kmers) ;
    res[i] = kmer_count ;
  }

  return res ;
}

/*** R
get_kmers(c("ATCGATGCTGATCGT",
            "ATGCTAGCTAGCTGATATATATCGATGTAGCTG"), 4)

get_kmers(c("TTTTAAAA"),
          n = 4)

get_kmers(c("TTTTAAAA"),
          n = 4,
          both_strands = TRUE)
*/
