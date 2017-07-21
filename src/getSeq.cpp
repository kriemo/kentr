#include "kentr.h"
#include <string>

//' Fetch DNA sequence from an indexed fasta file
//' @param df dataframe wtih columns chrom start and end
//' @param fapath path to indexed fasta file
//' @export
// [[Rcpp::export]]
DataFrame getSeq(DataFrame df, std::string fapath) {
  // convert to c string for htslib function call
  const char* cfapath = fapath.c_str();

  int nr = df.nrow() ;
  CharacterVector outseqs(nr) ;
  CharacterVector outseqnames(nr) ;

  CharacterVector chroms = df["chrom"] ;
  IntegerVector starts = df["start"] ;
  IntegerVector ends = df["end"] ;

  faidx_t *fai;
  fai = fai_load(cfapath);
  if (fai == NULL) stop("can't load faidx file " + fapath);

  for(int i = 0; i<nr ; i++) {
    auto chrom = as<std::string>(chroms[i]) ;
    auto start = std::to_string(starts[i]) ;
    auto end = std::to_string(ends[i]) ;

    auto region = std::string(chrom + std::string(":") + start + std::string("-") + end) ;
    Rcpp::Rcout << region << std::endl;
    auto reg = region.c_str() ;
    int seq_len ;
    char *seq = fai_fetch(fai, reg, &seq_len);
    if (seq_len >= 0) {

      std::string seq_out(seq, seq_len) ; // might not be necessary
      Rcpp::Rcout << seq << " sequence" ;
      outseqs[i] = seq ;
      outseqnames[i] = region ;
    } else {
      stop("sequence not found for " + region) ;
    }
    free(seq);
  }

  fai_destroy(fai) ;

  return DataFrame::create(_("header") = outseqnames,
                           _("seq") = outseqs,
                           _("stringsAsFactors") = false) ;

}

/*** R

getSeq(df, "/Users/kriemo/Projects/shared_dbases/genomes/GRCm38.p4.genome.fa")
*/
