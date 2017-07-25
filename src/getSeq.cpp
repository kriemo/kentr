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

  // load fasta index
  faidx_t *fai;
  fai = fai_load(cfapath);
  if (fai == NULL) stop("can't load faidx file " + fapath);

  // build sequence and header vectors for each bed entry
  for(int i = 0; i<nr ; i++) {
    auto chrom = as<std::string>(chroms[i]) ;
    auto start = std::to_string(starts[i] + 1) ; // htslib assumes 1 based if input as region
    auto end = std::to_string(ends[i] + 1) ;

    auto region = std::string(chrom + std::string(":") + start + std::string("-") + end) ;
    // pass region arg as c string
    auto reg = region.c_str() ;
    int seq_len ;
    char *seq = fai_fetch(fai, reg, &seq_len);
    if (seq_len >= 0) { // <0 is error
      std::string seq_out(seq, seq_len) ; // might not be necessary
      outseqs[i] = seq ;
      outseqnames[i] = region ;
    } else {
      stop("sequence not found for " + region) ;
    }
    free(seq); //necessary to prevent memory leak
  }

  fai_destroy(fai) ;

  return DataFrame::create(_("header") = outseqnames,
                           _("seq") = outseqs,
                           _("stringsAsFactors") = false) ;
}

/*** R
df <- tibble::tribble(
  ~chrom, ~start, ~end,
  "JH393279.1", 14559891, 14559892
)
getSeq(df, "/Users/kriemo/Projects/Martin/dbases/ensembl85/Ictidomys_tridecemlineatus.spetri2.dna.toplevel.fa")
*/
