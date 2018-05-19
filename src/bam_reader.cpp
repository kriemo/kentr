#include "kentr.h"

//if you forget to close the Bamfile you get a memory leak
BamReader::BamReader(const std::string& bampath,
                     bool check_idx,
                     int cache_size){
    const char* cbampath = bampath.c_str();
    in = sam_open(cbampath, "rb");
    if (in == NULL) {
      stop("Failed to open BAM file, check filepath " + bampath);
    }

    bz = in->fp.bgzf ; // bgzf file pointer
    idx = bam_index_load(cbampath); // load BAM index
    header = sam_hdr_read(in) ; // initialize header object

    if (idx == 0 && check_idx) {
      stop("BAM index file is not available for " + bampath);
    }
  }

// [[Rcpp::export]]
DataFrame read_bam_tags(std::string bampath,
                        std::vector<std::string> get_tags,
                        std::string region = ".") {
  // open bam
  BamReader bfile(bampath, false) ;

  // initialize empty alignment container
  bam1_t* read = bam_init1();

  // init vectors for bed data
  std::vector<std::string> chroms ;
  std::vector<std::size_t> starts ;
  std::vector<std::size_t> ends ;
  std::vector<std::string> strands ;
  std::vector<std::string> readnames ;

  int ntags ;
  ntags = get_tags.size() ;

  Rcpp::List out_df(5 + ntags) ;
  Rcpp::CharacterVector names(5 + ntags) ;

  std::vector<std::vector<std::string>> tags(ntags) ;

  hts_itr_t *iter = NULL;
  iter = sam_itr_querys(bfile.idx, bfile.header, region.c_str());
  if(bfile.header == NULL || iter == NULL) {
    Rcpp::stop("Unable to iterate to region within BAM.");
  }

  // iterate through bam, error code is -1
  while(sam_itr_next(bfile.in, iter, read) >= 0) {
    // extract out chrom, start, end, and strand
    std::string chrom = bfile.header->target_name[(read->core).tid] ;
    std::size_t start = (read->core).pos ;
    std::size_t end = bam_endpos(read) ;
    bool reversed = bam_is_rev(read) ;
    std::string strand = reversed ? "-" : "+" ;

    for (int i = 0; i < ntags; ++i){
      auto tag = std::string(bam_aux2Z(bam_aux_get(read, get_tags[i].c_str())))  ;
      tags[i].push_back(tag) ;
    }

    char *qname = bam_get_qname(read) ;

    chroms.push_back(chrom) ;
    starts.push_back(start) ;
    ends.push_back(end) ;
    strands.push_back(strand) ;
    readnames.push_back(qname) ;
  }

  out_df[0] = chroms ;
  out_df[1] = starts ;
  out_df[2] = ends ;
  out_df[3] = readnames ;
  out_df[4] = strands ;

  names[0] = "chrom" ;
  names[1] = "start" ;
  names[2] = "end" ;
  names[3] = "name" ;
  names[4] = "strand" ;

  for (int i = 0; i < ntags; ++i){
    int j ;
    j= i + 5 ;
    out_df[j] = tags[i] ;
    names[j] = get_tags[i] ;
  }

  Rcpp::DataFrame result(out_df);
  result.attr("names") = names;

  bam_destroy1(read);
  bam_hdr_destroy(bfile.header);

  return result;
}


DataFrame read_bam_default(std::string bampath,
                          std::string region = ".") {
  // open bam
  BamReader bfile(bampath, true) ;

  // initialize empty alignment container
  bam1_t* read = bam_init1();

  // init vectors for bed data
  std::vector<std::string> chroms ;
  std::vector<std::size_t> starts ;
  std::vector<std::size_t> ends ;
  std::vector<std::string> strands ;
  std::vector<std::string> readnames ;

  hts_itr_t *iter = NULL;
  iter = sam_itr_querys(bfile.idx, bfile.header, region.c_str());
  if(bfile.header == NULL || iter == NULL) {
    Rcpp::stop("Unable to iterate to region within BAM.");
  }

  // iterate through bam, error code is -1
  while(sam_itr_next(bfile.in, iter, read) >= 0) {
    // extract out chrom, start, end, and strand
    std::string chrom = bfile.header->target_name[(read->core).tid] ;
    std::size_t start = (read->core).pos ;
    std::size_t end = bam_endpos(read) ;
    bool reversed = bam_is_rev(read) ;
    std::string strand = reversed ? "-" : "+" ;
    char *qname = bam_get_qname(read) ;

    chroms.push_back(chrom) ;
    starts.push_back(start) ;
    ends.push_back(end) ;
    strands.push_back(strand) ;
    readnames.push_back(qname) ;
  }

  return DataFrame::create(_["chrom"] = chroms,
                           _["start"]= starts,
                           _["end"] = ends,
                           _["name"] = readnames,
                           _["strand"] = strands,
                           _("stringsAsFactors") = false) ;


  bam_destroy1(read);
  bam_hdr_destroy(bfile.header);

}

// [[Rcpp::export]]
DataFrame read_bam(std::string bampath,
                   std::string region,
                   std::vector<std::string> get_tags){
  DataFrame res ;
  if (get_tags[0] != "") {
    res = read_bam_tags(bampath, get_tags, region) ;
  } else {
    res = read_bam_default(bampath, region) ;
  }
  return res ;
}
