#include "kentr.h"

//if you forget to close the Bamfile you get a memory leak
BamReader::BamReader(const std::string& bampath){
    const char* cbampath = bampath.c_str();
    in = sam_open(cbampath, "rb");
    if (in == NULL) {
      stop("Failed to open BAM file, check filepath " + bampath);
    }

    bz = in->fp.bgzf ; // bgzf file pointer
    idx = bam_index_load(cbampath); // load BAM index
    header = sam_hdr_read(in) ; // initialize header object

    if (idx == 0) {
      stop("BAM index file is not available for " + bampath);
    }
  }


// [[Rcpp::export]]
DataFrame read_bam(std::string bampath) {
  // open bam
  BamReader bfile(bampath) ;


  // initialize empty alignment container
  bam1_t* read = bam_init1();

  // init vectors for bed data
  std::vector<std::string> chroms ;
  std::vector<std::size_t> starts ;
  std::vector<std::size_t> ends ;
  std::vector<std::string> strands ;

  // iterate through bam, error code is -1
  while (bam_read1(bfile.bz, read) >= 0){
    // extract out chrom, start, end, and strand
    std::string chrom = bfile.header->target_name[(read->core).tid] ;
    std::size_t start = (read->core).pos ;
    std::size_t end = bam_endpos(read) ;
    bool reversed = bam_is_rev(read) ;
    std::string strand = reversed ? "-" : "+" ;

    chroms.push_back(chrom) ;
    starts.push_back(start) ;
    ends.push_back(end) ;
    strands.push_back(strand) ;
  }

  return DataFrame::create(_["chrom"] = chroms,
                           _["start"]= starts,
                           _["end"] = ends,
                           _["strand"] = strands) ;


  bam_destroy1(read);
  bam_hdr_destroy(bfile.header);

}
