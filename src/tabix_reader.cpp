#include "kentr.h"

std::vector<std::string> tsv_values(std::string tsv) {

  std::vector<std::string> values ;
  std::stringstream ss(tsv) ;

  while (ss.good()) {
    std::string substr ;
    getline(ss, substr, '\t') ;

    if (substr.empty()) break ;

    values.push_back(substr) ;
  }

  return values ;
}

TabixReader::TabixReader(const std::string& tbxpath,
                        bool check_idx){
  const char* cpath = tbxpath.c_str();

  in = hts_open(cpath, "rz");
  if (in == NULL) {
    stop("Failed to open tabix file, check filepath " + tbxpath);
  }

  bz = in->fp.bgzf ; // bgzf file pointer
  idx = tbx_index_load(cpath); // load tbx index

  if (idx == 0 && check_idx) {
    stop("tabix index file is not available for " + tbxpath);
  }
}

// [[Rcpp::export]]
List read_tabix(std::string tbxpath,
                     std::string region = "."){
  // open tabix
  TabixReader tfile(tbxpath, true) ;

  hts_itr_t *iter = NULL;
  iter = tbx_itr_querys(tfile.idx, region.c_str());
  // initialize empty data record
  kstring_t str = {0,0,0} ;
  int n_fields ;

  std::vector<std::vector<std::string>> output;

  // determine # of output columns
  if (tbx_itr_next(tfile.in, tfile.idx, iter, &str) >= 0){
    std::vector<std::string> fields = tsv_values(str.s) ;
    n_fields = fields.size() ;
    output.resize(n_fields) ;

    for(int i = 0; i < n_fields; i++){
      output[i].push_back(fields[i]) ;
    }

  } else {
    // return empty df
    std::vector<std::string> c;
    std::vector<double> s;
    std::vector<double> e;

    return Rcpp::DataFrame::create(_("chrom") = c ,
                                   _("start") = s,
                                   _("end") = e) ;

  }

  while(tbx_itr_next(tfile.in, tfile.idx, iter, &str) >= 0) {
   std::vector<std::string> fields = tsv_values(str.s) ;
   if(fields.size() != n_fields){
      Rcpp::Rcout << str.s << std::endl ;
      Rcpp::Rcout << n_fields << " " << fields.size() << " "  << std::endl ;
      for(auto i:fields) {
        Rcpp::Rcout << i << std::endl ;
      }
      Rcpp::stop("mismatch field counts") ;
   }

    for(int i = 0; i < n_fields; i++){
      output[i].push_back(fields[i]) ;
    }

  }

  // figure out output names for chrom start end and possibly strand
  int n_no_names = 0 ;
  int n_strand = 0 ;
  std::vector<std::string> strand_vals = {"+", "-"};
  std::vector<std::string> names(n_fields);
  for(int i = 0; i < n_fields; i++){
    if(i == (tfile.idx->conf.sc - 1)){
      names[i] = "chrom" ;
    } else if (i == (tfile.idx->conf.bc - 1)){
      names[i] = "start";
    } else if (i == (tfile.idx->conf.ec - 1)) {
      names[i] = "end" ;
    } else {
      // check first element, if "+" or "-" assign column as strand
      if(std::count(strand_vals.begin(), strand_vals.end(), output[i][0])) {
        if(n_strand == 0) {
          names[i] = "strand" ;
          continue ;
        }
        n_strand += 1 ;
      }

      n_no_names += 1 ;
      names[i] = "X" + std::to_string(n_no_names) ;
    }
  }

  List res_lst(n_fields);
  for(int i = 0; i < n_fields; i++){
    res_lst[i] = output[i] ;
  }

  DataFrame res = DataFrame::create(res_lst,
                                    _("stringsAsFactors") = false);
  res.attr("names") = names ;

  return res;
}
