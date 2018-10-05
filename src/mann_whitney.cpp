#include "kentr.h"

// rcpp implementation of wilcox.test, with only normal approximation
MannWhitney mw_impl(NumericVector& vec,
                    IntegerVector& idx1,
                    IntegerVector& idx2,
                    int mu = 0,
                    std::string method = "average",
                    bool correct = true,
                    std::string alternative = "two.sided") {
  MannWhitney results;

  std::vector<double> vec_std = as<std::vector<double> >(vec) ;
  std::vector<double> vec1_std = as<std::vector<double> >(vec[idx1]) ;
  std::vector<double> vec2_std = as<std::vector<double> >(vec[idx2]) ;

  if(mu != 0) {
    for(double &val:vec1_std) val -= mu ;
  }

  std::vector<double> ranks ;
  rank(vec_std, ranks, method);

  int n_x = vec1_std.size() ;
  int n_y = vec2_std.size() ;

  int a = 0 ;
  for(int i = 0; i < vec1_std.size(); i++){
    a += ranks[i] ;
  }
  double w = a - n_x * (n_x + 1.0)  / 2.0 ;

  NumericVector ranks_rcpp = wrap(ranks) ;
  NumericVector NTIES = as<NumericVector>(Rcpp::table(ranks_rcpp)) ;

  double z = w - n_x * n_y / 2.0 ;

  double sigma = std::sqrt((n_x * n_y / 12.0) * ((n_x + n_y + 1.0) - sum( Rcpp::pow(NTIES, 3.0) - NTIES)
           / ((n_x + n_y) * (n_x + n_y - 1.0)))) ;

  double CORRECTION = 0.0 ;
  if(correct) {
    if(alternative == "two.sided") {
      CORRECTION  = ((z > 0) ? 1 : ((z < 0) ? -1 : 0)) * 0.5 ;
    } else if (alternative == "greater") {
      CORRECTION = 0.5 ;
    } else if (alternative == "less") {
      CORRECTION = -0.5 ;
    } else {
      stop("unknown alternative parameter") ;
    }
   }
   z = (z - CORRECTION) / sigma ;
   double PVAL ;

   if(alternative == "two.sided") {
    PVAL  = 2 * std::min(R::pnorm(z, 0, 1, 1, 0),
                         R::pnorm(z, 0, 1, 0, 0)) ;
   } else if (alternative == "greater") {
     PVAL =  R::pnorm(z, 0, 1, 0, 0) ;
   } else if (alternative == "less") {
     PVAL =  R::pnorm(z, 0, 1, 1, 0) ;
   } else {
     stop("unknown alternative parameter") ;
   }

  results.pval = PVAL ;
  results.stat = z ;
  results.w = w ;
  return(results) ;
}

// Apply Mann-Whitney test per column of a data.frame
// [[Rcpp::export]]
DataFrame mw_test_impl(DataFrame& df,
                       IntegerVector& idx1,
                       IntegerVector& idx2) {

  size_t ncols = df.ncol() ;
  CharacterVector df_names = df.names() ;;

  std::vector<MannWhitney> mw_stats ;
  for (int i = 0; i < ncols; i++) {

    String current_col = df_names[i] ;
    NumericVector vec = df[current_col] ;

    MannWhitney stats ;
    stats = mw_impl(vec, idx1, idx2) ;
    mw_stats.push_back(stats) ;
  }


  size_t nstats = mw_stats.size();
  NumericVector pvals(nstats) ;
  NumericVector w_stat(nstats);
  NumericVector stat(nstats);

  for(int i = 0; i < nstats; i++){
    pvals[i] = mw_stats[i].pval ;
    w_stat[i] = mw_stats[i].w ;
    stat[i] = mw_stats[i].stat ;
  }

  return DataFrame::create(_["id"] = df_names,
                           _["pval"] = pvals,
                           _["w_stat"] = w_stat,
                           _["test_stat"] = stat,
                           _("stringsAsFactors") = false);
}
