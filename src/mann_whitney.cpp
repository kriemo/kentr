#include "kentr.h"

// rcpp implementation of wilcox.test, with only normal approximation
MannWhitney mw_impl(NumericVector& vec1,
                    NumericVector& vec2,
                    int mu = 0,
                    std::string method = "average",
                    bool correct = true,
                    std::string alternative = "two.sided") {
  MannWhitney results;

 // std::vector<double> vec_std = as<std::vector<double> >(vec) ;
  std::vector<double> vec1_std = as<std::vector<double> >(vec1) ;
  std::vector<double> vec2_std = as<std::vector<double> >(vec2) ;

  if(mu != 0) {
    for(double &val:vec1_std) val -= mu ;
  }

  std::vector<double> ranks ;
  std::vector<double> vec_std(vec1_std);
  vec_std.insert(vec_std.end(), vec2_std.begin(), vec2_std.end()) ;
  rank(vec_std, ranks, method);

  int n_x = vec1_std.size() ;
  int n_y = vec2_std.size() ;
  bool exact_test = n_x < 50 && n_y < 50 ;

  double x_rank_sum = 0.0 ;
  for(int i = 0; i < vec1_std.size(); i++){
    x_rank_sum += ranks[i] ;
  }

  double statistic = x_rank_sum - n_x * (n_x + 1.0)  / 2.0 ;
  double PVAL ;

  std::vector<double> sortedranks(ranks) ;
  std::sort(sortedranks.begin(), sortedranks.end()) ;
  int unique_count = std::unique(sortedranks.begin(),
                                 sortedranks.end()) - sortedranks.begin() ;
  bool tied_ranks = ranks.size() != unique_count ;

  if(exact_test && !tied_ranks) {
    if(alternative == "two.sided") {
      double p ;
      if(statistic > (n_x * n_y / 2.0)) {
          p = R::pwilcox(statistic - 1.0, n_y, n_y, 0, 0) ;
        } else {
          p = R::pwilcox(statistic, n_x, n_y, 1, 0) ;
        }
        PVAL = std::min(2.0 * p, 1.0) ;
    } else if (alternative == "greater") {
        PVAL = R::pwilcox(statistic - 1.0, n_x, n_y, 0, 0) ;
    } else if (alternative == "less") {
        PVAL =  R::pwilcox(statistic, n_x, n_y, 1, 0) ;
    } else {
        stop("unknown alternative parameter") ;
    }
  } else {

    double z = statistic - n_x * n_y / 2.0 ;
    NumericVector ranks_rcpp = wrap(ranks) ;
    NumericVector NTIES = as<NumericVector>(Rcpp::table(ranks_rcpp)) ;

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

     if(alternative == "two.sided") {
       PVAL = 2 * std::min(R::pnorm(z, 0, 1, 1, 0),
                           R::pnorm(z, 0, 1, 0, 0)) ;
     } else if (alternative == "greater") {
       PVAL =  R::pnorm(z, 0, 1, 0, 0) ;
     } else if (alternative == "less") {
       PVAL =  R::pnorm(z, 0, 1, 1, 0) ;
     } else {
       stop("unknown alternative parameter") ;
     }
     if(exact_test && tied_ranks) {
       warning("cannot compute exact p-value with ties") ;
     }
  }
  results.pval = PVAL ;
  results.w = statistic ;
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
    NumericVector vec1 = vec[idx1] ;
    NumericVector vec2 = vec[idx2] ;

    vec1 = vec1[is_finite(vec1)] ;
    vec2 = vec2[is_finite(vec2)] ;
    MannWhitney stats ;

    if(vec1.size() < 1) {
      warning("not enough x values") ;
      stats.pval = NA_REAL ;
      stats.w = NA_INTEGER ;
      mw_stats.push_back(stats) ;
      continue ;
    }

    stats = mw_impl(vec1, vec2) ;
    mw_stats.push_back(stats) ;
  }

  size_t nstats = mw_stats.size();
  NumericVector pvals(nstats) ;
  IntegerVector w_stat(nstats);

  for(int i = 0; i < nstats; i++){
    pvals[i] = mw_stats[i].pval ;
    w_stat[i] = mw_stats[i].w ;
  }

  return DataFrame::create(_["id"] = df_names,
                           _["pval"] = pvals,
                           _["w_stat"] = w_stat,
                           _("stringsAsFactors") = false);
}
