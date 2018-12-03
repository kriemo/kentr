#' Mann whitney test
#' @param df dataframe
#' @param id_col column in df that contains the sample ids. Test only works for two
#' sample comparison
#' @param alternative one of "two.sided", "greater", or "less" see [stats::wilcox.test]
#' @param correct logical see [stats::wilcox.test]
#' @param mu logical see [stats::wilcox.test]
#' @examples
#' ids <- rep(c("A", "B"), each = 1000)
#' nums1 <- rnorm(1000)
#' set.seed(42)
#' offset <- sample(c(-100, 100), 1)
#' nums2 <- nums1 + offset
#' df <- data.frame(ids = ids,
#'                  num = c(nums1, nums2))
#'
#' df$num2 <- nums1 + offset + 1
#'
#'mann_whitney(df, "ids")
#'@export
mann_whitney <- function(df,
                         id_col,
                         alternative = "two.sided",
                         correct = TRUE,
                         mu = 0){

  if(!is.character(id_col)){
    stop("id_col must be a character")
  }

  if(!is.data.frame(df)){
    stop("df must be a data.frame")
  }

  if(!(is.numeric(mu) && length(mu) == 1)){
    stop("mu must be a numeric of length 1")
  }

  if(!(alternative %in% c("greater", "less", "two.sided"))){
    stop("alternative must be one of greater, less or two.sided")
  }

  if(!(is.logical(correct) && length(correct) == 1)){
    stop("correct must be a logical of length 1")
  }

  samples <- unique(df[[id_col]])

  if(length(samples) != 2) {
    stop("only two sample test is implemented")
  }

  sample_1_idx <- which(df[[id_col]] == samples[1])
  sample_2_idx <- which(df[[id_col]] == samples[2])
  df <- df[, colnames(df) != id_col, drop = F]

  if(!all(sapply(df, is.numeric))){
    stop("non-numeric column present in non-id columns")
  }

  # convert to 0 based indexes for cpp
  sample_1_idx <- sample_1_idx - 1
  sample_2_idx <- sample_2_idx - 1

  res <- mw_test_impl(df,
                   sample_1_idx,
                   sample_2_idx,
                   alternative,
                   correct,
                   mu)
  res
}
