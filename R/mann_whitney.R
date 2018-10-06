#' Mann whitney test
#' @param df dataframe
#' @param id_col column in df that contains the sample ids. Test only works for two
#' sample comparison
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
mann_whitney <- function(df, id_col){

  if(!is.character(id_col)){
    stop("id_col must be a character")
  }

  if(!is.data.frame(df)){
    stop("df must be a data.frame")
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

  mw_test_impl(df, sample_1_idx, sample_2_idx)
}
