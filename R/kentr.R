#' @importFrom Rcpp sourceCpp
#' @importFrom  purrr map_chr map
#' @import dplyr
#' @importFrom stringr str_split str_detect str_replace str_c
#' @import readr
#' @useDynLib kentr, .registration = TRUE
NULL

utils::globalVariables(c(".",
                         "header",
                         "seqnames",
                         "start",
                         "type",
                         "chrom",
                         "end",
                         "transcript_id",
                         "exon_number",
                         "strand"))
