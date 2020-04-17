#' Read in tabix indexed file as a data.frame
#' @param filename path to tabix file
#' @param region samtools region query string (i.e. chr1:100-1000)
#'
#' tbx_file <- system.file("extdata", "small_sorted.bam", package = "kentr")
#' res <- tabix_to_df(tbx_file)
#' head(res)
#'
#' res <- tabix_to_df(tbx_file,
#'                    region = "chr1:3e6-3.1e6")
#' res
#' @export
tabix_to_df <- function(filename = NULL,
                        region = "."){
  if(is.null(filename)) filename <- system.file("extdata", "tabix.bed.gz", package = "kentr")

  filename <- path.expand(filename)
  # returned as a list to avoid stringsAsFactors
  df <- read_tabix(filename, region)
  numeric_cols <- intersect(c("start", "end"), colnames(df))
  mutate_at(df, numeric_cols, as.numeric)

}
