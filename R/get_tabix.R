#' Read in tabix indexed file as a data.frame
#' @param filename path to tabix file
#' @param region samtools region query string (i.e. chr1:100-1000)
#'
#' @examples
#' tbx_file <- system.file("extdata", "tabix.bed.gz", package = "kentr")
#' res <- tabix_to_df(tbx_file)
#' head(res)
#'
#' res <- tabix_to_df(tbx_file,
#'                    region = "chr1:3e6-3.1e6")
#' res
#' @export
tabix_to_df <- function(filename,
                        region = "."){
  filename <- path.expand(filename)
  # returned as a list to avoid stringsAsFactors
  df <- read_tabix(filename, region)
  numeric_cols <- intersect(c("start", "end", "pos"),  colnames(df))
  mutate_at(df, numeric_cols, as.numeric)

}

#' List chromosomes in a tabix index
#' @param filename path to indexed tabix file
#'
#' @examples
#' tbx_file <- system.file("extdata", "tabix.bed.gz", package = "kentr")
#' res <- get_tabix_chroms(tbx_file)
#' head(res)
#' @export
get_tabix_chroms <- function(filename){
  list_tabix_chroms(path.expand(filename))
}
